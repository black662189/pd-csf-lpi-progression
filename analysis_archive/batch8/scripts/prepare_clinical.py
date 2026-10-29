import os
"""Prepare batch-8 clinical descriptions and frozen subpart joins; fit no models."""
from pathlib import Path
import hashlib
import json
import platform
import sys

import numpy as np
import pandas as pd


BASE = Path(os.environ['PD_EVIDENCE_ROOT'])
BATCH = Path(__file__).resolve().parents[1]
PREVIOUS = BATCH.parent / "batch6/results/lpi_v1"
OUT = BATCH / "results/descriptive_v1"
KEY = ["PATNO", "EVENT_ID"]
RATER_ITEMS = ["NP1COG", "NP1HALL", "NP1DPRS", "NP1ANXS", "NP1APAT", "NP1DDS"]
PATIENT_ITEMS = ["NP1SLPN", "NP1SLPD", "NP1PAIN", "NP1URIN", "NP1CNST", "NP1LTHD", "NP1FATG"]
ITEMS = RATER_ITEMS + PATIENT_ITEMS
SOURCES = {
    "anchors_CSF": PREVIOUS / "anchors_CSF.tsv",
    "anchors_plasma": PREVIOUS / "anchors_plasma.tsv",
    "paired_eligibility": PREVIOUS / "paired_person_eligibility.tsv",
    "paired_exposures": PREVIOUS / "paired_exposures.tsv",
    "clinical": BASE / "ppmi_clinical_analysis_20261023/cleaned/clinical_model_ready.tsv.gz",
    "curated": BASE / "ppmi_clinical_analysis_20261023/cleaned/curated_original_snapshot.tsv.gz",
    "rater": BASE / "ppmi_intake_20261023/inputs/MDS-UPDRS_Part_I_23Oct2026.csv",
    "patient": BASE / "ppmi_intake_20261023/inputs/MDS-UPDRS_Part_I_Patient_Questionnaire_23Oct2026.csv",
    "dictionary": BASE / "ppmi_clinical_analysis_20261023/reference/curated_dictionary.tsv",
    "plan": BATCH / "protocol/第八批分析计划_v1_20261028.txt",
    "planning_audit": BATCH / "metadata/planning_feasibility_20261028.json",
}
CONTINUOUS = {
    "age_visit": ("Age", "years"),
    "disease_years": ("Time since PD diagnosis", "years"),
    "updrs1_score": ("Baseline MDS-UPDRS Part I", "score, 0-52"),
    "updrs_main": ("Baseline MDS-UPDRS Part III OFF/untreated", "score, 0-132"),
    "ledd": ("Baseline total levodopa equivalent daily dose", "original daily dose value, not log-transformed"),
    "BMI": ("Baseline BMI", "kg/m^2"),
    "future_visits": ("Number of strictly future visits", "visits per person"),
    "first_future_year": ("Time to first future observation", "years since original anchor"),
    "last_future_year": ("Maximum future observation time", "years since original anchor"),
    "future_span_years": ("First-to-last future observation span", "years"),
}
CATEGORICAL = {
    "sex": "Sex (0=Female, 1=Male)",
    "genetic_group": "Genetic grouping",
    "subgroup": "Recorded genetic/enrollment subgroup",
    "hy_main": "Baseline H-Y OFF/untreated",
}
CHECKS = []


def check(name, condition, detail=""):
    if not bool(condition):
        raise AssertionError(f"{name}: {detail}")
    CHECKS.append({"check": name, "pass": True, "detail": detail})


def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype={"PATNO": str}, na_values=["NA"])


def write(frame, name):
    frame.to_csv(OUT / name, sep="\t", index=False, na_rep="NA", float_format="%.15g")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize_continuous(series):
    values = pd.to_numeric(series, errors="raise").dropna()
    return {
        "available": len(values), "missing": int(series.isna().sum()),
        "mean": values.mean(), "sd": values.std(ddof=1),
        "median": values.median(), "q25": values.quantile(.25), "q75": values.quantile(.75),
        "min": values.min(), "max": values.max(),
    }


def describe(groups):
    numeric, categorical = [], []
    for name, data in groups.items():
        for field, (label, unit) in CONTINUOUS.items():
            numeric.append({"cohort": name, "people": len(data), "field": field,
                            "label": label, "unit": unit, **summarize_continuous(data[field])})
        for field, label in CATEGORICAL.items():
            values = data[field].astype("string").fillna("Missing")
            for level, count in values.value_counts(dropna=False).sort_index().items():
                categorical.append({"cohort": name, "people": len(data), "field": field,
                                    "label": label, "level": level, "count": count,
                                    "proportion_all_people": count / len(data)})
            if "Missing" not in values.values:
                categorical.append({"cohort": name, "people": len(data), "field": field,
                                    "label": label, "level": "Missing", "count": 0,
                                    "proportion_all_people": 0.0})
    return pd.DataFrame(numeric), pd.DataFrame(categorical)


def future_for_anchors(clinical, anchors):
    eligible_anchor = anchors.loc[np.isfinite(anchors.x_frozen)
                                 & anchors.updrs1_score.between(0, 52)
                                 & np.isfinite(anchors.updrs_main)]
    future = clinical.loc[(clinical.pd_eligible == 1) & clinical.updrs1_score.between(0, 52),
                          KEY + ["visit_month", "updrs1_score"]]
    future = future.merge(eligible_anchor[["PATNO", "visit_month"]].rename(columns={"visit_month": "index_month"}),
                          on="PATNO", validate="many_to_one")
    future = future.loc[future.visit_month > future.index_month].copy()
    future = future.sort_values(["PATNO", "visit_month", "EVENT_ID"]).drop_duplicates(["PATNO", "visit_month"])
    future["time"] = (future.visit_month - future.index_month) / 12
    check("strict_future_unique_" + str(len(anchors)), (future.time > 0).all()
          and not future.duplicated(["PATNO", "time"]).any())
    return future


def with_followup(anchors, future):
    counts = future.groupby("PATNO").agg(future_visits=("time", "size"),
                                         first_future_year=("time", "min"), last_future_year=("time", "max"))
    data = anchors.merge(counts, on="PATNO", how="left", validate="one_to_one")
    data["future_visits"] = data.future_visits.fillna(0).astype(int)
    data["future_span_years"] = data.last_future_year - data.first_future_year
    return data


def join_forms(data, forms):
    audit = data.merge(forms, on=KEY, how="left", validate="one_to_one")
    for field in ["has_rater", "has_patient"]:
        audit[field] = audit[field].eq(True)
    item_values = audit[ITEMS].apply(pd.to_numeric, errors="raise")
    audit["both_forms"] = audit.has_rater & audit.has_patient
    audit["complete_valid_13"] = item_values.apply(lambda v: v.between(0, 4) & (v % 1 == 0)).all(axis=1)
    audit["rater_sum"] = item_values[RATER_ITEMS].sum(axis=1, min_count=6)
    audit["patient_sum"] = item_values[PATIENT_ITEMS].sum(axis=1, min_count=7)
    audit["subtotals_agree"] = (audit.rater_sum == audit.NP1RTOT) & (audit.patient_sum == audit.NP1PTOT)
    audit["frozen_total_agrees"] = (audit.rater_sum + audit.patient_sum) == audit.updrs1_score
    audit["date_same"] = audit.INFODT_rater.notna() & audit.INFODT_patient.notna() & audit.INFODT_rater.eq(audit.INFODT_patient)
    audit["eligible"] = audit.both_forms & audit.complete_valid_13 & audit.subtotals_agree & audit.frozen_total_agrees & audit.date_same
    audit["exclusion_reason"] = np.select(
        [~audit.both_forms, ~audit.complete_valid_13, ~audit.subtotals_agree,
         ~audit.frozen_total_agrees, ~audit.date_same],
        ["missing_form", "invalid_or_missing_item", "subtotal_disagreement", "frozen_total_disagreement", "form_date_disagreement"],
        default="included")
    return audit


def main():
    if (OUT / "preparation_summary.json").exists():
        raise FileExistsError("Completed output exists. Use a new version instead of overwriting it.")
    OUT.mkdir(parents=True, exist_ok=True)
    source_manifest = pd.DataFrame([{"role": role, "path": str(path), "bytes": path.stat().st_size, "sha256": sha(path)}
                                    for role, path in SOURCES.items()])
    planned = json.loads(SOURCES["planning_audit"].read_text(encoding="utf-8"))
    for source in planned["sources"]:
        role = "anchors_CSF" if source["role"] == "anchors" else source["role"]
        check("planning_source_hash_" + role, source_manifest.set_index("role").loc[role, "sha256"] == source["sha256"])
    clinical = read_tsv(SOURCES["clinical"])
    curated = pd.read_csv(SOURCES["curated"], sep="\t", dtype={"PATNO": str}, usecols=KEY + ["updrs1_score"])
    clinical = clinical.merge(curated, on=KEY, how="left", validate="one_to_one")
    clinical = clinical.loc[(clinical.source_curated == 1) & np.isfinite(clinical.visit_month)].copy()
    check("clinical_key_unique", not clinical.duplicated(KEY).any())
    check("curated_total_range", clinical.updrs1_score.dropna().between(0, 52).all())
    anchors, futures, people = {}, {}, {}
    for assay, expected_n, expected_rows in [("CSF", 200, 1389), ("plasma", 191, 1248)]:
        saved = read_tsv(SOURCES["anchors_" + assay])
        check(assay + "_anchor_person_unique", saved.PATNO.is_unique)
        data = saved.rename(columns={"visit_month": "saved_visit_month", "BMI": "saved_BMI"}).merge(
            clinical, on=KEY, how="left", validate="one_to_one")
        check(assay + "_anchor_date_matches", data.visit_month.eq(data.saved_visit_month).all())
        check(assay + "_anchor_BMI_matches", np.allclose(data.BMI, data.saved_BMI, equal_nan=True, rtol=0, atol=1e-10))
        full_future = future_for_anchors(clinical, data)
        data = with_followup(data, full_future)
        data["reconstructed_in_model"] = np.isfinite(data.x_frozen) & data.updrs1_score.between(0, 52) & np.isfinite(data.updrs_main) & data.future_visits.ge(2)
        check(assay + "_frozen_inclusion_matches", data.in_model.eq(data.reconstructed_in_model).all())
        future = full_future.loc[full_future.PATNO.isin(data.loc[data.in_model, "PATNO"])].copy()
        check(assay + "_frozen_counts", future.PATNO.nunique() == expected_n and len(future) == expected_rows)
        anchors[assay], futures[assay], people[assay] = data, future, data.loc[data.in_model].copy()
        write(future, "frozen_future_" + assay + ".tsv")
        write(data, "anchor_eligibility_" + assay + ".tsv")
    pair_keys = read_tsv(SOURCES["paired_eligibility"])
    pair_x = read_tsv(SOURCES["paired_exposures"])
    pair_ids = set(pair_keys.loc[pair_keys.same_anchor, "PATNO"])
    check("paired_person_keys", len(pair_ids) == 137 and pair_ids == set(pair_x.PATNO) and pair_x.PATNO.is_unique)
    common = futures["CSF"].loc[futures["CSF"].PATNO.isin(pair_ids)].merge(
        futures["plasma"].loc[futures["plasma"].PATNO.isin(pair_ids)], on=KEY + ["time"],
        validate="one_to_one", suffixes=("_CSF", "_plasma"))
    for field in ["visit_month", "index_month", "updrs1_score"]:
        check("paired_" + field + "_matches", common[field + "_CSF"].eq(common[field + "_plasma"]).all())
    paired_future = common[KEY + ["time", "visit_month_CSF", "index_month_CSF", "updrs1_score_CSF"]].rename(
        columns={field + "_CSF": field for field in ["visit_month", "index_month", "updrs1_score"]})
    check("paired_frozen_counts", paired_future.PATNO.nunique() == 137 and len(paired_future) == 923
          and paired_future.groupby("PATNO").size().ge(2).all())
    paired_anchors = people["CSF"].loc[people["CSF"].PATNO.isin(pair_ids)].drop(columns=["future_visits", "first_future_year", "last_future_year", "future_span_years"])
    people["paired"] = with_followup(paired_anchors, paired_future)
    futures["paired"] = paired_future
    write(paired_future, "frozen_future_paired.tsv")

    numeric, categorical = describe(people)
    write(numeric, "table1_continuous.tsv")
    write(categorical, "table1_categorical.tsv")
    wide = []
    for field, (label, unit) in CONTINUOUS.items():
        row = {"field": field, "label": label, "unit": unit}
        for group in people:
            val = numeric.loc[(numeric.cohort == group) & (numeric.field == field)].iloc[0]
            row[group] = f"{val['median']:.2f} [{val['q25']:.2f}, {val['q75']:.2f}]; available {int(val['available'])}/{int(val['people'])}"
        wide.append(row)
    for field, label in CATEGORICAL.items():
        for level in sorted(categorical.loc[categorical.field == field, "level"].unique()):
            row = {"field": field, "label": label + ": " + str(level), "unit": "n (% of all people)"}
            for group, data in people.items():
                found = categorical.loc[(categorical.cohort == group) & (categorical.field == field) & (categorical.level == level)]
                count = int(found.iloc[0]["count"]) if len(found) else 0
                row[group] = f"{count} ({100 * count / len(data):.1f}%)"
            wide.append(row)
    write(pd.DataFrame(wide), "table1_display.tsv")
    combined_people = pd.concat([data.assign(cohort=group) for group, data in people.items()], ignore_index=True)
    write(combined_people, "cohort_person_descriptions.tsv")

    csf = anchors["CSF"].copy()
    check("CSF_227_anchor_count", len(csf) == 227)
    csf["inclusion_reason"] = np.select(
        [~np.isfinite(csf.x_frozen), ~csf.updrs1_score.between(0, 52), ~np.isfinite(csf.updrs_main), csf.future_visits.eq(0), csf.future_visits.eq(1)],
        ["missing_LPI_at_original_anchor", "missing_baseline_Part_I", "missing_baseline_Part_III", "zero_eligible_future_visits", "one_eligible_future_visit"], default="included")
    check("CSF_200_27_flow", csf.in_model.sum() == 200 and (~csf.in_model).sum() == 27
          and csf.in_model.eq(csf.inclusion_reason.eq("included")).all())
    write(csf, "CSF_227_to_200_person_flow.tsv")
    flow = csf.groupby("inclusion_reason").agg(people=("PATNO", "size"), eligible_future_records=("future_visits", "sum")).reset_index()
    write(flow, "CSF_227_to_200_flow_summary.tsv")
    included, excluded = csf.loc[csf.in_model], csf.loc[~csf.in_model]
    ncomp, ccomp = describe({"included_200": included, "not_included_27": excluded})
    write(ncomp, "inclusion_comparison_continuous.tsv")
    write(ccomp, "inclusion_comparison_categorical.tsv")
    smd = []
    for field, (label, unit) in CONTINUOUS.items():
        if field in ["future_visits", "first_future_year", "last_future_year", "future_span_years"]:
            continue
        a, b = included[field].dropna(), excluded[field].dropna()
        denominator = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
        smd.append({"field": field, "level": "continuous", "included_available": len(a), "excluded_available": len(b),
                    "included_missing": int(included[field].isna().sum()), "excluded_missing": int(excluded[field].isna().sum()),
                    "included_mean_or_proportion": a.mean(), "excluded_mean_or_proportion": b.mean(),
                    "smd_included_minus_excluded": (a.mean() - b.mean()) / denominator if denominator > 0 else np.nan,
                    "definition": "mean difference / sqrt((sample variance included + sample variance excluded)/2)"})
    for field in CATEGORICAL:
        a, b = included[field].dropna().astype(str), excluded[field].dropna().astype(str)
        for level in sorted(set(a) | set(b)):
            pa, pb = a.eq(level).mean(), b.eq(level).mean()
            denominator = np.sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
            smd.append({"field": field, "level": level, "included_available": len(a), "excluded_available": len(b),
                        "included_missing": int(included[field].isna().sum()), "excluded_missing": int(excluded[field].isna().sum()),
                        "included_mean_or_proportion": pa, "excluded_mean_or_proportion": pb,
                        "smd_included_minus_excluded": (pa - pb) / denominator if denominator > 0 else np.nan,
                        "definition": "one-level-versus-other proportion difference / sqrt(mean Bernoulli variance); descriptive only"})
    write(pd.DataFrame(smd), "inclusion_SMD.tsv")

    coverage, spans = [], []
    for group, data in futures.items():
        year = np.ceil(data.time).astype(int)
        for end in range(1, int(year.max()) + 1):
            subset = data.loc[year == end]
            n_people = subset.PATNO.nunique()
            coverage.append({"cohort": group, "interval_start_exclusive": end - 1, "interval_end_inclusive": end,
                             "interval": f"({end-1},{end}]", "people_with_observation": n_people, "records": len(subset),
                             "original_model_people": len(people[group]), "fraction_original_people": n_people / len(people[group]),
                             "at_least_25pct_coverage": n_people / len(people[group]) >= .25,
                             "interpretation": "actual observed people and records; not a survival risk set"})
        per = data.groupby("PATNO").agg(future_visits=("time", "size"), first_future_year=("time", "min"), last_future_year=("time", "max")).reset_index()
        per["future_span_years"] = per.last_future_year - per.first_future_year
        per["cohort"] = group
        spans.append(per)
    coverage = pd.DataFrame(coverage)
    write(coverage, "annual_observation_coverage.tsv")
    write(pd.concat(spans, ignore_index=True), "individual_followup_spans.tsv")
    write(coverage.loc[coverage.interval_end_inclusive.isin([1, 3, 5])], "prespecified_1_3_5_year_display_eligibility.tsv")

    rater = pd.read_csv(SOURCES["rater"], dtype={"PATNO": str})
    patient = pd.read_csv(SOURCES["patient"], dtype={"PATNO": str})
    for name, raw in [("rater", rater), ("patient", patient)]:
        check(name + "_unique_visit_key", not raw.duplicated(KEY).any())
        check(name + "_REC_ID_unique_nonmissing", raw.REC_ID.is_unique and raw.REC_ID.notna().all())
    forms = rater[KEY + ["REC_ID", "INFODT", "NUPSOURC", "NP1RTOT"] + RATER_ITEMS].assign(has_rater=True).merge(
        patient[KEY + ["REC_ID", "INFODT", "NUPSOURC", "NP1PTOT"] + PATIENT_ITEMS].assign(has_patient=True),
        on=KEY, how="outer", suffixes=("_rater", "_patient"), validate="one_to_one")
    baseline = join_forms(people["CSF"][KEY + ["visit_month", "updrs1_score"]], forms)
    future_audit = join_forms(futures["CSF"], forms)
    check("subpart_all_200_baselines_qualified", len(baseline) == 200 and baseline.eligible.all())
    for scope, audit in [("baseline", baseline), ("future", future_audit)]:
        check(scope + "_no_valid_total_disagreement", not (audit.complete_valid_13 & ~audit.frozen_total_agrees).any())
        check(scope + "_no_valid_subtotal_disagreement", not (audit.complete_valid_13 & ~audit.subtotals_agree).any())
    included_future = future_audit.loc[future_audit.eligible]
    check("subpart_future_200_1379", len(included_future) == 1379 and included_future.PATNO.nunique() == 200
          and included_future.groupby("PATNO").size().ge(2).all())
    check("subpart_three_missing_forms_seven_date_disagreements", future_audit.exclusion_reason.eq("missing_form").sum() == 3
          and future_audit.exclusion_reason.eq("form_date_disagreement").sum() == 7
          and (~future_audit.eligible).sum() == 10)
    write(baseline, "subpart_baseline_eligibility_ledger.tsv")
    write(future_audit, "subpart_future_eligibility_ledger.tsv")
    write(future_audit.loc[~future_audit.eligible], "subpart_excluded_future_records.tsv")
    write(baseline[KEY + ["NP1RTOT", "NP1PTOT", "updrs1_score"]], "subpart_baseline.tsv")
    write(included_future[KEY + ["NP1RTOT", "NP1PTOT", "updrs1_score"]], "subpart_future.tsv")
    write(future_audit.groupby("exclusion_reason").agg(records=("PATNO", "size"), people=("PATNO", "nunique")).reset_index(),
          "subpart_eligibility_summary.tsv")
    item_description = []
    for scope, data in [("baseline_200", baseline), ("frozen_future_1389", future_audit), ("qualified_future_1379", included_future)]:
        for item in ITEMS:
            values = pd.to_numeric(data[item], errors="raise")
            valid = values.between(0, 4) & values.mod(1).eq(0)
            record = {"scope": scope, "unit_of_description": "visit records", "item": item,
                      "official_form": "rater_6_items" if item in RATER_ITEMS else "patient_7_items",
                      "records": len(data), "people": data.PATNO.nunique(), "missing": int(values.isna().sum()),
                      "invalid_nonmissing": int((values.notna() & ~valid).sum()), "valid_observations": int(valid.sum()),
                      "people_with_valid_observation": data.loc[valid, "PATNO"].nunique(),
                      "zero_fraction_valid_observations": values.loc[valid].eq(0).mean()}
            for score in range(5):
                record[f"score_{score}_count"] = int(values.eq(score).sum())
            item_description.append(record)
    write(pd.DataFrame(item_description), "13_item_distribution_missingness.tsv")
    check("no_CSF_baseline_scale_change", people["CSF"].x_frozen.notna().all())
    check("CSF_genetic_frozen_counts", people["CSF"].genetic_group.value_counts().to_dict() == {"genetic": 126, "sporadic": 74})
    for role, source in SOURCES.items():
        check("input_unchanged_" + role, sha(source) == source_manifest.set_index("role").loc[role, "sha256"])
    write(source_manifest, "input_sha256.tsv")
    write(pd.DataFrame(CHECKS), "preparation_checks.tsv")
    summary = {
        "status": "executed_and_verified_no_models_fitted", "CSF_people": 200, "CSF_future_rows": 1389,
        "plasma_people": 191, "plasma_future_rows": 1248, "paired_people": 137, "paired_future_rows": 923,
        "CSF_anchors": 227, "CSF_not_in_model": 27,
        "CSF_exclusion_reasons": csf.loc[~csf.in_model, "inclusion_reason"].value_counts().to_dict(),
        "subpart_baseline_people": len(baseline), "subpart_future_people": included_future.PATNO.nunique(),
        "subpart_future_rows": len(included_future), "subpart_future_exclusion_reasons": future_audit.loc[~future_audit.eligible, "exclusion_reason"].value_counts().to_dict(),
        "checks_passed": len(CHECKS), "models_fitted": 0,
        "python": sys.version, "pandas": pd.__version__, "numpy": np.__version__, "platform": platform.platform(),
        "script_sha256": sha(Path(__file__)),
        "notes": ["No raw input changed; no downloads, SSH or model fitting in this module.",
                  "LEDD is the original nonnegative visit total daily equivalent dose; descriptive table does not log-transform it.",
                  "SMD uses observed values and describes inclusion; no P values, weighting or imputation.",
                  "Annual bins are actual observation counts, not survival risk sets.",
                  "Two form dates compared at the supplied precision; no claim of simultaneous sampling.",
                  "All output files are derived from frozen anchors; no later anchor substitution."]
    }
    (OUT / "preparation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    note = f"""第八批临床输入准备与描述模块\n状态：已运行并验证；未拟合统计模型。\n脚本：{Path(__file__)}\n结果目录：{OUT}\n\n已核对：CSF200人/1389次，血浆191人/1248次，配对137人/923次。CSF227个临床合格原锚点中200人纳入，27人未纳入。\n27人的具体原因：{summary['CSF_exclusion_reasons']}。少于2次合格未来观察不直接等同于失访。\n两个分表连接：200人基线全部合格；未来1379次合格。排除3次缺分表、7次分表INFODT不一致。未挪动锚点或重标准化暴露；原1389次总分结果保持不变。\n\n供R使用：subpart_baseline.tsv和subpart_future.tsv，列为PATNO、EVENT_ID、NP1RTOT、NP1PTOT、updrs1_score。PATNO按字符串连接；两表PATNO+EVENT_ID唯一。\n全量资格台账：subpart_baseline_eligibility_ledger.tsv（200行）、subpart_future_eligibility_ledger.tsv（1389行）；排除记录另存subpart_excluded_future_records.tsv。\n\nTable1：table1_continuous.tsv、table1_categorical.tsv含机器可读数字；table1_display.tsv为中位数[IQR]或n(%)展示。年龄/病程以年，LEDD为原始每日总当量值而非log1p值。原字典未明确写物理剂量单位，不另擅标mg/day。BMI以原kg/m²尺度展示。\n纳排SMD只描述基线；连续变量用均值差除以两组样本方差均值的平方根，分类变量按每个类别的一对其余比例差计算。缺失列明，不插补、不增加P值。\n年度覆盖以(0,1]、(1,2]等真实观察人数及次数汇总，不称风险集；13条目按访视记录描述，0分比例分母为该条目合法观察数。\n输入SHA256、{len(CHECKS)}项通过的核对、软件版本和脚本哈希均在结果目录留存。\n\n运行方式：python batch8/scripts/prepare_clinical.py\n完成结果不覆盖；如需变更，先建立新的输出版本。\n"""
    (BATCH / "metadata/执行准备说明.txt").write_text(note, encoding="utf-8")
    print(json.dumps({"output": str(OUT), **summary}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
