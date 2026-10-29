"""Freeze an exploratory protocol, then build parent-lipid scores without clinical outcomes."""
from pathlib import Path
import datetime
import hashlib
import json
import re
import shutil
import numpy as np
import pandas as pd

B=Path(__file__).resolve().parents[1]; P=B.parent
for sub in ['protocol','input/parent_v2','results/preparation_v2','jobs','logs']:
    (B/sub).mkdir(parents=True,exist_ok=True)
O=B/'results/preparation_v2'; I=B/'input/parent_v2'
lock=B/'protocol/analysis_lock_v2.json'
assert not lock.exists(), 'Protocol and completed preparation must not be overwritten'
protocol={
 'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'stage':'post-selection exploratory extension; fixed before new parent-conditional models, not prospective preregistration',
 'question':'Does the frozen lysophospholipid association persist conditional on its corresponding parent phospholipid signal?',
 'pairs':{'LPI':'PI','LPC':'PC','LPE':'PE'},
 'supersedes':'v1 preparation stopped before any model: PI uses chain-level names, not sum-composition. v2 corrects annotation only; eligibility and statistical rules unchanged',
 'parent_annotation':'PI: unmodified two-chain names ^PI\\(\\d+:\\d+[_/]\\d+:\\d+\\)$; PC/PE: exact unmodified diacyl sum-composition ^(PC|PE)\\(\\d+:\\d+\\)$; exclude ether, oxidized and isotope standards; no sn-position interpretation',
 'feature_eligibility':'at least 95% finite positive adjusted_area_ratio on the frozen molecular anchors, positive sample SD of log2 signal',
 'parent_score':'mean member z(log2 adjusted_area_ratio), require all eligible members; standardize mean on complete frozen anchors, ddof=1',
 'lyso_score':'reuse batch9 frozen x and all clinical covariate transforms unchanged',
 'formula':'y ~ time*(x+parent+y0+up0+age0+sex+duration0+ledd0+genetic_group)+(1+time||PATNO)',
 'same_rows_comparator':'identical rows and formula without parent or time:parent',
 'inference':'REML, lmerTest Satterthwaite, bobyqa maxfun=100000; convergence/rank/singular/nonfinite fits not qualified; no automatic model substitution',
 'primary_family':'six time:x tests: 3 classes x 2 fluids; BH n=6, alongside nominal 95% CI',
 'paired_design':'reuse frozen 137 participants and 923 shared follow-up rows, restrict all 6 datasets together if any parent incomplete; standardize parents within that common paired cohort; frozen x unchanged',
 'paired_bootstrap':{'replicates':2000,'seed':20261030,'workers':8,'unit':'participant, same draw for all classes and both fluids, new cluster ID per resampled copy','pilot':20,'max_failure_fraction':0.05,'estimand':'conditional CSF time:x minus conditional plasma time:x','ci':[0.025,0.975,0.05/6,1-0.05/6],'ci_method':'percentile; fixed 3-contrast Bonferroni family interval 98.333333%; no test from separate fluid p-values'},
 'secondary':'one preselected CSF LPI + PI model for ESS using batch10 fixed ESS rows; descriptive exploratory p-value; not independent endpoint validation',
 'measurement_boundary':'scores index relative normalized signals; no cross-class molar ratios, enzyme activity, metabolic flux or causality claims',
 'stopping_rule':'complete fixed six models, three paired contrasts, one ESS model and bounded external qualification; no new endpoint or lipid search according to significance',
 'resources':{'combined_thread_ceiling':48,'server200_reserved_threads':9,'server230_compute_threads':0,'blas_threads':1,'planned_memory_GiB':16},
 'validation':'reproduce original six beta and SE to 1e-7; verify parent score arithmetic independently; exact participant/visit keys; pilot before full; all source hashes recorded'
}
lock.write_text(json.dumps(protocol,indent=2,ensure_ascii=False),encoding='utf-8')
(B/'protocol/LOCK_v2.sha256').write_text(hashlib.sha256(lock.read_bytes()).hexdigest()+'\n')
dtype={'PATNO':str,'EVENT_ID':str,'CLINICAL_EVENT':str}
read=lambda p:pd.read_csv(p,sep='\t',dtype=dtype)
manifest=[]; audit=[]; scales=[]; long=[]; scores=[]
def register(path):
    manifest.append(dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
for assay,fn in [('CSF','Metabolomic_Analysis_of_LRRK2_PD__CSF__27Oct2026.csv'),('plasma','Metabolomic_Analysis_of_LRRK2_PD_5_of_5_27Oct2026.csv')]:
    ap=P/f'batch9/results/classes_v1/anchor_scores_{assay}.tsv'
    anchors=read(ap)[['PATNO','EVENT_ID']];register(ap)
    assert not anchors.PATNO.duplicated().any()
    rp=P/'source_recheck_20261028_v2/inputs'/fn
    raw=pd.read_csv(rp,dtype=dtype,usecols=['PATNO','CLINICAL_EVENT','TYPE','TESTNAME','TESTVALUE','UNITS','RUNDATE']).rename(columns={'CLINICAL_EVENT':'EVENT_ID'})
    register(rp);assert raw.TYPE.nunique()==1
    selected=raw[raw.TESTNAME.str.fullmatch(r'PI\(\d+:\d+[_/]\d+:\d+\)|(PC|PE)\(\d+:\d+\)')].merge(anchors,on=['PATNO','EVENT_ID'],validate='many_to_one')
    assert not selected.duplicated(['PATNO','EVENT_ID','TESTNAME','UNITS']).any()
    selected.to_csv(I/f'parent_raw_units_{assay}.tsv',sep='\t',index=False)
    selected=selected[selected.UNITS.eq('adjusted_area_ratio')]
    s=anchors.copy();s['assay']=assay
    for parent in ['PI','PC','PE']:
        z={}
        for feature in sorted(f for f in selected.TESTNAME.unique() if f.startswith(parent+'(')):
            d=anchors.merge(selected[selected.TESTNAME.eq(feature)][['PATNO','EVENT_ID','TESTVALUE']],on=['PATNO','EVENT_ID'],how='left',validate='one_to_one')
            v=pd.to_numeric(d.TESTVALUE,errors='raise').to_numpy(float)
            ok=np.isfinite(v)&(v>0);lv=np.full(len(v),np.nan);lv[ok]=np.log2(v[ok]);mu=np.nanmean(lv);sd=np.nanstd(lv,ddof=1)
            eligible=ok.mean()>=.95 and np.isfinite(sd) and sd>0
            audit.append(dict(assay=assay,parent=parent,feature=feature,anchors=len(d),positive_n=int(ok.sum()),coverage=float(ok.mean()),eligible=bool(eligible),log2_mean=mu,log2_sd=sd))
            if eligible:
                z[feature]=(lv-mu)/sd
                q=d[['PATNO','EVENT_ID']].copy();q['assay']=assay;q['parent']=parent;q['feature']=feature;q['signal']=v;q['log2_signal']=lv;q['z_feature']=z[feature];long.append(q)
        assert z, f'No qualified features: {assay} {parent}'
        zz=np.column_stack(list(z.values()));complete=np.isfinite(zz).all(axis=1);mean=np.full(len(anchors),np.nan);mean[complete]=zz[complete].mean(axis=1)
        mu=np.nanmean(mean);sd=np.nanstd(mean,ddof=1);assert complete.sum()>=30 and sd>0
        s[parent]=(mean-mu)/sd
        scales.append(dict(assay=assay,parent=parent,members=';'.join(z),member_n=len(z),anchor_n=len(anchors),complete_anchor_n=int(complete.sum()),mean_member_z=mu,sd_member_z=sd))
    scores.append(s)
pd.DataFrame(audit).to_csv(O/'feature_eligibility.tsv',sep='\t',index=False)
pd.DataFrame([x for x in audit if x['eligible']]).to_csv(B/'protocol/parent_membership_v2.tsv',sep='\t',index=False)
pd.DataFrame(scales).to_csv(I/'parent_score_scaling.tsv',sep='\t',index=False)
pd.concat(long,ignore_index=True).to_csv(I/'parent_feature_values.tsv',sep='\t',index=False)
pd.concat(scores,ignore_index=True).to_csv(I/'parent_scores.tsv',sep='\t',index=False)
register(B/'protocol/parent_membership_v2.tsv')
for path in sorted(I.glob('*.tsv')): register(path)
(O/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
for name in ['job_runner.py','submit_job.py','launch.py']:
    shutil.copyfile(P/'batch10/scripts'/name,B/'scripts'/name)
(O/'COMPLETE.txt').write_text('Outcome-blind parent score preparation completed; model fitting not yet run.\n')
print(json.dumps(scales,indent=2))
