"""Correct isotope-name semantics before inspecting any new association result."""
from pathlib import Path
import datetime,json,hashlib
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parents[1];P=B.parent;O=B/'results/measurement_background_v2';O.mkdir(parents=True,exist_ok=False)
amend={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'supersedes':'Background feature eligibility in analysis_lock_v1 only; all endpoints, models and thresholds unchanged','reason':'Natural sphingoid-base d18:1 was erroneously excluded by the isotope pattern \\(d[0-9]. Correct isotope marking is a standalone (d7)/(d9)/etc, not the natural d18:1 chain prefix. Found by reviewing feature names. No new background-model coefficients or P values inspected before this correction. Prior background_v1 must not be used for final inference.','corrected_isotope_exclusion_regex':r'\(d[0-9]+\)|-d[0-9]+(?:_|$)','new_outputs':['measurement_background_v2','background_v2'],'source_lock_sha256':(B/'protocol/LOCK.sha256').read_text().strip()}
out=B/'protocol/background_annotation_amendment_v2.json';assert not out.exists();out.write_text(json.dumps(amend,indent=2))
raw=pd.read_csv(P/'source_recheck_20261028_v2/inputs/Metabolomic_Analysis_of_LRRK2_PD__CSF__27Oct2026.csv',dtype={'PATNO':str,'CLINICAL_EVENT':str}).rename(columns={'CLINICAL_EVENT':'EVENT_ID'})
a=pd.read_csv(B/'input/clinical_anchors.tsv',sep='\t',dtype={'PATNO':str})[['PATNO','EVENT_ID']]
feats=pd.read_csv(B/'protocol/feature_membership.tsv',sep='\t').feature.tolist()
classes='BMP|CE|Cer|DG|GalCer|GlcCer|GM1|GM2|GM3|GB3|HexCer|LacCer|MG|PA|PC|PE|PG|PI|PS|SM|TG|LPC|LPE|LPG|LPI'
bg=raw[raw.UNITS.eq('adjusted_area_ratio')].merge(a,on=['PATNO','EVENT_ID'],validate='many_to_one')
bg=bg[bg.TESTNAME.str.match(r'^('+classes+r')\(')&~bg.TESTNAME.isin(feats)&~bg.TESTNAME.str.contains(amend['corrected_isotope_exclusion_regex'],regex=True)]
assert not bg.duplicated(['PATNO','EVENT_ID','TESTNAME']).any()
bg=bg[['PATNO','EVENT_ID','TESTNAME','TESTVALUE']];bg.to_csv(B/'input/background_candidates_v2.tsv',sep='\t',index=False)
bg.loc[~np.isfinite(bg.TESTVALUE)|(bg.TESTVALUE<=0),'TESTVALUE']=np.nan
v=bg.pivot(index=['PATNO','EVENT_ID'],columns='TESTNAME',values='TESTVALUE').reindex(pd.MultiIndex.from_frame(a))
z=np.log2(v);mean=z.mean();sd=z.std(ddof=1);cov=z.notna().mean();ok=(cov>=.95)&sd.gt(0)&np.isfinite(sd)
qual=pd.DataFrame(dict(feature=z.columns,coverage=cov.values,mean_log2=mean.values,sd_log2=sd.values,eligible=ok.values))
qual.to_csv(O/'background_feature_eligibility.tsv',sep='\t',index=False)
z=(z.loc[:,ok]-mean[ok])/sd[ok];assert z.shape[1]>=10 and not set(feats).intersection(z.columns)
observed=z.notna().sum(axis=1);valid=(observed/z.shape[1]>=.95);bgmean=z.mean(axis=1).where(valid);mu=bgmean.mean();sigma=bgmean.std(ddof=1)
score=pd.DataFrame({'background':(bgmean-mu)/sigma,'background_members_observed':observed}).reset_index()
score.to_csv(B/'input/background_scores_v2.tsv',sep='\t',index=False)
f=pd.read_csv(P/'batch9/results/classes_v1/model_data_CSF__LPI.tsv',sep='\t',dtype={'PATNO':str}).drop_duplicates('PATNO')
link=f[['PATNO','x']].merge(score,on='PATNO',validate='one_to_one');rho=link.x.corr(link.background)
summary=dict(anchor_n=len(a),candidate_features=len(qual),eligible_features=int(ok.sum()),eligible_anchor_scores=int(valid.sum()),mean_before_final_standardization=float(mu),sd_before_final_standardization=float(sigma),original_model_people=len(link),model_people_with_background=int(link.background.notna().sum()),rho_LPI_background=float(rho),spearman_LPI_background=float(link.x.corr(link.background,method='spearman')),two_predictor_VIF=float(1/(1-rho*rho)),scope='Fixed structurally named non-target lipid background; mixture of biological and technical variation; not a pure technical covariate.')
(O/'background_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
