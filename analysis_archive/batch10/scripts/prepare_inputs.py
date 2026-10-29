import os
"""Outcome-blind measurement extraction plus exact visit-key clinical qualification."""
from pathlib import Path
import hashlib,json,re
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parents[1];P=B.parent;O=B/'results/preparation_v1';O.mkdir(parents=True,exist_ok=True)
assert not (O/'source_manifest.json').exists(), 'Completed preparation is immutable'
lock=json.loads((B/'protocol/analysis_lock_v1.json').read_text(encoding='utf-8'))
assert hashlib.sha256((B/'protocol/analysis_lock_v1.json').read_bytes()).hexdigest()==(B/'protocol/LOCK.sha256').read_text().strip()
dtype={'PATNO':str,'EVENT_ID':str,'CLINICAL_EVENT':str}
read=lambda p:pd.read_csv(p,sep='\t',dtype=dtype)
members=read(B/'protocol/feature_membership.tsv');feats=members.feature.tolist()
paired=read(P/'batch9/results/classes_v1/paired_anchor_scores.tsv')[['PATNO','EVENT_ID']]
clinical=Path(os.path.join(os.environ['PD_EVIDENCE_ROOT'], 'ppmi_clinical_analysis_20261023'))
c=read(clinical/'cleaned/curated_original_snapshot.tsv.gz')
m=read(clinical/'cleaned/clinical_model_ready.tsv.gz')
f=read(P/'batch9/results/classes_v1/model_data_CSF__LPI.tsv')
a=read(P/'batch9/results/classes_v1/anchor_scores_CSF.tsv')
assert len(f)==1389 and f.PATNO.nunique()==200 and len(a)==227
assert not c.duplicated(['PATNO','EVENT_ID']).any() and not m.duplicated(['PATNO','EVENT_ID']).any()
months=pd.to_datetime(c.visit_date,format='%m/%Y',errors='coerce')
c['curated_month']=months.dt.year*12+months.dt.month-1
ff=f.merge(c[['PATNO','EVENT_ID','curated_month','moca','gds','scopa','ess','LEDD']],on=['PATNO','EVENT_ID'],how='left',validate='one_to_one')
aa=a[['PATNO','EVENT_ID','visit_month']].merge(c[['PATNO','EVENT_ID','curated_month','moca','gds','scopa','ess','LEDD']],on=['PATNO','EVENT_ID'],how='left',validate='one_to_one')
assert np.all(ff.visit_month==ff.curated_month) and np.all(aa.visit_month==aa.curated_month)
fm=f[['PATNO','EVENT_ID']].merge(m[['PATNO','EVENT_ID','ledd']],on=['PATNO','EVENT_ID'],validate='one_to_one')
assert np.allclose(ff.LEDD, fm.ledd,equal_nan=True)
baseline=aa[aa.PATNO.isin(f.PATNO)].set_index('PATNO')
assert np.allclose(f.ledd0,np.log1p(f.PATNO.map(baseline.LEDD)))
ff.rename(columns={'LEDD':'ledd_visit'},inplace=True)
ff.to_csv(B/'input/clinical_future.tsv',sep='\t',index=False)
aa.rename(columns={k:k+'0' for k in ['moca','gds','scopa','ess']}).to_csv(B/'input/clinical_anchors.tsv',sep='\t',index=False)
qa=[]
for k,s in lock['clinical_endpoints'].items():
    lo,hi=s['range']
    for label,d in [('baseline_200',aa[aa.PATNO.isin(f.PATNO)]),('future_1389',ff)]:
        x=d[k];bad=x.notna()&((x<lo)|(x>hi)|(x%1!=0))
        qa.append(dict(endpoint=k,scope=label,rows=len(d),nonmissing=int(x.notna().sum()),invalid=int(bad.sum()),minimum=x.min(),maximum=x.max()))
pd.DataFrame(qa).to_csv(O/'clinical_range_audit.tsv',sep='\t',index=False)
(O/'clinical_join_audit.json').write_text(json.dumps({'exact_event_joins':True,'baseline_months_equal':True,'future_months_equal':True,'visit_LEDD_matches_cleaned':True,'baseline_LEDD_transform_reproduced':True,'LEDD_min':float(ff.ledd_visit.min()),'LEDD_max':float(ff.ledd_visit.max()),'date_precision':'month; outcome and molecular anchor paired by original visit; no day-level timing claim'},indent=2))
sources=[]
for assay,name in [('CSF','Metabolomic_Analysis_of_LRRK2_PD__CSF__27Oct2026.csv'),('plasma','Metabolomic_Analysis_of_LRRK2_PD_5_of_5_27Oct2026.csv')]:
    path=P/'source_recheck_20261028_v2/inputs'/name
    raw=pd.read_csv(path,dtype=dtype).rename(columns={'CLINICAL_EVENT':'EVENT_ID'})
    assert raw.TYPE.nunique()==1
    targets=raw[raw.TESTNAME.isin(feats)].merge(paired,on=['PATNO','EVENT_ID'],validate='many_to_one')
    assert not targets.duplicated(['PATNO','EVENT_ID','TESTNAME','UNITS']).any()
    targets[['PATNO','EVENT_ID','TESTNAME','TESTVALUE','UNITS','RUNDATE']].to_csv(B/f'input/target_units_{assay}.tsv',sep='\t',index=False)
    # Direct standard rows are inspected as a measurement cross-check, never outcome variables.
    standard=raw[raw.TESTNAME.str.contains(r'LPE\(18:1\(d7\)\)|LPC\(18:1\(d7\)\)',regex=True)].merge(paired,on=['PATNO','EVENT_ID'],validate='many_to_one')
    standard[['PATNO','EVENT_ID','TESTNAME','TESTVALUE','UNITS','RUNDATE']].to_csv(B/f'input/direct_standards_{assay}.tsv',sep='\t',index=False)
    if assay=='CSF':
        bg=raw[raw.UNITS.eq('adjusted_area_ratio')].merge(a[['PATNO','EVENT_ID']],on=['PATNO','EVENT_ID'],validate='many_to_one')
        classes='BMP|CE|Cer|DG|GalCer|GlcCer|GM1|GM2|GM3|GB3|HexCer|LacCer|MG|PA|PC|PE|PG|PI|PS|SM|TG|LPC|LPE|LPG|LPI'
        bg=bg[bg.TESTNAME.str.match(r'^('+classes+r')\(')&~bg.TESTNAME.isin(feats)&~bg.TESTNAME.str.contains(r'\(d[0-9]|d[0-9]_',regex=True)]
        assert not bg.duplicated(['PATNO','EVENT_ID','TESTNAME']).any()
        bg[['PATNO','EVENT_ID','TESTNAME','TESTVALUE']].to_csv(B/'input/background_candidates.tsv',sep='\t',index=False)
    sources.append({'path':str(path),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':len(raw)})
for path in [clinical/'cleaned/curated_original_snapshot.tsv.gz',clinical/'cleaned/clinical_model_ready.tsv.gz',P/'batch9/results/classes_v1/model_data_CSF__LPI.tsv',P/'batch9/results/classes_v1/anchor_scores_CSF.tsv']:
    sources.append({'path':str(path),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
(O/'source_manifest.json').write_text(json.dumps(sources,indent=2))
print('Prepared exact joins, raw unit targets, and outcome-blind background candidates')
