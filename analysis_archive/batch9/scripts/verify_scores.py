"""Independently reconstruct frozen score scaling and analysis membership in Python."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

B=Path(__file__).resolve().parents[1]
O=B/'results/classes_v1'
checks=[]
def check(name,ok):
    assert ok,name
    checks.append(name)
def read(n):
    return pd.read_csv(O/(n+'.tsv'),sep='\t',dtype={'PATNO':str,'EVENT_ID':str})
def bh(p):
    p=np.array(p,float);valid=np.isfinite(p);out=np.full(len(p),np.nan)
    ix=np.flatnonzero(valid);order=ix[np.argsort(p[ix])]
    out[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(order)+1))[::-1])[::-1])
    return out
membership=pd.read_csv(B/'protocol/feature_membership.tsv',sep='\t')
features=list(membership.feature)
fv=read('anchor_feature_values');sc=read('score_scaling')
check('ten fixed members',len(features)==10 and membership['class'].value_counts().to_dict()=={'LPC':6,'LPE':2,'LPI':2})
for a in ['CSF','plasma']:
    v=fv[fv.assay.eq(a)]
    scores=read('anchor_scores_'+a).set_index('PATNO')
    mat=v.pivot(index='PATNO',columns='feature',values='z_feature').reindex(scores.index)[features]
    check(a+' unique anchor feature keys',not v.duplicated(['PATNO','EVENT_ID','feature']).any())
    for f in features:
        z=mat[f].dropna()
        check(a+f+' feature standardization',abs(z.mean())<1e-10 and abs(z.std(ddof=1)-1)<1e-10)
    for row in sc[sc.assay.eq(a)].itertuples():
        fs=row.members.split(';')
        complete=mat[features if row.score=='LYSO_without_LPI16' else fs].notna().all(axis=1)
        mean=mat[fs].mean(axis=1).where(complete)
        target=(mean-mean.mean())/mean.std(ddof=1)
        check(a+row.score+' frozen complete anchor count',int(complete.sum())==row.complete_anchor_people)
        check(a+row.score+' score reconstruction',np.allclose(target,scores[row.score],atol=1e-10,rtol=1e-10,equal_nan=True))
        check(a+row.score+' scale reconstruction',abs(mean.mean()-row.mean_of_member_z)<1e-10 and abs(mean.std(ddof=1)-row.sd_of_member_z)<1e-10)
        if a=='plasma' and row.score=='LYSO_without_LPI16':continue
        d=read('model_data_'+a+'__'+row.score)
        check(a+row.score+' model exposure equals anchor score',np.allclose(d.x,scores[row.score].reindex(d.PATNO),atol=1e-10))
        check(a+row.score+' future independent subject visits',not d.duplicated(['PATNO','time']).any() and d.time.gt(0).all() and d.groupby('PATNO').time.nunique().ge(2).all())
    u=read('model_data_'+a+'__LYSO_ALL')
    if a=='CSF':
        minus=read('model_data_CSF__LYSO_without_LPI16')
        pd.testing.assert_frame_equal(u.drop(columns='x'),minus.drop(columns='x'))
        checks.append('leave-LPI16-out uses identical analysis rows')
pair=read('paired_anchor_scores').set_index('PATNO')
check('paired people unique',pair.index.is_unique and len(pair)>=50)
for cl in ['LPI','LPC','LPE','LYSO_ALL']:
    ds=[]
    for a in ['CSF','plasma']:
        d=read('paired_model_data_'+a+'_'+cl)
        s=pair[cl+'_'+a];z=(s-s.mean())/s.std(ddof=1)
        check(a+cl+' paired frozen exposure reconstruction',np.allclose(d.x,z.reindex(d.PATNO),atol=1e-10))
        check(a+cl+' paired subject and rows',d.PATNO.nunique()==len(pair) and len(d)==923)
        ds.append(d)
    pd.testing.assert_frame_equal(ds[0].drop(columns='x'),ds[1].drop(columns='x'))
    checks.append(cl+' exact shared clinical rows')
r=read('primary_score_results')
check('nine planned score models',len(r)==9 and r.status.eq('tested').all())
check('BH family of nine independently reproduced',np.allclose(r.q_family9,bh(r.p),atol=1e-12))
pilot=read('paired_pilot')
check('pilot has 20 draws per class',len(pilot)==80 and not pilot.duplicated(['iter','score']).any())
check('pilot <=5 percent failures for every class',pilot.groupby('score').qualified.mean().ge(.95).all())
out={'state':'passed','checks_passed':len(checks),'checks':checks}
(O/'python_score_verification.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps({'checks_passed':len(checks)}))
