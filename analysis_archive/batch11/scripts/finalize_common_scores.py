"""Require identical qualified parent species in both fluids before any new fit."""
from pathlib import Path
import datetime,hashlib,json
import pandas as pd
import numpy as np
B=Path(__file__).resolve().parents[1];I=B/'input/parent_v3';I.mkdir(exist_ok=False)
protocol=json.loads((B/'protocol/analysis_lock_v2.json').read_text(encoding='utf-8'))
protocol['created_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
protocol['supersedes']='v2 preparation only, no outcomes modeled; one PI species fails CSF coverage. Require intersection of eligible members across fluids to make parent definitions identical.'
protocol['feature_eligibility']+='; use intersection of qualifying species across CSF and plasma for every full and paired model'
lock=B/'protocol/analysis_lock_v3.json';assert not lock.exists()
lock.write_text(json.dumps(protocol,indent=2,ensure_ascii=False),encoding='utf-8')
(B/'protocol/LOCK_v3.sha256').write_text(hashlib.sha256(lock.read_bytes()).hexdigest()+'\n')
prev=B/'input/parent_v2';dtype={'PATNO':str,'EVENT_ID':str}
read=lambda p:pd.read_csv(p,sep='\t',dtype=dtype)
members=read(B/'protocol/parent_membership_v2.tsv')
common=members.groupby(['parent','feature']).assay.nunique();common=common[common.eq(2)].index
members=members.set_index(['parent','feature']).loc[common].reset_index()
fv=read(prev/'parent_feature_values.tsv').merge(members[['assay','parent','feature']],on=['assay','parent','feature'],validate='many_to_one')
all_scores=[];scales=[]
for assay in ['CSF','plasma']:
    a=read(prev/'parent_scores.tsv').query('assay == @assay')[['PATNO','EVENT_ID','assay']]
    for parent in ['PI','PC','PE']:
        fs=sorted(members.query('assay == @assay and parent == @parent').feature)
        wide=fv.query('assay == @assay and parent == @parent').pivot(index=['PATNO','EVENT_ID'],columns='feature',values='z_feature').reindex(columns=fs)
        wide=wide.reindex(pd.MultiIndex.from_frame(a[['PATNO','EVENT_ID']]))
        val=wide.mean(axis=1,skipna=False);mu=val.mean();sd=val.std(ddof=1)
        a[parent]=((val-mu)/sd).to_numpy()
        scales.append(dict(assay=assay,parent=parent,members=';'.join(fs),member_n=len(fs),anchor_n=len(a),complete_anchor_n=int(val.notna().sum()),mean_member_z=mu,sd_member_z=sd))
    all_scores.append(a)
members.to_csv(B/'protocol/parent_membership_v3.tsv',sep='\t',index=False)
fv.to_csv(I/'parent_feature_values.tsv',sep='\t',index=False)
pd.concat(all_scores).to_csv(I/'parent_scores.tsv',sep='\t',index=False)
pd.DataFrame(scales).to_csv(I/'parent_score_scaling.tsv',sep='\t',index=False)
paths=[lock,B/'protocol/parent_membership_v3.tsv',*sorted(I.glob('*.tsv'))]
manifest=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in paths]
(I/'manifest.json').write_text(json.dumps(manifest,indent=2))
print(pd.DataFrame(scales).to_string(index=False))
