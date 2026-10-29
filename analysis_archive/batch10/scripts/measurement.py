"""Fixed three-unit descriptive diagnostics and one outcome-blind background score."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
B=Path(__file__).resolve().parents[1];P=B.parent;O=B/'results/measurement_v1';O.mkdir(parents=True,exist_ok=False)
dtype={'PATNO':str,'EVENT_ID':str}
read=lambda p:pd.read_csv(p,sep='\t',dtype=dtype)
members=read(B/'protocol/feature_membership.tsv');feats=members.feature.tolist()
ismap={f:('LPC_d7' if f.startswith('LPC(') else 'LPE_d7') for f in feats}
summaries=[];pairs=[];unit_agreement=[];is_audits=[];is_links=[]
for assay in ['CSF','plasma']:
    t=read(B/f'input/target_units_{assay}.tsv')
    t.loc[~np.isfinite(t.TESTVALUE)|(t.TESTVALUE<=0),'TESTVALUE']=np.nan
    w=t.pivot(index=['PATNO','EVENT_ID'],columns=['UNITS','TESTNAME'],values='TESTVALUE')
    cols=pd.MultiIndex.from_product([['area','area_ratio','adjusted_area_ratio'],feats])
    w=w.reindex(columns=cols);n_before=len(w);w=w.dropna();assert len(w)>50
    standard=read(B/f'input/direct_standards_{assay}.tsv')
    inferred=np.log2(w['area']/w['area_ratio'])
    is_response={}
    for group in ['LPC_d7','LPE_d7']:
        fs=[f for f in feats if ismap[f]==group];z=inferred[fs]
        central=z.median(axis=1);is_response[group]=central
        spread=z.max(axis=1)-z.min(axis=1)
        iso=is_audits.append(dict(assay=assay,IS_group=group,n=len(w),max_within_sample_log2_spread=float(spread.max()),median_within_sample_log2_spread=float(spread.median()),direct_names=';'.join(standard.TESTNAME.unique())))
        direct=standard[standard.UNITS.eq('area')&standard.TESTNAME.str.startswith(group[:3]+'(')]
        if direct.TESTNAME.nunique()==1 and not direct.duplicated(['PATNO','EVENT_ID']).any():
            dv=direct.set_index(['PATNO','EVENT_ID']).TESTVALUE.reindex(w.index)
            valid=dv.gt(0)&dv.notna()
            is_audits[-1].update(direct_matched=int(valid.sum()),max_abs_direct_log2_difference=float((np.log2(dv[valid])-central[valid]).abs().max()))
    for unit in ['area','area_ratio','adjusted_area_ratio']:
        z=np.log2(w[unit][feats]);mat=z.corr(method='spearman')
        records=[]
        for i,f in enumerate(feats):
            for j,g in enumerate(feats[:i]):
                rec=dict(assay=assay,unit=unit,people=len(w),feature1=f,feature2=g,rho=mat.loc[f,g],same_internal_standard=ismap[f]==ismap[g]);pairs.append(rec);records.append(rec)
            is_links.append(dict(assay=assay,unit=unit,feature=f,n=len(w),rho_log_feature_vs_inferred_IS=z[f].corr(is_response[ismap[f]],method='spearman')))
        rr=pd.DataFrame(records)
        for scope,d in [('all_pairs',rr),('same_IS',rr[rr.same_internal_standard]),('different_IS',rr[~rr.same_internal_standard])]:
            summaries.append(dict(assay=assay,unit=unit,scope=scope,people=len(w),initial_paired_people=n_before,pairs=len(d),median_rho=d.rho.median(),min_rho=d.rho.min(),max_rho=d.rho.max()))
        # Descriptive linear partial correlation in log2 signals, controlling both inferred standards.
        X=np.column_stack([np.ones(len(w))]+[v.to_numpy() for v in is_response.values()])
        residual=z.to_numpy()-X@np.linalg.lstsq(X,z.to_numpy(),rcond=None)[0]
        cr=np.corrcoef(residual,rowvar=False)[np.tril_indices(len(feats),-1)]
        summaries.append(dict(assay=assay,unit=unit,scope='log2_Pearson_partial_both_inferred_IS',people=len(w),initial_paired_people=n_before,pairs=len(cr),median_rho=np.median(cr),min_rho=np.min(cr),max_rho=np.max(cr)))
    for f in feats:
        for u1,u2 in [('area','area_ratio'),('area_ratio','adjusted_area_ratio'),('area','adjusted_area_ratio')]:
            unit_agreement.append(dict(assay=assay,feature=f,unit1=u1,unit2=u2,people=len(w),spearman=w[u1][f].corr(w[u2][f],method='spearman')))
for name,d in [('correlation_summary',summaries),('all_pair_correlations',pairs),('within_feature_unit_concordance',unit_agreement),('inferred_IS_audit',is_audits),('inferred_IS_correlations',is_links)]:
    pd.DataFrame(d).to_csv(O/(name+'.tsv'),sep='\t',index=False)
# All 227 anchors define feature eligibility and scale before future-outcome filtering.
a=read(B/'input/clinical_anchors.tsv')[['PATNO','EVENT_ID']]
bg=read(B/'input/background_candidates.tsv')
bg.loc[~np.isfinite(bg.TESTVALUE)|(bg.TESTVALUE<=0),'TESTVALUE']=np.nan
v=bg.pivot(index=['PATNO','EVENT_ID'],columns='TESTNAME',values='TESTVALUE').reindex(pd.MultiIndex.from_frame(a))
z=np.log2(v);mean=z.mean();sd=z.std(ddof=1);coverage=z.notna().mean()
eligible=(coverage>=.95)&sd.gt(0)&np.isfinite(sd)
qual=pd.DataFrame(dict(feature=z.columns,coverage=coverage.values,mean_log2=mean.values,sd_log2=sd.values,eligible=eligible.values))
qual.to_csv(O/'background_feature_eligibility.tsv',sep='\t',index=False)
z=(z.loc[:,eligible]-mean[eligible])/sd[eligible]
assert z.shape[1]>=10 and not set(z.columns).intersection(feats)
observed=z.notna().sum(axis=1);valid=(observed/z.shape[1]>=.95)
bgmean=z.mean(axis=1).where(valid);mu=bgmean.mean();sigma=bgmean.std(ddof=1)
score=pd.DataFrame({'background':(bgmean-mu)/sigma,'background_members_observed':observed}).reset_index()
score.to_csv(B/'input/background_scores.tsv',sep='\t',index=False)
f=read(P/'batch9/results/classes_v1/model_data_CSF__LPI.tsv').drop_duplicates('PATNO')
link=f[['PATNO','x']].merge(score,on='PATNO',validate='one_to_one')
rho=link.x.corr(link.background);sp=link.x.corr(link.background,method='spearman')
summary=dict(anchor_n=len(a),candidate_features=len(qual),eligible_features=int(eligible.sum()),eligible_anchor_scores=int(valid.sum()),mean_before_final_standardization=float(mu),sd_before_final_standardization=float(sigma),original_model_people=len(link),model_people_with_background=int(link.background.notna().sum()),rho_LPI_background=float(rho),spearman_LPI_background=float(sp),two_predictor_VIF=float(1/(1-rho*rho)),interpretation='Common IS response can contribute to correlation; diagnostic cannot partition biological and technical sources of common variation. Background-adjusted association is conditional on a composite of other lipids, not proof of independence from all technical confounding.')
(O/'background_summary.json').write_text(json.dumps(summary,indent=2))
(O/'COMPLETE.txt').write_text('Descriptive measurement diagnostics and one fixed background score completed.\n')
print(json.dumps(summary,indent=2))
