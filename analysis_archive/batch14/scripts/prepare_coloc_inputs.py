import os
"""Harmonize dense MBOAT7 regional statistics without significance/IV filtering."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,re
from scipy.stats import norm
B=Path(os.path.join(os.environ['PD_PROJECT_ROOT'], 'batch14'))
O=B/'results/harmonized_v1';O.mkdir(parents=True,exist_ok=False)
assert (B/'raw/remote_regions_v1/COMPLETE.json').exists()
assert (B/'raw/lpi_region_v1/COMPLETE.json').exists()
plan=json.loads((B/'protocol/analysis_plan_v1.json').read_text());lo,hi=plan['primary_window']
specs=[('LPI18','lipid',B/'raw/lpi_region_v1/GCST90318378.region.tsv.gz'),
       ('LPI16','lipid',B/'raw/lpi_region_v1/GCST90318519.region.tsv.gz'),
       ('cortex','expression',B/'raw/remote_regions_v1/QTD000171.MBOAT7.tsv.gz'),
       ('nigra','expression',B/'raw/remote_regions_v1/QTD000206.MBOAT7.tsv.gz'),
       ('PD','cc',B/'raw/remote_regions_v1/GCST009325.region.tsv.gz')]
summaries=[];manifests=[];data={}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name,kind,path in specs:
 d=pd.read_csv(path,sep='\t');nr=len(d);assert nr>0
 if kind=='expression':
  assert d.gene_id.str.split('.').str[0].eq('ENSG00000125505').all()
  x=pd.DataFrame({'chromosome':d.chromosome.astype(str),'position':d.position,'effect_allele_original':d.alt,'other_allele_original':d.ref,
    'beta_original':d.beta,'se':d.se,'eaf_original':d.ac/d.an,'p_original':d.pvalue,'n':d.an/2,'rsid':d.rsid,'source_variant_id':d.variant})
  x['source_ref']=d.ref;x['source_alt']=d.alt;x['cohorts']=np.nan
 else:
  x=pd.DataFrame({'chromosome':d.chromosome.astype(str),'position':d.base_pair_location,'effect_allele_original':d.effect_allele,'other_allele_original':d.other_allele,
    'beta_original':d.beta,'se':d.standard_error,'eaf_original':d.effect_allele_frequency,'p_original':d.p_value,
    'n':d.n if kind=='lipid' else d.N_cases+d.N_controls,'rsid':d.rs_id if kind=='lipid' else d.rsid,'source_variant_id':d.variant_id})
  x['cohorts']=d.Direction.str.count('[+-]') if kind=='lipid' else np.nan
  if kind=='lipid':x['HetISq']=d.HetISq;x['Direction']=d.Direction
  else:
   x['N_cases']=d.N_cases;x['N_controls']=d.N_controls;x['hm_code']=d.hm_code;x['hm_coordinate_conversion']=d.hm_coordinate_conversion
 flags={
  'invalid_coordinate':~(x.chromosome.eq('19')&x.position.between(lo,hi)&x.position.mod(1).eq(0)),
  'invalid_beta_se':~(np.isfinite(x.beta_original)&np.isfinite(x.se)&x.se.gt(0)),
  'invalid_alleles':~(x.effect_allele_original.str.fullmatch('[ACGT]+',na=False)&x.other_allele_original.str.fullmatch('[ACGT]+',na=False)&x.effect_allele_original.ne(x.other_allele_original)),
  'invalid_frequency_or_N': (~(np.isfinite(x.eaf_original)&x.eaf_original.gt(0)&x.eaf_original.lt(1)&np.isfinite(x.n)&x.n.gt(0))) if kind!='cc' else pd.Series(False,index=x.index)
 }
 reject=np.logical_or.reduce([v.to_numpy() for v in flags.values()]);bad=x.loc[reject].copy()
 if len(bad):
  bad['reason']=[';'.join(k for k,v in flags.items() if v.loc[i]) for i in bad.index]
  bad.to_csv(O/f'{name}.excluded.tsv',sep='\t',index=False)
 x=x.loc[~reject].copy()
 assert len(x)>0
 # Lexicographically larger explicit allele is the harmonized effect allele.
 order=x.effect_allele_original.ge(x.other_allele_original)
 x['effect_allele']=np.where(order,x.effect_allele_original,x.other_allele_original)
 x['other_allele']=np.where(order,x.other_allele_original,x.effect_allele_original)
 flipped=x.effect_allele.ne(x.effect_allele_original)
 x['beta']=np.where(flipped,-x.beta_original,x.beta_original)
 x['eaf']=np.where(flipped,1-x.eaf_original,x.eaf_original)
 x['maf']=np.minimum(x.eaf,1-x.eaf)
 x['snp']='19:'+x.position.astype(int).astype(str)+':'+x.other_allele+':'+x.effect_allele
 x['palindromic']=(x.other_allele+x.effect_allele).isin(['AT','CG'])
 x['is_snv']=x.effect_allele.str.len().eq(1)&x.other_allele.str.len().eq(1)
 x['p_from_beta_se']=2*norm.sf(np.abs(x.beta/x.se))
 x['allele_swap']=flipped
 dup=x[x.duplicated('snp',keep=False)]
 if len(dup):
  conflict=dup.groupby('snp')[['beta','se','eaf','n','p_original']].nunique(dropna=False).max(axis=1)>1
  if conflict.any():
   dup[dup.snp.isin(conflict[conflict].index)].to_csv(O/f'{name}.conflicting_duplicates.tsv',sep='\t',index=False)
   raise ValueError(name+' has conflicting duplicate variants; inspect before proceeding')
 ndups=int(x.duplicated('snp').sum());x=x.drop_duplicates('snp').sort_values(['position','snp']).reset_index(drop=True)
 x.to_csv(O/f'{name}.tsv.gz',sep='\t',index=False)
 data[name]=x
 inv=1/x.se.to_numpy()**2;ng=2*x.n.to_numpy()*x.maf.to_numpy()*(1-x.maf.to_numpy())
 est_sdy=float(np.sqrt(np.sum(inv*ng)/np.sum(inv**2))) if kind!='cc' else None
 scale_quant=np.quantile(np.sqrt(ng*x.se.to_numpy()**2),[.05,.25,.5,.75,.95]).tolist() if kind!='cc' else None
 best=x.sort_values('p_from_beta_se').iloc[0]
 summaries.append({'trait':name,'kind':kind,'raw_rows':nr,'valid_unique_rows':len(x),'excluded_rows':int(reject.sum()),'duplicate_rows_collapsed':ndups,
   'positions':[int(x.position.min()),int(x.position.max())],'n_min':float(x.n.min()),'n_median':float(x.n.median()),'n_max':float(x.n.max()),
   'n_unique':int(x.n.nunique()),'snv_count':int(x.is_snv.sum()),'palindromic_count':int(x.palindromic.sum()),
   'min_p_original':float(x.p_original.min()),'min_p_beta_se':float(best.p_from_beta_se),'lead':best.snp,
   'sdY_estimated':est_sdy,'sdY_main':1 if kind=='expression' else est_sdy,'scale_quantiles':scale_quant,
   'exclusion_reason_counts':{k:int(v.sum()) for k,v in flags.items()},'significance_filter_used':False,'MR_N_cohort_gate_used':False})
 manifests.append({'trait':name,'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size})
pairs=[('cortex','LPI18'),('cortex','LPI16'),('nigra','LPI18'),('nigra','LPI16'),('cortex','PD'),('nigra','PD'),('LPI18','PD'),('LPI16','PD')]
coverage=[]
for a,b in pairs:
 m=data[a].merge(data[b],on='snp',suffixes=('_a','_b'));assert len(m)>0
 m['eaf_absolute_difference']=abs(m.eaf_a-m.eaf_b)
 m.to_csv(O/f'{a}__{b}.joined.tsv.gz',sep='\t',index=False)
 coverage.append({'trait1':a,'trait2':b,'n_overlap':len(m),'fraction_trait1':len(m)/len(data[a]),'fraction_trait2':len(m)/len(data[b]),
   'freq_diff_gt_0_2':int(m.eaf_absolute_difference.gt(.2).sum()),'position_min':int(m.position_a.min()),'position_max':int(m.position_a.max())})
pd.DataFrame(coverage).to_csv(O/'pair_coverage.tsv',sep='\t',index=False)
(O/'trait_qc.json').write_text(json.dumps(summaries,indent=2,allow_nan=False))
(O/'source_manifest.json').write_text(json.dumps(manifests,indent=2))
(O/'COMPLETE.json').write_text(json.dumps({'traits':5,'pairs':8,'complete':True}))
print(json.dumps({'traits':summaries,'coverage':coverage},indent=2,allow_nan=False))
