import os
"""Read-only raw inputs; one dataset per expression phenotype; no significance filter."""
from pathlib import Path
import json,hashlib,gc
import numpy as np
import pandas as pd
from scipy.stats import norm
B=Path(os.path.join(os.environ['PD_PROJECT_ROOT'], 'batch15'))
O=B/'results/harmonized_v1';O.mkdir(exist_ok=False)
assert (B/'raw/lpi_regions_v1/COMPLETE.json').exists() and (B/'raw/remote_regions_v3/COMPLETE.json').exists()
plan=json.loads((B/'protocol/analysis_plan_v1.json').read_text()); manifest=[];qc=[];sources=[];exclusions=[]

def harmonize(d,kind,reg,identifier):
 if kind=='expression':
  x=pd.DataFrame({'chromosome':d.chromosome.astype(str),'position':d.position,'effect_allele_original':d.alt,'other_allele_original':d.ref,'beta_original':d.beta,'se':d.se,'eaf_original':d.ac/d.an,'n':d.an/2,'p_original':d.pvalue,'rsid':d.rsid})
  x['median_tpm']=d.median_tpm;x['minor_allele_samples']=d.ma_samples;x['imputation_r2']=d.r2;x['allele_number']=d.an;x['alt_count']=d.ac
 else:
  x=pd.DataFrame({'chromosome':d.chromosome.astype(str),'position':d.base_pair_location,'effect_allele_original':d.effect_allele,'other_allele_original':d.other_allele,'beta_original':d.beta,'se':d.standard_error,'eaf_original':d.effect_allele_frequency,'n':d.n if kind=='lipid' else d.N_cases+d.N_controls,'p_original':d.p_value,'rsid':d.rs_id if kind=='lipid' else d.rsid})
  if kind=='lipid':
   x['Direction']=d.Direction;x['HetISq']=d.HetISq
  else:
   x['N_cases']=d.N_cases;x['N_controls']=d.N_controls;x['hm_code']=d.hm_code
 flags={
  'coordinate':~(x.chromosome.eq(reg['chromosome'])&x.position.between(reg['start'],reg['end'])&x.position.mod(1).eq(0)),
  'beta_se':~(np.isfinite(x.beta_original)&np.isfinite(x.se)&x.se.gt(0)),
  'p_value':~(np.isfinite(x.p_original)&x.p_original.gt(0)&x.p_original.le(1)),
  'alleles':~(x.effect_allele_original.str.fullmatch('[ACGT]+',na=False)&x.other_allele_original.str.fullmatch('[ACGT]+',na=False)&x.effect_allele_original.ne(x.other_allele_original)),
  'frequency_n':~(np.isfinite(x.eaf_original)&x.eaf_original.gt(0)&x.eaf_original.lt(1)&np.isfinite(x.n)&x.n.gt(0)) if kind!='cc' else pd.Series(False,index=x.index)
 }
 reject=np.logical_or.reduce([v.to_numpy() for v in flags.values()]);bad=x.loc[reject].copy()
 if len(bad):
  bad['reason']=[';'.join(k for k,v in flags.items() if v.loc[i]) for i in bad.index];bad['dataset']=identifier;exclusions.append(bad)
 x=x.loc[~reject].copy()
 if x.empty:return x,{'raw_rows':len(d),'excluded_rows':len(d),'duplicate_rows':0,'valid_unique_rows':0}
 order=x.effect_allele_original.ge(x.other_allele_original)
 x['effect_allele']=np.where(order,x.effect_allele_original,x.other_allele_original);x['other_allele']=np.where(order,x.other_allele_original,x.effect_allele_original)
 flip=x.effect_allele.ne(x.effect_allele_original);x['beta']=np.where(flip,-x.beta_original,x.beta_original);x['eaf']=np.where(flip,1-x.eaf_original,x.eaf_original)
 x['maf']=np.minimum(x.eaf,1-x.eaf);x['snp']=x.chromosome+':'+x.position.astype(int).astype(str)+':'+x.other_allele+':'+x.effect_allele
 x['p_from_beta_se']=2*norm.sf(np.abs(x.beta/x.se));x['allele_swap']=flip;x['is_snv']=x.effect_allele.str.len().eq(1)&x.other_allele.str.len().eq(1)
 dup=x[x.duplicated('snp',keep=False)]
 if len(dup):
  conflict=dup.groupby('snp')[['beta','se','eaf','n','p_original']].nunique(dropna=False).max(axis=1)>1
  if conflict.any():
   dup[dup.snp.isin(conflict[conflict].index)].to_csv(O/f'{identifier}.conflicting_duplicates.tsv',sep='\t',index=False)
   raise ValueError(identifier+' conflicting duplicate statistics')
 nd=int(x.duplicated('snp').sum());x=x.drop_duplicates('snp').sort_values(['position','snp']).reset_index(drop=True)
 sd=None
 if kind!='cc':
  inv=1/x.se.to_numpy()**2;ng=2*x.n.to_numpy()*x.maf.to_numpy()*(1-x.maf.to_numpy());sd=float(np.sqrt(np.dot(inv,ng)/np.dot(inv,inv)))
 result={'raw_rows':len(d),'excluded_rows':int(reject.sum()),'duplicate_rows':nd,'valid_unique_rows':len(x),'sdY_estimated':sd,
  'min_p_original':float(x.p_original.min()),'min_p_beta_se':float(x.p_from_beta_se.min()),'n_min':float(x.n.min()),'n_max':float(x.n.max()),
  'position_min':int(x.position.min()),'position_max':int(x.position.max()),'fixed_lipid_lead_present':bool(x.snp.eq(reg['lead_key']).any())}
 return x,result

def save_dataset(d,kind,reg,identifier,source_id,gene='',phenotype=''):
 x,summary=harmonize(d,kind,reg,identifier)
 path=O/f'{identifier}.tsv.gz';x.to_csv(path,sep='\t',index=False)
 rec={'dataset':identifier,'region_id':reg['id'],'source_id':source_id,'kind':kind,'gene_id':gene,'phenotype':phenotype,'path':str(path),
      'sdY_main':reg['sdY_main'] if kind=='lipid' else 1 if kind=='expression' else None,**summary}
 qc.append(rec);manifest.append(rec);return x

def read_raw(path):
 sources.append({'path':str(path),'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 return pd.read_csv(path,sep='\t',dtype={'chromosome':str})

for reg in plan['regions']:
 d=read_raw(B/'raw/lpi_regions_v1'/f"{reg['id']}.tsv.gz")
 save_dataset(d,'lipid',reg,reg['id']+'__lipid',reg['study']);del d
 d=read_raw(B/'raw/remote_regions_v3'/f"{reg['id']}__GCST009325.tsv.gz")
 save_dataset(d,'cc',reg,reg['id']+'__PD','GCST009325');del d
 for src in plan['expression']:
  d=read_raw(B/'raw/remote_regions_v3'/f"{reg['id']}__{src}.tsv.gz")
  assert not d[['gene_id','molecular_trait_id']].isna().any().any()
  for (gene,phen),part in d.groupby(['gene_id','molecular_trait_id'],sort=True):
   assert part.gene_id.nunique()==1 and part.molecular_trait_id.nunique()==1
   tag=hashlib.sha256(str(phen).encode()).hexdigest()[:10]
   identifier=reg['id']+'__'+src+'__'+str(gene)+'__'+tag
   save_dataset(part,'expression',reg,identifier,src,str(gene),str(phen))
  del d;gc.collect()
pd.DataFrame(manifest).to_csv(O/'dataset_manifest.tsv',sep='\t',index=False)
(O/'dataset_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False))
(O/'source_manifest.json').write_text(json.dumps(sources,indent=2))
if exclusions:pd.concat(exclusions).to_csv(O/'excluded_rows.tsv.gz',sep='\t',index=False)
result={'complete':True,'datasets':len(manifest),'expression_phenotypes':sum(x['kind']=='expression' for x in manifest),'raw_sources':len(sources),'excluded_rows':sum(x['excluded_rows'] for x in manifest),'exact_duplicates_collapsed':sum(x['duplicate_rows'] for x in manifest)}
(O/'COMPLETE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
print(pd.DataFrame(manifest).groupby(['region_id','source_id']).agg(phenotypes=('dataset','size'),variants=('valid_unique_rows','sum')).to_string())
