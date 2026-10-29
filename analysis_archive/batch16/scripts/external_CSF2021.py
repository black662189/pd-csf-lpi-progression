import os
from pathlib import Path
import json, hashlib, sys
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(os.environ['PD_PROJECT_ROOT']);B=P/'batch16';O=B/'external_analysis_v1';O.mkdir(exist_ok=False)
plan=B/'protocol/external_CSF2021_plan_v1.json';src=B/'external_search_v1/CSF2021_targets_metadata.tsv'
(O/'LOCK.sha256').write_text(hashlib.sha256(plan.read_bytes()).hexdigest()+'\n')
manifest={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [plan,src,B/'external_search_v1/CSF2021_qualification.json']}
(O/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
d=pd.read_csv(src,sep='\t');assert len(d)==30 and d.sample_id.nunique()==30 and d.group.value_counts().to_dict()=={'PD':20,'Control':10}
assert d[['age','sex','LPI16_0','LPI18_0']].notna().all().all() and d[['LPI16_0','LPI18_0']].gt(0).all().all()
d['PD']=d.group.eq('PD').astype(int);d['sex_M']=d.sex.map({'M':1,'F':0});d['age10']=(d.age-d.age.mean())/10
scaling=[]
for f in ['LPI16_0','LPI18_0']:
 v=np.log2(d[f]);d[f+'_z']=(v-v.mean())/v.std(ddof=1);scaling.append({'variable':f,'log2_mean':v.mean(),'log2_SD':v.std(ddof=1),'n':30})
v=d[['LPI16_0_z','LPI18_0_z']].mean(axis=1);d['LPI_score']=(v-v.mean())/v.std(ddof=1);scaling.append({'variable':'mean_species_z','mean':v.mean(),'SD':v.std(ddof=1),'n':30})
d['stage']=d.neuropathology.map({'LBDL':1,'LBDE':2,'LBDN':3})
d.to_csv(O/'adapted_scores.csv',index=False);pd.DataFrame(scaling).to_csv(O/'scaling.csv',index=False)
rows=[];verification=[]
for comparison,dd,predictor in [('PD_vs_control',d,'PD'),('PD_pathology_ordinal',d[d.PD.eq(1)],'stage')]:
 for target in ['LPI_score','LPI16_0_z','LPI18_0_z']:
  X=sm.add_constant(dd[[predictor,'age10','sex_M']].astype(float));y=dd[target].to_numpy();assert np.linalg.matrix_rank(X)==4
  fit=sm.OLS(y,X).fit(cov_type='HC3',use_t=True);a=fit.conf_int().loc[predictor]
  fam='primary_auxiliary_1' if comparison=='PD_vs_control' and target=='LPI_score' else 'species_PD_family2' if comparison=='PD_vs_control' else 'pathology_family3'
  rows.append({'comparison':comparison,'outcome':target,'predictor':predictor,'people':len(dd),'beta':fit.params[predictor],'se_HC3':fit.bse[predictor],'df':fit.df_resid,'lo':a.iloc[0],'hi':a.iloc[1],'p':fit.pvalues[predictor],'family':fam,'q':np.nan,'adjustment':'age10 + sex_M','scale_reference':'all30 postmortem donors; external-adapted score'})
  xx=X.to_numpy();bread=np.linalg.inv(xx.T@xx);coef=bread@xx.T@y;err=y-xx@coef;h=np.einsum('ij,jk,ik->i',xx,bread,xx);meat=xx.T@np.diag((err/(1-h))**2)@xx;cov=bread@meat@bread
  ce=float(np.max(abs(coef-fit.params.to_numpy())));ve=float(np.max(abs(cov-fit.cov_params().to_numpy())));assert ce<1e-10 and ve<1e-10
  verification.append({'comparison':comparison,'outcome':target,'coefficient_max_error':ce,'HC3_covariance_max_error':ve,'pass':True})
r=pd.DataFrame(rows)
for fam,n in [('species_PD_family2',2),('pathology_family3',3)]:
 mask=r.family.eq(fam);assert mask.sum()==n;r.loc[mask,'q']=multipletests(r.loc[mask,'p'],method='fdr_bh')[1]
assert len(r)==6 and np.isfinite(r[['beta','se_HC3','df','lo','hi','p']]).all().all()
r.to_csv(O/'external_results.csv',index=False);r.to_csv(B/'delivery_v1/tables/Table_S34_external_postmortem_CSF.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(1,2,figsize=(10.5,4.6));fig.subplots_adjust(wspace=.38,bottom=.28,top=.88)
rng=np.random.default_rng(20261029);colors=['#67727D','#138D90'];labels=['Control','PD']
for i,name in enumerate(labels):
 vals=d.loc[d.group.eq(name),'LPI_score'];ax[0].scatter(i+rng.uniform(-.12,.12,len(vals)),vals,s=34,color=colors[i],alpha=.8,edgecolor='white',linewidth=.5);ax[0].plot([i-.15,i+.15],[vals.median()]*2,color='black',lw=2)
ax[0].set(xticks=[0,1],xticklabels=['Control\nn=10','PD\nn=20'],ylabel='External-adapted LPI score (SD)',title='A  Independent postmortem CSF cohort')
cols=['#D09253','#138D90','#244C66']
for i,stage in enumerate(['LBDL','LBDE','LBDN']):
 vals=d.loc[d.neuropathology.eq(stage),'LPI_score'];ax[1].scatter(i+rng.uniform(-.12,.12,len(vals)),vals,s=34,color=cols[i],alpha=.8,edgecolor='white',linewidth=.5);ax[1].plot([i-.15,i+.15],[vals.median()]*2,color='black',lw=2)
ax[1].set(xticks=[0,1,2],xticklabels=['Limbic\nn=7','Early neocortical\nn=5','Neocortical\nn=8'],ylabel='External-adapted LPI score (SD)',title='B  Neuropathological stage among PD donors')
main=r[(r.comparison=='PD_vs_control')&r.outcome.eq('LPI_score')].iloc[0];path=r[(r.comparison=='PD_pathology_ordinal')&r.outcome.eq('LPI_score')].iloc[0]
fig.text(.08,.13,f'Age/sex-adjusted PD contrast: {main.beta:.2f} SD (95% CI {main.lo:.2f} to {main.hi:.2f}); P={main.p:.3f}.',fontsize=9)
fig.text(.08,.08,f'Adjusted ordinal-stage estimate: {path.beta:.2f} SD/step (95% CI {path.lo:.2f} to {path.hi:.2f}); q={path.q:.3f}.',fontsize=9)
fig.text(.08,.025,'Points are individual donors; bars are medians. This cross-sectional analysis does not test longitudinal Part I change.',fontsize=9)
for ext in ['png','pdf']:fig.savefig(B/f'delivery_v1/figures/Figure_S6_external_postmortem_CSF.{ext}',dpi=240,bbox_inches='tight')
plt.close(fig)
report={'status':'completed_and_verified','tests':6,'primary_auxiliary':main.to_dict(),'pathology_score':path.to_dict(),'independent_matrix_checks':verification,'source_files_unchanged':all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in manifest.items()),'software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':sm.__version__},'interpretation':'Independent postmortem cross-sectional context, not external longitudinal progression replication'}
# Primary nominal p has no adjusted q because it is the single auxiliary primary test.
report['primary_auxiliary']['q']=None
(O/'COMPLETE.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8');print(r.to_string(index=False))
