import os
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(os.environ['PD_PROJECT_ROOT']);B=P/'batch16';R=B/'server200/results/diagnostics_v1';C=B/'server200/results/contrasts_v2'
O=B/'delivery_v1';(O/'figures').mkdir(parents=True,exist_ok=True);(O/'tables').mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':240})
navy='#244C66';teal='#138D90';orange='#CB7040';gray='#67727D'
d=pd.read_csv(P/'batch9/results/classes_v1/model_data_CSF__LPI.tsv',sep='\t')
persons=d.drop_duplicates('PATNO').copy();persons['tertile']=pd.qcut(persons.x,3,labels=['Lower','Middle','Higher'])
assert len(persons)==200 and persons.PATNO.nunique()==200
long=d.merge(persons[['PATNO','tertile']],on='PATNO',validate='many_to_one')
long['year']=np.ceil(long.time-1e-9).astype(int)
annual=long.groupby(['PATNO','tertile','year'],observed=True).y.mean().reset_index()
baseline=persons[['PATNO','tertile','y0']].rename(columns={'y0':'y'}).assign(year=0)
annual=pd.concat([baseline,annual],ignore_index=True)
group=annual.groupby(['tertile','year'],observed=True).y.agg(['count','mean','std']).reset_index();group['se']=group['std']/np.sqrt(group['count']);group['halfwidth']=stats.t.ppf(.975,group['count']-1)*group.se
group.to_csv(O/'tables/descriptive_trajectory_summary.csv',index=False)
c=pd.read_csv(C/'clinical_contrasts.tsv',sep='\t');assert c[['Estimate','Std. Error','lower','upper','df']].notna().all().all()
c.to_csv(O/'tables/Table_S31_clinical_contrasts.csv',index=False)
fig,ax=plt.subplots(2,2,figsize=(11.5,9.3));fig.subplots_adjust(hspace=.43,wspace=.32,bottom=.16,top=.94)
ax[0,0].hist(persons.x,bins=18,color=navy,edgecolor='white',alpha=.9)
for q in persons.x.quantile([1/3,2/3]):ax[0,0].axvline(q,color=orange,ls='--',lw=1.3)
ax[0,0].set(xlabel='Baseline CSF LPI score\n(frozen eligible-index SD scale)',ylabel='Participants',title='A  Exposure distribution in 200 participants')
for name,col in zip(['Lower','Middle','Higher'],[orange,gray,teal]):
 g=group[(group.tertile==name)&group.year.le(10)&group['count'].ge(10)]
 ax[0,1].plot(g.year,g['mean'],'o-',color=col,label=name,lw=1.6,ms=3.5)
 ax[0,1].fill_between(g.year,g['mean']-g.halfwidth,g['mean']+g.halfwidth,color=col,alpha=.11)
ax[0,1].set(xlabel='Years after baseline sampling\n(annual observation intervals)',ylabel='Observed mean MDS-UPDRS Part I',title='B  Descriptive trajectories by baseline tertile',xticks=[0,1,3,5,8,10]);ax[0,1].legend(title='LPI score',frameon=False,fontsize=9)
for panel,context,title in [(ax[1,0],'score','C  LPI score'),(ax[1,1],'PI_conditional','D  LPI score conditional on measured PI')]:
 for structure,col,offset,label in [('diagonal',navy,-.12,'Original'),('correlated',teal,0,'Correlated random effects'),('quadratic',orange,.12,'Quadratic time')]:
  g=c[(c.model==context+'__'+structure)&(c.contrast=='change_difference_per_1SD')].sort_values('year')
  panel.errorbar(g.year+offset,g.Estimate,yerr=[g.Estimate-g.lower,g.upper-g.Estimate],fmt='o-',color=col,ms=4,lw=1.4,capsize=3,label=label)
 panel.axhline(0,color=gray,ls='--',lw=1);panel.set(xlabel='Years after baseline sampling',ylabel='Difference in Part I change per 1 SD\n(model estimate and 95% CI)',title=title,xticks=[1,3,5],ylim=(-2.9,.35))
ax[1,0].legend(frameon=False,fontsize=8.2,loc='lower left')
counts=[]
for name in ['Lower','Middle','Higher']:
 v=group[(group.tertile==name)&group.year.isin([0,1,3,5,8,10])].set_index('year')['count']
 counts.append(f"{name}: "+'/'.join(str(int(v.loc[k])) for k in [0,1,3,5,8,10]))
fig.text(.07,.095,'B: observed participant counts at years 0/1/3/5/8/10; '+'; '.join(counts),fontsize=8.2)
fig.text(.07,.063,'B uses one participant-level mean per annual interval; shading is a descriptive 95% interval. Observation sets vary over time.',fontsize=8.2)
fig.text(.07,.035,'C–D use the continuous score, unchanged participants and covariates. Contrasts describe associations, not treatment effects.',fontsize=8.2)
for ext in ['png','pdf']:fig.savefig(O/f'figures/Figure_6_exposure_trajectories_and_contrasts.{ext}',bbox_inches='tight')
plt.close(fig)

res=pd.read_csv(R/'residuals_PRIVATE.tsv',sep='\t')
fig,ax=plt.subplots(2,3,figsize=(12,7.5));fig.subplots_adjust(hspace=.44,wspace=.32,bottom=.12,top=.93)
summary=[]
for row,(m,lab) in enumerate([('score__diagonal','LPI score'),('PI_conditional__diagonal','LPI score plus PI')]):
 g=res[res.model.eq(m)];rr=g.scaled_residual.to_numpy()
 for col,x,xlab in [(0,g.fitted_conditional,'Conditional fitted Part I'),(1,g.time,'Years after sampling')]:
  a=ax[row,col];a.scatter(x,rr,s=8,color=navy,alpha=.24,rasterized=True);a.axhline(0,color=orange,lw=1)
  bins=pd.qcut(x,10,duplicates='drop');b=pd.DataFrame({'x':x.to_numpy(),'r':rr,'bin':bins}).groupby('bin',observed=True)[['x','r']].mean()
  a.plot(b.x,b.r,'o-',color=teal,lw=1.4,ms=3);a.set(xlabel=xlab,ylabel='Residual / residual SD',title=f'{lab}: residual pattern')
 (theoretical,observed),(slope,intercept,r)=stats.probplot(rr,dist='norm')
 a=ax[row,2];a.scatter(theoretical,observed,s=9,color=navy,alpha=.55,rasterized=True);a.plot(theoretical,intercept+slope*theoretical,color=orange,lw=1);a.set(xlabel='Theoretical normal quantiles',ylabel='Ordered scaled residuals',title=f'{lab}: normal Q–Q')
 summary.append({'model':m,'n':len(g),'max_scaled_residual':float(rr.max()),'min_scaled_residual':float(rr.min()),'fraction_abs_over3':float((abs(rr)>3).mean()),'outcome_zero_fraction':float(g.y.eq(0).mean())})
fig.text(.07,.035,'Original REML fits. Green lines show descriptive bin means. Residual upper tails depart from the normal reference; no observations were removed.',fontsize=9)
for ext in ['png','pdf']:fig.savefig(O/f'figures/Figure_S5_model_diagnostics.{ext}',bbox_inches='tight')
plt.close(fig)
pd.DataFrame(summary).to_csv(O/'tables/residual_summary.csv',index=False)
st=pd.read_csv(R/'model_status.tsv',sep='\t');co=pd.read_csv(R/'fixed_effects.tsv',sep='\t');metrics=pd.read_csv(R/'model_metrics.tsv',sep='\t');lrt=pd.read_csv(R/'likelihood_comparisons.tsv',sep='\t')
co.merge(st[['model','status','people','observations','REML']],on='model',validate='many_to_one').to_csv(O/'tables/Table_S29a_model_coefficients.csv',index=False)
st.merge(metrics,on=['model','REML'],validate='one_to_one').to_csv(O/'tables/Table_S29b_model_fit_and_status.csv',index=False)
lrt.to_csv(O/'tables/Table_S30_model_comparisons.csv',index=False)
record={'status':'generated_from_verified_model_outputs','models':len(st),'qualified':int(st.status.eq('tested').sum()),'participants':200,'rows':1389,'source_contrasts':'contrasts_v2; v1 contrast table superseded','figures':[p.name for p in (O/'figures').glob('*')],'source_hashes':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [C/'clinical_contrasts.tsv',R/'fixed_effects.tsv',R/'model_metrics.tsv',P/'batch9/results/classes_v1/model_data_CSF__LPI.tsv']}}
(O/'figure_and_table_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf-8');print(json.dumps(record,indent=2))
