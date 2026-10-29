#!/usr/bin/env Rscript
# Fixed exploratory parent-conditional extension. All old inputs are read-only.
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
P<-Sys.getenv('PD_PROJECT_ROOT');B<-file.path(P,'batch11')
.libPaths(c(Sys.getenv('PD_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lmerTest)})
setDTthreads(1);options(warn=1)
mode<-commandArgs(TRUE)[1];stopifnot(mode%in%c('pilot','full'))
O<-file.path(B,'results/parent_v1');I<-file.path(B,'input/parent_v3')
if(mode=='pilot'){stopifnot(!dir.exists(O));dir.create(O,recursive=TRUE)}else stopifnot(file.exists(file.path(O,'PILOT_PASS.txt')))
wt<-function(x,n)fwrite(x,file.path(O,paste0(n,'.tsv')),sep='\t')
rt<-function(p)fread(p,colClasses=c(PATNO='character',EVENT_ID='character'),na.strings=c('NA',''))
classes<-c('LPI','LPC','LPE');parent_map<-c(LPI='PI',LPC='PC',LPE='PE')
F0<-y~time*(x+y0+up0+age0+sex+duration0+ledd0+genetic_group)+(1+time||PATNO)
F1<-update(F0,.~.+parent+time:parent)
fit<-function(d,label,adjust=TRUE,export=FALSE){
 r<-data.table(analysis=label,people=uniqueN(d$PATNO),n=nrow(d),status='not_tested',beta=NA_real_,se=NA_real_,df=NA_real_,p=NA_real_,lo=NA_real_,hi=NA_real_,seconds=NA_real_,detail='')
 stopifnot(uniqueN(d$PATNO)>=30,nrow(d)>=60,all(d$time>0),!anyDuplicated(d[,.(PATNO,time)]),all(d[,uniqueN(time),by=PATNO]$V1>=2))
 t0<-proc.time()[['elapsed']];w<-character()
 m<-tryCatch(withCallingHandlers(lmer(if(adjust)F1 else F0,data=d,REML=TRUE,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
 r$seconds<-proc.time()[['elapsed']]-t0
 if(inherits(m,'error')){r$status<-'model_error';r$detail<-conditionMessage(m);return(r)}
 co<-coef(summary(m));conv<-m@optinfo$conv$lme4$messages;drop<-attr(getME(m,'X'),'col.dropped')
 if(!'time:x'%in%rownames(co)){r$status<-'term_not_estimable';return(r)}
 a<-co['time:x',];crit<-qt(.975,a['df'])
 r[,`:=`(status='tested',beta=unname(a['Estimate']),se=unname(a['Std. Error']),df=unname(a['df']),p=unname(a['Pr(>|t|)']),lo=unname(a['Estimate']-crit*a['Std. Error']),hi=unname(a['Estimate']+crit*a['Std. Error']),detail=paste(unique(c(w,conv)),collapse='; '))]
 if(length(drop))r$status<-'rank_deficient' else if(length(conv))r$status<-'convergence_warning' else if(isSingular(m,tol=1e-4))r$status<-'singular_fit'
 if(!all(is.finite(unlist(r[,.(beta,se,df,p,lo,hi)]))))r$status<-'nonfinite_inference'
 if(r$status!='tested')r$p<-NA_real_
 if(export){wt(as.data.table(co,keep.rownames='term'),paste0(label,'_coefficients'));saveRDS(m,file.path(O,paste0(label,'_fit.rds')))}
 r
}
if(mode=='pilot'){
 scores<-rt(file.path(I,'parent_scores.tsv'));members<-fread(file.path(B,'protocol/parent_membership_v3.tsv'))
 stopifnot(!anyDuplicated(scores[,.(assay,PATNO)]))
 # Independent score arithmetic: recompute directly from raw unit rows in R.
 checks<-list()
 for(a in c('CSF','plasma')){
  aa<-scores[assay==a];raw<-rt(file.path(B,paste0('input/parent_v2/parent_raw_units_',a,'.tsv')))[UNITS=='adjusted_area_ratio']
  for(pa in unname(parent_map)){
   fs<-members[assay==a&parent==pa]$feature
   z<-sapply(fs,function(f){v<-merge(aa[,.(PATNO,EVENT_ID)],raw[TESTNAME==f,.(PATNO,EVENT_ID,TESTVALUE)],by=c('PATNO','EVENT_ID'),all.x=TRUE,sort=FALSE);v<-v[match(aa$PATNO,v$PATNO)];x<-v$TESTVALUE;x[!is.finite(x)|x<=0]<-NA_real_;x<-log2(x);(x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)})
   v<-rowMeans(z);v<-(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE);expected<-aa[[pa]]
   err<-max(abs(v-expected),na.rm=TRUE);ok<-identical(is.na(v),is.na(expected))&&err<1e-10
   checks[[paste(a,pa)]]<-data.table(assay=a,parent=pa,members=length(fs),anchors=nrow(aa),max_abs_error=err,pass=ok)
  }
 }
 checks<-rbindlist(checks);wt(checks,'score_verification');stopifnot(all(checks$pass))
 full<-paired<-list();repro<-main<-pairbase<-flow<-corr<-list()
 old<-fread(file.path(P,'batch9/results/classes_v1/primary_score_results.tsv'))
 for(a in c('CSF','plasma'))for(cl in classes){
  label<-paste(a,cl,sep='__');pa<-parent_map[[cl]]
  d0<-rt(file.path(P,paste0('batch9/results/classes_v1/model_data_',label,'.tsv')))
  ref<-fit(d0,label,FALSE);prior<-old[analysis==label]
  ref[,pass:=status=='tested'&abs(beta-prior$beta)<1e-7&abs(se-prior$se)<1e-7&n==prior$n&people==prior$people];repro[[label]]<-ref
  ps<-scores[assay==a&is.finite(get(pa)),.(PATNO,parent=get(pa))]
  d<-merge(d0,ps,by='PATNO');setorder(d,PATNO,time);full[[label]]<-d
  wt(d,paste0('full_model_data_',label))
  r0<-fit(d,paste0(label,'_same_rows'),FALSE,TRUE);r0[,model:='same_rows_reference']
  r1<-fit(d,paste0(label,'_conditional'),TRUE,TRUE);r1[,model:='parent_conditional']
  rr<-rbind(r0,r1);rr[,`:=`(assay=a,score=cl,parent_class=pa)];main[[label]]<-rr
  pp<-unique(d[,.(PATNO,x,parent)]);rho<-cor(pp$x,pp$parent)
  corr[[label]]<-data.table(assay=a,score=cl,parent=pa,people=nrow(pp),pearson=rho,two_variable_VIF=1/(1-rho^2),SE_inflation=r1$se/r0$se)
  flow[[label]]<-data.table(assay=a,score=cl,original_people=uniqueN(d0$PATNO),original_rows=nrow(d0),retained_people=uniqueN(d$PATNO),retained_rows=nrow(d),excluded_ids=paste(setdiff(d0$PATNO,d$PATNO),collapse=';'))
  dp<-rt(file.path(P,paste0('batch9/results/classes_v1/paired_model_data_',a,'_',cl,'.tsv')))
  paired[[label]]<-merge(dp,ps,by='PATNO')
 }
 repro<-rbindlist(repro);wt(repro,'baseline_reproduction');stopifnot(all(repro$pass))
 main<-rbindlist(main);main[model=='parent_conditional',q_family6:=p.adjust(p,'BH',n=6)];wt(main,'full_results')
 wt(rbindlist(flow),'sample_flow');wt(rbindlist(corr),'parent_correlations')
 common<-Reduce(intersect,lapply(paired,function(d)unique(d$PATNO)))
 paired_scale<-list()
 for(label in names(paired)){
  d<-paired[[label]][PATNO%in%common];setorder(d,PATNO,time)
  pp<-unique(d[,.(PATNO,parent)]);mu<-mean(pp$parent);ss<-sd(pp$parent);d[,parent:=(parent-mu)/ss]
  paired[[label]]<-d;wt(d,paste0('paired_model_data_',label))
  r<-fit(d,paste0('paired_',label),TRUE,TRUE);r[,`:=`(assay=strsplit(label,'__',fixed=TRUE)[[1]][1],score=strsplit(label,'__',fixed=TRUE)[[1]][2])];pairbase[[label]]<-r
  paired_scale[[label]]<-data.table(analysis=label,people=nrow(pp),parent_mean=mu,parent_sd=ss)
 }
 cv<-setdiff(names(paired[[1]]),c('x','parent'))
 stopifnot(all(vapply(paired,function(d)identical(d[,..cv],paired[[1]][,..cv]),logical(1))))
 pairbase<-rbindlist(pairbase);wt(pairbase,'paired_baseline');wt(rbindlist(paired_scale),'paired_parent_scaling')
 stopifnot(all(pairbase$status=='tested'))
 # One already selected clinical domain, explicitly exploratory.
 ess0<-rt(file.path(P,'batch10/results/clinical_v1/ess_model_data.tsv'))
 ess<-merge(ess0,scores[assay=='CSF'&is.finite(PI),.(PATNO,parent=PI)],by='PATNO')
 e0<-fit(ess,'ESS_same_rows',FALSE,TRUE);e1<-fit(ess,'ESS_conditional',TRUE,TRUE)
 wt(ess,'ESS_model_data');wt(rbind(e0,e1),'ESS_results')
 cache<-list(full=full,paired=paired,paired_people=length(common),paired_rows=nrow(paired[[1]]))
 saveRDS(cache,file.path(O,'parent_inputs.rds'));writeLines(unname(tools::md5sum(file.path(O,'parent_inputs.rds'))),file.path(O,'cache_md5.txt'))
}else{
 stopifnot(!file.exists(file.path(O,'paired_bootstrap.tsv')),identical(readLines(file.path(O,'cache_md5.txt')),unname(tools::md5sum(file.path(O,'parent_inputs.rds')))))
 cache<-readRDS(file.path(O,'parent_inputs.rds'))
}
ids<-sort(unique(cache$paired[[1]]$PATNO));parts<-lapply(cache$paired,function(d)split(d,d$PATNO))
set.seed(20261030);draws<-replicate(2000,sample(ids,length(ids),replace=TRUE),simplify=FALSE)
boot<-function(j){
 draw<-draws[[j]]
 rbindlist(lapply(classes,function(cl){
  rr<-lapply(c('CSF','plasma'),function(a){
   ps<-parts[[paste(a,cl,sep='__')]]
   d<-rbindlist(lapply(seq_along(draw),function(k){v<-copy(ps[[draw[k]]]);v[,PATNO:=paste0('b',k)];v}))
   stopifnot(uniqueN(d$PATNO)==length(ids));fit(d,paste(a,cl))
  })
  data.table(iter=j,score=cl,CSF_beta=rr[[1]]$beta,plasma_beta=rr[[2]]$beta,delta=rr[[1]]$beta-rr[[2]]$beta,qualified=rr[[1]]$status=='tested'&rr[[2]]$status=='tested',CSF_status=rr[[1]]$status,plasma_status=rr[[2]]$status,seconds=rr[[1]]$seconds+rr[[2]]$seconds)
 }))
}
if(mode=='pilot'){
 ans<-rbindlist(parallel::mclapply(1:20,boot,mc.cores=8,mc.preschedule=FALSE));wt(ans,'paired_pilot')
 su<-ans[,.(planned=20,qualified=sum(qualified),failure_fraction=mean(!qualified),model_seconds=sum(seconds)),by=score];wt(su,'paired_pilot_summary')
 stopifnot(all(su$failure_fraction<=.05));writeLines('pass',file.path(O,'PILOT_PASS.txt'));print(su)
}else{
 dir.create(file.path(O,'bootstrap_parts'));pilot<-fread(file.path(O,'paired_pilot.tsv'))
 rest<-parallel::mclapply(21:2000,function(j){r<-boot(j);fwrite(r,file.path(O,'bootstrap_parts',sprintf('rep_%04d.tsv',j)),sep='\t');r},mc.cores=8,mc.preschedule=FALSE)
 ans<-rbindlist(c(list(pilot),rest));setorder(ans,score,iter)
 stopifnot(nrow(ans)==6000,!anyDuplicated(ans[,.(iter,score)]));wt(ans,'paired_bootstrap')
 base<-fread(file.path(O,'paired_baseline.tsv'))
 contrasts<-rbindlist(lapply(classes,function(cl){
  d<-ans[score==cl];frac<-mean(!d$qualified)
  ci<-if(frac<=.05)quantile(d[qualified==TRUE]$delta,c(.025,.975,.05/6,1-.05/6))else rep(NA_real_,4)
  bc<-base[score==cl&assay=='CSF']$beta;bp<-base[score==cl&assay=='plasma']$beta
  data.table(score=cl,people=cache$paired_people,rows=cache$paired_rows,beta_CSF=bc,beta_plasma=bp,delta=bc-bp,lo95=unname(ci[1]),hi95=unname(ci[2]),lo_family98_333=unname(ci[3]),hi_family98_333=unname(ci[4]),bootstrap_planned=2000,bootstrap_qualified=sum(d$qualified),failure_fraction=frac,status=if(frac<=.05)'estimated_exploratory'else 'unqualified_failure_rate')
 }))
 contrasts[,family_supported:=status=='estimated_exploratory'&(hi_family98_333<0|lo_family98_333>0)];wt(contrasts,'paired_contrasts')
 writeLines('All fixed models and paired resampling finished. Inspect qualification before inference.',file.path(O,'COMPLETE.txt'));print(contrasts)
}
capture.output(sessionInfo(),file=file.path(O,paste0(mode,'_sessionInfo.txt')))
