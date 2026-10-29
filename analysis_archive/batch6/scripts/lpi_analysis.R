#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
mode<-args[1];O<-args[2];stopifnot(mode%in%c('prepare','stability','paired_pilot','paired_full'))
ROOT<-Sys.getenv('PD_EVIDENCE_ROOT')
OLD<-file.path(Sys.getenv('PD_PROJECT_ROOT'), 'batch5/results/tg_v1')
.libPaths(c(file.path(ROOT,'Rlib_r453'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lmerTest)})
setDTthreads(1);set.seed(20261028);options(warn=1)
dir.create(O,recursive=TRUE,showWarnings=FALSE)
write_t<-function(d,n)fwrite(d,file.path(O,paste0(n,'.tsv')),sep='\t')
rgz<-function(p)fread(cmd=paste('gzip -dc',shQuote(p)),na.strings=c('NA','','NaN'))
cache<-file.path(O,'lpi_inputs.rds')
fit<-function(d,f,label){
 r<-data.table(analysis=label,n=nrow(d),people=uniqueN(d$PATNO),status='not_tested',beta=NA_real_,se=NA_real_,df=NA_real_,p=NA_real_,lo=NA_real_,hi=NA_real_,detail='',seconds=NA_real_)
 if(nrow(d)<60||uniqueN(d$PATNO)<30||uniqueN(d$y)<3){r$status<-'insufficient_data';return(r)}
 w<-character();t0<-proc.time()[['elapsed']]
 m<-tryCatch(withCallingHandlers(lmer(f,data=d,REML=TRUE,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
 r$seconds<-proc.time()[['elapsed']]-t0
 if(inherits(m,'error')){r$status<-'model_error';r$detail<-conditionMessage(m);return(r)}
 co<-coef(summary(m));conv<-m@optinfo$conv$lme4$messages;drop<-attr(getME(m,'X'),'col.dropped')
 if(!'time:x'%in%rownames(co)){r$status<-'term_not_estimable';return(r)}
 z<-co['time:x',];crit<-qt(.975,z['df'])
 r[,`:=`(status='tested',beta=unname(z['Estimate']),se=unname(z['Std. Error']),df=unname(z['df']),p=unname(z['Pr(>|t|)']),lo=unname(z['Estimate']-crit*z['Std. Error']),hi=unname(z['Estimate']+crit*z['Std. Error']),detail=paste(unique(c(w,conv)),collapse='; '))]
 if(length(drop))r$status<-'rank_deficient'else if(length(conv))r$status<-'convergence_warning'else if(isSingular(m,tol=1e-4))r$status<-'singular_fit'
 if(!all(is.finite(unlist(r[,.(beta,se,df,p,lo,hi)]))))r$status<-'nonfinite_inference'
 if(r$status!='tested')r$p<-NA_real_
 r
}
if(mode=='prepare'){
 if(file.exists(cache))stop('New output directory required')
 inp<-readRDS(file.path(OLD,'tg_frozen_inputs.rds'));ref<-fread(file.path(OLD,'tg_full547_UPDRS1_randomslope.tsv'))
 stopifnot(identical(readLines(file.path(OLD,'tg_smoke_pass_cache_md5.txt')),unname(tools::md5sum(file.path(OLD,'tg_frozen_inputs.rds')))))
 ids<-c(CSF=438L,plasma=166L);stopifnot(all(inp$slots[slot_id%in%ids]$feature=='LPI(16:0)'))
 ds<-lapply(ids,function(i)copy(inp$datasets[[i]]));fs<-lapply(ids,function(i)inp$formulas[[i]])
 repro<-rbindlist(lapply(names(ids),function(a){r<-fit(ds[[a]],fs[[a]],a);q<-ref[slot_id==ids[[a]]];r[,`:=`(reference_beta=q$beta,reference_se=q$se,reference_p=q$p,original_q_full547=q$q_full547,pass=status=='tested'&n==q$n&people==q$people&abs(beta-q$beta)<1e-7&abs(se-q$se)<1e-7&abs(p-q$p)<1e-8)];r}))
 write_t(repro,'reproduction');stopifnot(all(repro$pass))
 cl<-rgz(file.path(ROOT,'ppmi_clinical_analysis_20261023/cleaned/clinical_model_ready.tsv.gz'))
 cl[,PATNO:=as.character(PATNO)];cl<-cl[source_curated==1&is.finite(visit_month)]
 long<-rgz(file.path(ROOT,'metabolomics_v2/cleaned/PPMI_named_adjusted_all_visits_cleaned_long.tsv.gz'))
 long[,`:=`(PATNO=as.character(PATNO),EVENT_ID=CLINICAL_EVENT,assay=ifelse(TYPE=='Plasma','plasma','CSF'))]
 stopifnot(!anyDuplicated(cl[,.(PATNO,EVENT_ID)]),!anyDuplicated(long[,.(PATNO,EVENT_ID,assay,TESTNAME)]))
 lt<-long[TESTNAME=='LPI(16:0)'];write_t(lt[,.(records=.N,people=uniqueN(PATNO),missing=sum(!is.finite(TESTVALUE)),nonpositive=sum(TESTVALUE<=0,na.rm=TRUE),log2_disagreement=sum(is.finite(TESTVALUE)&TESTVALUE>0&abs(log2(TESTVALUE)-log2_value)>1e-9,na.rm=TRUE)),by=.(assay,UNITS)],'measurement_values')
 stopifnot(all(lt$UNITS=='adjusted_area_ratio'),all(is.finite(lt$log2_value)),all(lt$TESTVALUE>0),all(abs(log2(lt$TESTVALUE)-lt$log2_value)<1e-9))
 finite<-function(d,vs)Reduce(`&`,lapply(vs,function(v)if(is.numeric(d[[v]]))is.finite(d[[v]])else !is.na(d[[v]])))
 anchors<-list();flow<-list()
 for(a in names(ids)){
  sm<-unique(long[assay==a,.(PATNO,EVENT_ID)])
  aa<-merge(cl,sm,by=c('PATNO','EVENT_ID'))
  aa<-aa[pd_eligible==1&disease_years>=0&hy_main>=1&hy_main<=5&finite(aa,c('age10','sex','duration5','logledd','genetic_group'))]
  setorder(aa,PATNO,visit_month,EVENT_ID);aa<-aa[!duplicated(PATNO)]
  aa<-merge(aa,lt[assay==a,.(PATNO,EVENT_ID,log2_value,RUNDATE)],by=c('PATNO','EVENT_ID'),all.x=TRUE)
  s<-inp$scales[slot_id==ids[[a]]];aa[,x_frozen:=(log2_value-s$log2_mean)/s$log2_sd]
  d<-ds[[a]];ix<-match(d$PATNO,aa$PATNO)
  stopifnot(!anyNA(ix),all(abs(d$x-aa$x_frozen[ix])<1e-9),all(d$index_month==aa$visit_month[ix]),all(d$time>0),!anyDuplicated(d[,.(PATNO,time)]))
  stopifnot(all(d[,uniqueN(time),by=PATNO]$V1>=2))
  aa[,in_model:=PATNO%in%d$PATNO];anchors[[a]]<-aa
  flow[[a]]<-data.table(assay=a,clinical_anchor_people=nrow(aa),detected_anchor_people=sum(is.finite(aa$log2_value)),model_people=uniqueN(d$PATNO),model_rows=nrow(d),BMI_available_model=sum(is.finite(aa[in_model==TRUE]$BMI)))
  write_t(aa[,.(PATNO,EVENT_ID,visit_month,log2_value,x_frozen,BMI,RUNDATE,in_model)],paste0('anchors_',a))
 }
 write_t(rbindlist(flow),'sample_flow')
 pd<-merge(anchors$CSF[in_model==TRUE,.(PATNO,anchor_c=EVENT_ID,month_c=visit_month)],anchors$plasma[in_model==TRUE,.(PATNO,anchor_p=EVENT_ID,month_p=visit_month)],by='PATNO')
 pd[,same_anchor:=anchor_c==anchor_p&month_c==month_p];write_t(pd,'paired_person_eligibility')
 people<-pd[same_anchor==TRUE]$PATNO
 dc<-copy(ds$CSF[PATNO%in%people]);dp<-copy(ds$plasma[PATNO%in%people])
 keys<-c('PATNO','EVENT_ID','time');common<-merge(dc[,..keys],dp[,..keys],by=keys)
 dc<-merge(dc,common,by=keys);dp<-merge(dp,common,by=keys)
 keep<-dc[,.(visits=.N),by=PATNO][visits>=2]$PATNO;dc<-dc[PATNO%in%keep];dp<-dp[PATNO%in%keep]
 setorderv(dc,keys);setorderv(dp,keys)
 cv<-setdiff(names(dc),'x');stopifnot(identical(dc[,..cv],dp[,..cv]))
 pairscale<-list()
 if(nrow(dc))for(a in names(ids)){
  d<-if(a=='CSF')dc else dp
  u<-unique(d[,.(PATNO,x)]);mu<-mean(u$x);ss<-sd(u$x);stopifnot(is.finite(ss),ss>0)
  d[,x:=(x-mu)/ss];pairscale[[a]]<-data.table(assay=a,people=nrow(u),mean_original_z=mu,sd_original_z=ss)
  if(a=='CSF')dc<-d else dp<-d
 }
 write_t(rbindlist(pairscale),'paired_scales')
 pf<-data.table(stage=c('CSF_model','plasma_model','both_models','same_original_anchor','common_future_eligible'),people=c(uniqueN(ds$CSF$PATNO),uniqueN(ds$plasma$PATNO),nrow(pd),sum(pd$same_anchor),uniqueN(dc$PATNO)),rows=c(nrow(ds$CSF),nrow(ds$plasma),NA,NA,nrow(dc)))
 write_t(pf,'paired_sample_flow')
 saveRDS(list(ds=ds,fs=fs,anchors=anchors,paired=list(CSF=dc,plasma=dp),reproduction=repro,pair_n=uniqueN(dc$PATNO)),cache)
 cat('PREPARE_VERIFIED paired_people=',uniqueN(dc$PATNO),'paired_rows=',nrow(dc),'\n')
}else{
 z<-readRDS(cache)
 if(mode=='stability'){
  stopifnot(!file.exists(file.path(O,'lpi_LOO.tsv')))
  d<-z$ds$CSF;f<-z$fs$CSF;ids<-sort(unique(d$PATNO))
  ans<-rbindlist(parallel::mclapply(ids,function(id){r<-fit(d[PATNO!=id],f,paste0('omit_',id));r[,omitted_person:=id];r},mc.cores=8,mc.preschedule=FALSE))
  stopifnot(nrow(ans)==length(ids),!anyDuplicated(ans$omitted_person))
  write_t(ans,'lpi_LOO');good<-ans[status=='tested'];bb<-z$reproduction[analysis=='CSF']$beta
  write_t(data.table(planned=length(ids),tested=nrow(good),failed=nrow(ans)-nrow(good),baseline_beta=bb,min_beta=if(nrow(good))min(good$beta)else NA_real_,max_beta=if(nrow(good))max(good$beta)else NA_real_,direction_reversals=sum(sign(good$beta)!=sign(bb))),'lpi_LOO_summary')
  a<-z$anchors$CSF;d[,BMI0:=a$BMI[match(PATNO,a$PATNO)]]
  persons<-unique(d[,.(PATNO,BMI0)]);coverage<-mean(is.finite(persons$BMI0))
  stopifnot(!any(is.finite(persons$BMI0)&persons$BMI0<=0))
  if(coverage>=.8){dd<-d[is.finite(BMI0)];bmi<-rbindlist(list(fit(dd,f,'BMI_complete_original'),fit(dd,update(f,.~.+BMI0+time:BMI0),'BMI_complete_adjusted')));write_t(bmi,'BMI_sensitivity')}
  write_t(data.table(model_people=nrow(persons),BMI_available=sum(is.finite(persons$BMI0)),coverage=coverage,status=if(coverage>=.8)'executed'else 'not_qualified'),'BMI_eligibility')
  write_t(a[in_model==TRUE,.(people=.N),by=RUNDATE],'RUNDATE_description')
  writeLines('RUNDATE is not authenticated as an experimental batch. No batch-adjusted model fitted; no arbitrary date binning. Baseline LEDD already adjusted, residual medication confounding remains possible.',file.path(O,'batch_medication_scope.txt'))
  cat('STABILITY_FINISHED\n')
 }else{
  if(z$pair_n<50){write_t(data.table(status='not_qualified',people=z$pair_n,reason='fewer_than_50_same_original_anchor_paired_people'),'paired_status');cat('PAIRED_CLOSED_INSUFFICIENT_N\n');quit(save='no')}
  pc<-z$paired$CSF;pp<-z$paired$plasma
  baseline<-rbindlist(list(fit(pc,z$fs$CSF,'CSF'),fit(pp,z$fs$plasma,'plasma')))
  write_t(baseline,paste0(mode,'_baseline'))
  if(any(baseline$status!='tested')){write_t(data.table(status='not_qualified',people=z$pair_n,reason='paired_baseline_model_failed'),'paired_status');cat('PAIRED_CLOSED_MODEL\n');quit(save='no')}
  ids<-sort(unique(pc$PATNO));cparts<-split(pc,pc$PATNO);pparts<-split(pp,pp$PATNO)
  set.seed(20261028);draws<-replicate(1000,sample(ids,length(ids),replace=TRUE),simplify=FALSE)
  boot<-function(j){draw<-draws[[j]]
   dc<-rbindlist(lapply(seq_along(draw),function(k){d<-copy(cparts[[draw[k]]]);d[,PATNO:=paste0('b',k)];d}))
   dp<-rbindlist(lapply(seq_along(draw),function(k){d<-copy(pparts[[draw[k]]]);d[,PATNO:=paste0('b',k)];d}))
   stopifnot(uniqueN(dc$PATNO)==length(ids),identical(dc[,.(PATNO,EVENT_ID,time,y)],dp[,.(PATNO,EVENT_ID,time,y)]))
   rc<-fit(dc,z$fs$CSF,'CSF');rp<-fit(dp,z$fs$plasma,'plasma')
   data.table(iter=j,CSF_status=rc$status,plasma_status=rp$status,CSF_beta=rc$beta,plasma_beta=rp$beta,delta=rc$beta-rp$beta,qualified=rc$status=='tested'&rp$status=='tested',seconds=rc$seconds+rp$seconds,CSF_detail=rc$detail,plasma_detail=rp$detail)
  }
  if(mode=='paired_pilot'){
   ans<-rbindlist(parallel::mclapply(1:20,boot,mc.cores=8));write_t(ans,'paired_pilot')
   write_t(data.table(planned=20,qualified=sum(ans$qualified),failure_fraction=mean(!ans$qualified),model_seconds=sum(ans$seconds),projected_1000_model_seconds=sum(ans$seconds)*50),'paired_pilot_summary')
   if(mean(!ans$qualified)>.05){write_t(data.table(status='pilot_failed',people=z$pair_n,reason='pilot_model_failure_fraction_over_5pct'),'paired_status');cat('PAIRED_PILOT_NOT_QUALIFIED\n')}else writeLines('pass',file.path(O,'paired_pilot_pass.txt'))
  }else{
   if(!file.exists(file.path(O,'paired_pilot_pass.txt'))){cat('PAIRED_FULL_NOT_RUN_NO_PILOT_PASS\n');quit(save='no')}
   stopifnot(!file.exists(file.path(O,'paired_bootstrap.tsv')))
   pilot<-fread(file.path(O,'paired_pilot.tsv'));ans<-rbindlist(c(list(pilot),parallel::mclapply(21:1000,boot,mc.cores=8,mc.preschedule=FALSE)))
   stopifnot(nrow(ans)==1000,!anyDuplicated(ans$iter));setorder(ans,iter);write_t(ans,'paired_bootstrap')
   fraction<-mean(!ans$qualified);ci<-if(fraction<=.05)quantile(ans[qualified==TRUE]$delta,c(.025,.975))else c(NA_real_,NA_real_)
   write_t(data.table(status=if(fraction<=.05)'estimated_exploratory'else 'unqualified_failure_rate',people=z$pair_n,rows=nrow(pc),planned=1000,qualified=sum(ans$qualified),failure_fraction=fraction,beta_CSF=baseline[analysis=='CSF']$beta,beta_plasma=baseline[analysis=='plasma']$beta,delta=baseline[analysis=='CSF']$beta-baseline[analysis=='plasma']$beta,lo=unname(ci[1]),hi=unname(ci[2])),'paired_contrast')
   v<-merge(unique(pc[,.(PATNO,CSF=x)]),unique(pp[,.(PATNO,plasma=x)]),by='PATNO');write_t(v,'paired_exposures')
   write_t(data.table(people=nrow(v),pearson=cor(v$CSF,v$plasma),spearman=cor(v$CSF,v$plasma,method='spearman'),interpretation='descriptive_not_transport_or_origin'),'paired_correlation')
   cat('PAIRED_FULL_FINISHED qualified=',sum(ans$qualified),'\n')
  }
 }
}
capture.output(sessionInfo(),file=file.path(O,paste0(mode,'_sessionInfo.txt')))
