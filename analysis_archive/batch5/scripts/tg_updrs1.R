#!/usr/bin/env Rscript
# Batch 5: existing-data exploratory UPDRS I full-family analysis.
# Exact anchor/scale/outcome logic from 30_clinical_extensions.R; no model fallback.
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3L)stop('Usage: tg_updrs1.R MODE ROOT OUTPUT [WORKERS]; MODE=prepare|smoke|full|loo')
mode<-args[1];ROOT<-args[2];O<-args[3];workers<-if(length(args)>=4L)as.integer(args[4])else 1L
stopifnot(mode%in%c('prepare','smoke','full','loo'),is.finite(workers),workers>=1L,workers<=8L)
.libPaths(c(file.path(ROOT,'Rlib_r453'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lmerTest)})
setDTthreads(1);set.seed(20261028);options(warn=1)
dir.create(O,recursive=TRUE,showWarnings=FALSE)
B<-file.path(ROOT,'ppmi_clinical_analysis_20261023');M<-file.path(ROOT,'metabolomics_v2')
N<-file.path(ROOT,'public_followup_20261025')
rgz<-function(p)fread(cmd=paste('gzip -dc',shQuote(p)),na.strings=c('NA','','NaN'))
finite<-function(d,vs)Reduce(`&`,lapply(vs,function(v)if(is.numeric(d[[v]]))is.finite(d[[v]])else !is.na(d[[v]])))
firstrow<-function(d){d<-copy(d);setorder(d,PATNO,visit_month,EVENT_ID);d[!duplicated(PATNO)]}
eligible<-function(d)d[pd_eligible==1&disease_years>=0&hy_main>=1&hy_main<=5&finite(d,c('age10','sex','duration5','logledd','genetic_group'))]
safe<-function(expr){w<-character();v<-tryCatch(withCallingHandlers(expr,warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e);list(value=v,warnings=unique(w))}
make_formula<-function(d){
 cv<-c('y0','up0','age0','sex','duration0','ledd0','genetic_group')
 cv<-cv[vapply(cv,function(v)uniqueN(d[[v]][!is.na(d[[v]])])>1,logical(1))]
 as.formula(paste('y~time*(',paste(c('x',cv),collapse='+'),')+(1+time||PATNO)'))
}
cachefile<-file.path(O,'tg_frozen_inputs.rds')
if(mode=='prepare'){
 if(file.exists(cachefile))stop('Input cache already exists; use a new OUTPUT directory to reprepare.')
 files<-c(clinical=file.path(B,'cleaned/clinical_model_ready.tsv.gz'),raw=file.path(B,'cleaned/curated_original_snapshot.tsv.gz'),metabolites=file.path(M,'cleaned/PPMI_named_adjusted_all_visits_cleaned_long.tsv.gz'),family=file.path(B,'results/all_model_results.tsv.gz'),reference=file.path(N,'results/30_randomslope_results.tsv'))
 stopifnot(all(file.exists(files)))
 manifest<-data.table(role=names(files),path=unname(files),bytes=file.info(files)$size,md5=unname(tools::md5sum(files)))
 fwrite(manifest,file.path(O,'tg_input_manifest.tsv'),sep='\t')
 prior<-rgz(files['family']);slots<-prior[family=='Metabolite'&endpoint=='UPDRS_future_slope'&model=='primary',.(assay,feature)]
 stopifnot(nrow(slots)==547L,!anyDuplicated(slots),sum(slots$assay=='CSF')==249L,sum(slots$assay=='plasma')==298L)
 slots[,slot_id:=seq_len(.N)];setcolorder(slots,c('slot_id','assay','feature'));fwrite(slots,file.path(O,'tg_frozen_547_slots.tsv'),sep='\t')
 clinical<-rgz(files['clinical']);raw<-rgz(files['raw']);clinical[,PATNO:=as.character(PATNO)];raw[,PATNO:=as.character(PATNO)]
 stopifnot(!anyDuplicated(raw[,.(PATNO,EVENT_ID)]),!anyDuplicated(clinical[,.(PATNO,EVENT_ID)]))
 clinical<-merge(clinical,raw[,.(PATNO,EVENT_ID,updrs1_score)],by=c('PATNO','EVENT_ID'),all.x=TRUE)
 clinical<-clinical[source_curated==1&is.finite(visit_month)]
 clinical[,sex:=factor(sex,levels=c(0,1),labels=c('Female','Male'))];clinical[,genetic_group:=factor(genetic_group)]
 bad<-is.finite(clinical$updrs1_score)&(clinical$updrs1_score<0|clinical$updrs1_score>52)
 fwrite(data.table(endpoint='updrs1_score',out_of_range=sum(bad)),file.path(O,'tg_endpoint_quality.tsv'),sep='\t')
 set(clinical,which(bad),'updrs1_score',NA_real_)
 long<-rgz(files['metabolites']);long[,`:=`(PATNO=as.character(PATNO),EVENT_ID=CLINICAL_EVENT,assay=ifelse(TYPE=='Plasma','plasma','CSF'))]
 stopifnot(!anyDuplicated(long[,.(PATNO,EVENT_ID,assay,TESTNAME)]))
 anchors<-lapply(c('plasma','CSF'),function(a){sm<-unique(long[assay==a,.(PATNO,EVENT_ID)]);firstrow(eligible(merge(clinical,sm,by=c('PATNO','EVENT_ID'))))});names(anchors)<-c('plasma','CSF')
 fut<-clinical[pd_eligible==1&is.finite(updrs1_score)&is.finite(visit_month),.(PATNO,EVENT_ID,visit_month,y=updrs1_score)]
 datasets<-vector('list',547L);flows<-vector('list',547L);forms<-vector('list',547L);scale_info<-vector('list',547L)
 for(i in seq_len(nrow(slots))){
  as0<-slots$assay[i];ft<-slots$feature[i]
  d<-merge(anchors[[as0]],long[assay==as0&TESTNAME==ft,.(PATNO,EVENT_ID,x=log2_value,RUNDATE)],by=c('PATNO','EVENT_ID'),all.x=TRUE)
  d<-d[is.finite(x)];xm<-mean(d$x);xs<-sd(d$x)
  scale_info[[i]]<-data.table(slot_id=i,log2_mean=xm,log2_sd=xs,scale_n=nrow(d))
  # EXACT ancestral scale: after frozen clinical anchor and feature availability,
  # before endpoint completeness or future follow-up eligibility. Never rescale later.
  d[,x:=(x-xm)/xs]
  dd<-d[is.finite(updrs1_score)]
  aa<-dd[is.finite(updrs_main),.(PATNO,index_month=visit_month,x,y0=updrs1_score,up0=updrs_main,hy0=hy_main,age0=age10,duration0=duration5,ledd0=logledd,sex,genetic_group)]
  jj<-merge(fut,aa,by='PATNO',allow.cartesian=TRUE);jj<-jj[visit_month>index_month];jj[,time:=(visit_month-index_month)/12]
  setorder(jj,PATNO,time,EVENT_ID);jj<-jj[!duplicated(jj[,.(PATNO,time)])];jj<-jj[,if(.N>=2).SD,by=PATNO]
  if(nrow(jj))stopifnot(all(jj$time>0),!anyDuplicated(jj[,.(PATNO,time)]))
  datasets[[i]]<-jj;forms[[i]]<-make_formula(jj)
  flows[[i]]<-data.table(slot_id=i,assay_anchor_people=nrow(anchors[[as0]]),detected_anchor_people=nrow(d),endpoint_complete_people=nrow(dd),baseline_UPDRS3_complete_people=nrow(aa),future_people=uniqueN(jj$PATNO),future_rows=nrow(jj),median_max_followup_years=if(nrow(jj))median(jj[,max(time),by=PATNO]$V1)else NA_real_)
  if(as0=='plasma'&&ft=='TG(20:4_34:3)'){
   fwrite(d,file.path(O,'tg_target_exposure_anchors.tsv'),sep='\t');fwrite(jj,file.path(O,'tg_target_future_inputs.tsv'),sep='\t')
   fu<-jj[,.(future_visits=.N,first_future_year=min(time),last_future_year=max(time)),by=PATNO]
   retained<-merge(d[,.(PATNO,EVENT_ID,x,updrs1_score,updrs_main,BMI,RUNDATE)],fu,by='PATNO',all.x=TRUE);retained[,retained_future:=!is.na(future_visits)]
   fwrite(retained,file.path(O,'tg_target_followup_availability.tsv'),sep='\t')
   av<-rbindlist(lapply(c('BMI','RUNDATE','logledd','treatment'),function(v)data.table(field=v,anchor_people=nrow(d),available=sum(!is.na(d[[v]])&as.character(d[[v]])!=''),source='existing frozen exposure anchor',interpretation=if(v=='RUNDATE')'Assay run date; not automatically a biological collection batch'else 'Availability only; no additional adjustment fitted')))
   av<-rbind(av,data.table(field=c('fasting_status','lipid_lowering_medication'),anchor_people=nrow(d),available=NA_integer_,source='not a validated field in clinical_model_ready or cleaned metabolite long table',interpretation='Not inferred from absent fields; broader source qualification required'))
   fwrite(av,file.path(O,'tg_covariate_availability.tsv'),sep='\t')
  }
 }
 refs<-fread(files['reference']);refs<-refs[family=='Selected_metabolite_extensions'&endpoint=='updrs1_score'&design=='future_slope']
 saveRDS(list(slots=slots,datasets=datasets,formulas=forms,flows=rbindlist(flows),scales=rbindlist(scale_info),references=refs,manifest=manifest),cachefile,compress=TRUE)
 fwrite(merge(slots,rbindlist(flows),by='slot_id'),file.path(O,'tg_sample_flow.tsv'),sep='\t')
 fwrite(merge(slots,rbindlist(scale_info),by='slot_id'),file.path(O,'tg_anchor_scale.tsv'),sep='\t')
 writeLines(c('Outcome: UPDRS I; exploratory after endpoint inspection.','547 original fluid-feature slots retained including failures. BH n=547.','Fixed anchor before feature and endpoint availability. Exposure log2 z-score fixed before outcome filtering.','Strictly future deduplicated visits; at least two/person. Uncorrelated random intercept/time slope; no fallback.','LOO freezes exposure scale and full target formula; it is influence diagnosis, not independent replication.'),file.path(O,'tg_analysis_lock.txt'))
 capture.output(sessionInfo(),file=file.path(O,'tg_prepare_sessionInfo.txt'));cat('TG_PREPARE_FINISHED slots=547\n');quit(save='no')
}
stopifnot(file.exists(cachefile));inp<-readRDS(cachefile)
fit_one<-function(i,omit=NULL){
 d<-copy(inp$datasets[[i]]);if(!is.null(omit))d<-d[PATNO!=omit]
 r<-copy(inp$slots[i]);r[,`:=`(endpoint='updrs1_score',term='time:x',model='uncorrelated_random_intercept_time_slope',n=nrow(d),people=uniqueN(d$PATNO),omitted_person=if(is.null(omit))''else as.character(omit),status='not_tested',beta=NA_real_,se=NA_real_,df=NA_real_,p=NA_real_,p_diagnostic=NA_real_,lo=NA_real_,hi=NA_real_,formula=paste(deparse(inp$formulas[[i]]),collapse=' '),detail='',elapsed_seconds=NA_real_)]
 t0<-proc.time()[['elapsed']]
 if(nrow(d)<60L||uniqueN(d$PATNO)<30L||uniqueN(d$y)<3L){r$status<-'insufficient_longitudinal_n_or_variation';return(r)}
 if(!all(is.finite(d$x))||!is.finite(sd(d$x))||sd(d$x)==0){r$status<-'nonfinite_or_constant_exposure';return(r)}
 f<-safe(lmer(inp$formulas[[i]],data=d,REML=TRUE,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))))
 r$elapsed_seconds<-proc.time()[['elapsed']]-t0
 if(inherits(f$value,'error')){r$status<-'model_error';r$detail<-conditionMessage(f$value);return(r)}
 fit<-f$value;co<-coef(summary(fit));conv<-fit@optinfo$conv$lme4$messages;drop<-attr(getME(fit,'X'),'col.dropped')
 if(!'time:x'%in%rownames(co)){r$status<-'term_not_estimable';return(r)}
 z<-co['time:x',];crit<-qt(.975,z['df'])
 r[,`:=`(status='tested',beta=unname(z['Estimate']),se=unname(z['Std. Error']),df=unname(z['df']),p=unname(z['Pr(>|t|)']),p_diagnostic=unname(z['Pr(>|t|)']),lo=unname(z['Estimate']-crit*z['Std. Error']),hi=unname(z['Estimate']+crit*z['Std. Error']),detail=paste(c(f$warnings,conv),collapse='; '))]
 if(length(drop))r$status<-'rank_deficient'else if(length(conv))r$status<-'convergence_warning'else if(isSingular(fit,tol=1e-4))r$status<-'singular_fit'
 if(!all(is.finite(unlist(r[,.(beta,se,df,p,lo,hi)])))&&r$status=='tested')r$status<-'nonfinite_coefficient_inference'
 if(r$status!='tested')r$p<-NA_real_
 r
}
target<-which(inp$slots$assay=='plasma'&inp$slots$feature=='TG(20:4_34:3)');stopifnot(length(target)==1L)
if(mode=='smoke'){
 ids<-c(target,which(inp$slots$assay=='CSF'&inp$slots$feature=='Ornithine'));stopifnot(length(ids)==2L)
 ans<-rbindlist(lapply(ids,fit_one));ref<-inp$references[,.(assay,feature,reference_status=status,reference_n=n,reference_people=people,reference_beta=beta,reference_p=p,reference_se=se)]
 z<-merge(ans,ref,by=c('assay','feature'),all.x=TRUE)
 z[,pass:=status=='tested'&reference_status=='tested'&n==reference_n&people==reference_people&abs(beta-reference_beta)<1e-5&abs(se-reference_se)<1e-5&abs(p-reference_p)<1e-6]
 fwrite(z,file.path(O,'tg_smoke_reproduction.tsv'),sep='\t');print(z[,.(assay,feature,n,people,beta,p,pass)])
 if(!all(z$pass))stop('Smoke reproduction failed; full run is prohibited until the discrepancy is resolved.')
 writeLines(unname(tools::md5sum(cachefile)),file.path(O,'tg_smoke_pass_cache_md5.txt'))
}else{
 gate<-file.path(O,'tg_smoke_pass_cache_md5.txt');stopifnot(file.exists(gate),identical(readLines(gate),unname(tools::md5sum(cachefile))))
 if(.Platform$OS.type!='unix'&&workers>1L)stop('Parallel execution uses Linux fork; set WORKERS=1 on other systems.')
 if(mode=='full'){
  parts<-file.path(O,'full_parts');dir.create(parts,showWarnings=FALSE)
  run<-function(i){r<-fit_one(i);fwrite(r,file.path(parts,sprintf('slot_%03d.tsv',i)),sep='\t');r}
  ans<-rbindlist(parallel::mclapply(seq_len(547L),run,mc.cores=workers,mc.preschedule=FALSE),fill=TRUE)
  stopifnot(nrow(ans)==547L,!anyDuplicated(ans$slot_id),setequal(ans$slot_id,seq_len(547L)))
  setorder(ans,slot_id);ans[,`:=`(planned_family_n=547L,q_full547=p.adjust(p,'BH',n=547L))]
  fwrite(ans,file.path(O,'tg_full547_UPDRS1_randomslope.tsv'),sep='\t')
  fwrite(ans[,.(planned=.N,tested=sum(status=='tested'),nominal=sum(p<.05,na.rm=TRUE),FDR=sum(q_full547<.05,na.rm=TRUE)),by=assay],file.path(O,'tg_full547_summary.tsv'),sep='\t')
  cat('TG_FULL_FINISHED',nrow(ans),'tested',sum(ans$status=='tested'),'FDR',sum(ans$q_full547<.05,na.rm=TRUE),'\n')
 }else{
  baseline<-fit_one(target);stopifnot(baseline$status=='tested')
  ids<-sort(unique(inp$datasets[[target]]$PATNO));parts<-file.path(O,'loo_parts');dir.create(parts,showWarnings=FALSE)
  run<-function(j){r<-fit_one(target,ids[j]);fwrite(r,file.path(parts,sprintf('loo_%03d.tsv',j)),sep='\t');r}
  ans<-rbindlist(parallel::mclapply(seq_along(ids),run,mc.cores=workers,mc.preschedule=FALSE),fill=TRUE)
  stopifnot(nrow(ans)==length(ids),!anyDuplicated(ans$omitted_person))
  ans[,`:=`(baseline_beta=baseline$beta,delta_beta=beta-baseline$beta)];ans[,relative_change:=delta_beta/abs(baseline_beta)]
  fwrite(ans,file.path(O,'tg_target_LOO.tsv'),sep='\t');fwrite(baseline,file.path(O,'tg_target_LOO_baseline.tsv'),sep='\t')
  good<-ans[status=='tested'];summary<-data.table(planned=length(ids),tested=nrow(good),failed=sum(ans$status!='tested'),baseline_beta=baseline$beta,min_beta=if(nrow(good))min(good$beta)else NA_real_,max_beta=if(nrow(good))max(good$beta)else NA_real_,direction_reversals=sum(sign(good$beta)!=sign(baseline$beta)),max_absolute_delta=if(nrow(good))max(abs(good$delta_beta))else NA_real_)
  fwrite(summary,file.path(O,'tg_target_LOO_summary.tsv'),sep='\t');cat('TG_LOO_FINISHED',nrow(ans),'tested',nrow(good),'\n')
 }
}
capture.output(sessionInfo(),file=file.path(O,paste0('tg_',mode,'_sessionInfo.txt')))
