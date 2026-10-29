#!/usr/bin/env Rscript
# Prespecified single-covariate BMI sensitivity on identical complete-case data.
# No exposure restandardization, no anchor advancement, no RUNDATE adjustment.
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop('Usage: tg_bmi_sensitivity.R ROOT TG_RUN_DIRECTORY NEW_OUTPUT_DIRECTORY')
ROOT<-args[1];I<-args[2];O<-args[3]
.libPaths(c(file.path(ROOT,'Rlib_r453'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lmerTest)})
setDTthreads(1);set.seed(20261028);options(warn=1)
if(dir.exists(O)&&length(list.files(O)))stop('Output must be a new or empty directory; prior results are not overwritten.')
dir.create(O,recursive=TRUE,showWarnings=FALSE)
files<-c(cache=file.path(I,'tg_frozen_inputs.rds'),anchor=file.path(I,'tg_target_exposure_anchors.tsv'),full=file.path(I,'tg_full547_UPDRS1_randomslope.tsv'),gate=file.path(I,'tg_smoke_pass_cache_md5.txt'))
stopifnot(all(file.exists(files)),identical(readLines(files['gate']),unname(tools::md5sum(files['cache']))))
fwrite(data.table(role=names(files),path=unname(files),bytes=file.info(files)$size,md5=unname(tools::md5sum(files))),file.path(O,'tg_bmi_input_manifest.tsv'),sep='\t')
inp<-readRDS(files['cache']);target<-which(inp$slots$assay=='plasma'&inp$slots$feature=='TG(20:4_34:3)');stopifnot(length(target)==1L)
d0<-copy(inp$datasets[[target]]);d0[,PATNO:=as.character(PATNO)]
a<-fread(files['anchor'],na.strings=c('NA','','NaN'));a[,PATNO:=as.character(PATNO)]
stopifnot(!anyDuplicated(a$PATNO),all(c('BMI','RUNDATE','x')%in%names(a)),is.numeric(a$BMI),all(d0$PATNO%in%a$PATNO))
full<-fread(files['full']);reference<-full[assay=='plasma'&feature=='TG(20:4_34:3)']
stopifnot(nrow(reference)==1L,reference$status=='tested',reference$n==nrow(d0),reference$people==uniqueN(d0$PATNO))
index<-match(d0$PATNO,a$PATNO)
# Preserve all original row order and covariate/factor encodings; add one column.
stopifnot(all(abs(d0$x-a$x[index])<1e-10))
d0[,BMI0:=a$BMI[index]]
if(any(is.finite(d0$BMI0)&d0$BMI0<=0))stop('Nonpositive finite BMI requires input clarification; no automatic exclusion.')
d<-d0[is.finite(BMI0)]
f0<-inp$formulas[[target]]
f1<-update(f0,.~.+BMI0+time:BMI0)
writeLines(c('One planned BMI sensitivity; frozen anchor, exposure scale and original formula retained.',
 'Original-model and BMI-adjusted fits use identical finite-baseline-BMI participants and rows.',
 'BMI0 uses recorded kg/m2 value; no rescaling, imputation, outlier exclusion, covariate search or random-effect fallback.',
 'Full-sample original estimate is read from the completed 547-family output, not refitted.',
 'Sensitivity P values are descriptive and do not replace original q_full547 or constitute a new discovery family.',
 paste('Base formula:',paste(deparse(f0),collapse=' ')),paste('BMI formula:',paste(deparse(f1),collapse=' '))),file.path(O,'tg_bmi_analysis_lock.txt'))
safe<-function(expr){w<-character();v<-tryCatch(withCallingHandlers(expr,warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e);list(value=v,warnings=unique(w))}
fit<-function(data,formula,label){
 r<-data.table(analysis=label,source='new_fit_same_BMI_complete_case_sample',n=nrow(data),people=uniqueN(data$PATNO),status='not_tested',beta=NA_real_,se=NA_real_,df=NA_real_,p=NA_real_,p_diagnostic=NA_real_,lo=NA_real_,hi=NA_real_,formula=paste(deparse(formula),collapse=' '),detail='',elapsed_seconds=NA_real_)
 if(nrow(data)<60L||uniqueN(data$PATNO)<30L||uniqueN(data$y)<3L){r$status<-'insufficient_longitudinal_n_or_variation';return(r)}
 if(!all(is.finite(data$x))||sd(data$x)==0){r$status<-'nonfinite_or_constant_exposure';return(r)}
 if(label=='BMI_adjusted_complete_cases'&&uniqueN(data$BMI0)<2L){r$status<-'constant_BMI';return(r)}
 t0<-proc.time()[['elapsed']]
 f<-safe(lmer(formula,data=data,REML=TRUE,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))))
 r$elapsed_seconds<-proc.time()[['elapsed']]-t0
 if(inherits(f$value,'error')){r$status<-'model_error';r$detail<-conditionMessage(f$value);return(r)}
 m<-f$value;co<-coef(summary(m));conv<-m@optinfo$conv$lme4$messages;drop<-attr(getME(m,'X'),'col.dropped')
 if(!'time:x'%in%rownames(co)){r$status<-'term_not_estimable';return(r)}
 z<-co['time:x',];crit<-qt(.975,z['df'])
 r[,`:=`(status='tested',beta=unname(z['Estimate']),se=unname(z['Std. Error']),df=unname(z['df']),p=unname(z['Pr(>|t|)']),p_diagnostic=unname(z['Pr(>|t|)']),lo=unname(z['Estimate']-crit*z['Std. Error']),hi=unname(z['Estimate']+crit*z['Std. Error']),detail=paste(c(f$warnings,conv),collapse='; '))]
 if(length(drop))r$status<-'rank_deficient'else if(length(conv))r$status<-'convergence_warning'else if(isSingular(m,tol=1e-4))r$status<-'singular_fit'
 if(!all(is.finite(unlist(r[,.(beta,se,df,p,lo,hi)])))&&r$status=='tested')r$status<-'nonfinite_coefficient_inference'
 if(r$status!='tested')r$p<-NA_real_
 r
}
base<-fit(d,f0,'original_model_BMI_complete_cases')
adjusted<-fit(d,f1,'BMI_adjusted_complete_cases')
baseline<-reference[,.(analysis='original_model_full_sample',source='completed_full547_output_not_refitted',n,people,status,beta,se,df,p,p_diagnostic,lo,hi,formula,detail,elapsed_seconds)]
ans<-rbindlist(list(baseline,base,adjusted),fill=TRUE)
ans[,`:=`(original_discovery_q_full547=reference$q_full547,interpretation='sensitivity_not_new_FDR_discovery')]
fwrite(ans,file.path(O,'tg_bmi_model_comparison.tsv'),sep='\t')
delta<-data.table(contrast=c('BMI_missingness_subset_effect','BMI_adjustment_same_sample_effect'),from_analysis=c(baseline$analysis,base$analysis),to_analysis=c(base$analysis,adjusted$analysis),from_status=c(baseline$status,base$status),to_status=c(base$status,adjusted$status),people_lost=c(baseline$people-base$people,base$people-adjusted$people),rows_lost=c(baseline$n-base$n,base$n-adjusted$n),delta_beta=c(base$beta-baseline$beta,adjusted$beta-base$beta),reference_beta=c(baseline$beta,base$beta))
delta[,relative_change:=delta_beta/abs(reference_beta)]
delta[,comparison_qualified:=from_status=='tested'&to_status=='tested']
fwrite(delta,file.path(O,'tg_bmi_effect_decomposition.tsv'),sep='\t')
persons<-unique(d0[,.(PATNO,x,y0,up0,BMI0)]);persons[,BMI_complete:=is.finite(BMI0)]
fwrite(persons,file.path(O,'tg_bmi_person_eligibility.tsv'),sep='\t')
flow<-data.table(stage=c('all_target_exposure_anchors','original_future_model','BMI_complete_future_model','BMI_incomplete_future_model'),people=c(nrow(a),uniqueN(d0$PATNO),uniqueN(d$PATNO),uniqueN(d0[!is.finite(BMI0)]$PATNO)),rows=c(nrow(a),nrow(d0),nrow(d),nrow(d0[!is.finite(BMI0)])),BMI_available_people=c(sum(is.finite(a$BMI)),sum(persons$BMI_complete),uniqueN(d$PATNO),0L))
fwrite(flow,file.path(O,'tg_bmi_sample_flow.tsv'),sep='\t')
fwrite(persons[,.(people=.N,mean_x=mean(x),mean_baseline_UPDRS1=mean(y0),mean_baseline_UPDRS3=mean(up0),BMI_min=if(any(is.finite(BMI0)))min(BMI0[is.finite(BMI0)])else NA_real_,BMI_max=if(any(is.finite(BMI0)))max(BMI0[is.finite(BMI0)])else NA_real_),by=BMI_complete],file.path(O,'tg_bmi_missingness_description.tsv'),sep='\t')
# RUNDATE is an assay date, not proven biological collection batch. Description only.
dates<-a[,.(PATNO,RUNDATE)];dates[,`:=`(in_original_future=PATNO%in%d0$PATNO,in_BMI_complete_future=PATNO%in%d$PATNO)]
fwrite(dates[,.(anchor_people=.N,original_future_people=sum(in_original_future),BMI_complete_future_people=sum(in_BMI_complete_future)),by=RUNDATE],file.path(O,'tg_RUNDATE_levels_description.tsv'),sep='\t')
capture.output(sessionInfo(),file=file.path(O,'tg_bmi_sessionInfo.txt'))
print(ans[,.(analysis,n,people,status,beta,lo,hi,p)]);print(delta)
cat('TG_BMI_SENSITIVITY_FINISHED',sum(ans$status=='tested'),'of',nrow(ans),'model_records\n')
