#!/usr/bin/env Rscript
# Batch 8: fixed three-hypothesis clinical interpretation family.
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
B<-args[1];O<-args[2]
ROOT<-Sys.getenv('PD_EVIDENCE_ROOT')
.libPaths(c(file.path(ROOT,'Rlib_r453'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lmerTest)})
setDTthreads(1);set.seed(20261028);options(warn=1)
if(dir.exists(O))stop('New versioned output directory required')
dir.create(O,recursive=TRUE)
write_t<-function(d,n)fwrite(d,file.path(O,paste0(n,'.tsv')),sep='\t')
zpath<-file.path(dirname(B),'batch6/results/lpi_v1/lpi_inputs.rds')
z<-readRDS(zpath);d<-copy(z$ds$CSF);f<-z$fs$CSF
stopifnot(nrow(d)==1389L,uniqueN(d$PATNO)==200L,all(d$time>0),!anyDuplicated(d[,.(PATNO,EVENT_ID)]))
stopifnot(identical(levels(d$genetic_group),c('genetic','sporadic')))
fit_one<-function(dd,ff,id,term,role){
 r<-data.table(analysis=id,role=role,term=term,n=nrow(dd),people=uniqueN(dd$PATNO),status='not_tested',beta=NA_real_,se=NA_real_,df=NA_real_,p=NA_real_,lo=NA_real_,hi=NA_real_,detail='',seconds=NA_real_,formula=paste(deparse(ff),collapse=' '))
 warnings<-character();t0<-proc.time()[['elapsed']]
 m<-tryCatch(withCallingHandlers(lmer(ff,data=dd,REML=TRUE,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))),warning=function(e){warnings<<-c(warnings,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
 r$seconds<-proc.time()[['elapsed']]-t0
 if(inherits(m,'error')){r$status<-'model_error';r$detail<-conditionMessage(m);return(list(result=r,model=NULL))}
 co<-coef(summary(m));conv<-m@optinfo$conv$lme4$messages;drop<-attr(getME(m,'X'),'col.dropped')
 if(!term%in%rownames(co)){r$status<-'term_not_estimable';return(list(result=r,model=m))}
 q<-co[term,];crit<-qt(.975,q['df'])
 r[,`:=`(status='tested',beta=unname(q['Estimate']),se=unname(q['Std. Error']),df=unname(q['df']),p=unname(q['Pr(>|t|)']),lo=unname(q['Estimate']-crit*q['Std. Error']),hi=unname(q['Estimate']+crit*q['Std. Error']),detail=paste(unique(c(warnings,conv)),collapse='; '))]
 if(length(drop))r$status<-'rank_deficient'else if(length(conv))r$status<-'convergence_warning'else if(isSingular(m,tol=1e-4))r$status<-'singular_fit'
 if(!all(is.finite(unlist(r[,.(beta,se,df,p,lo,hi)]))))r$status<-'nonfinite_inference'
 if(r$status!='tested')r$p<-NA_real_
 list(result=r,model=m)
}
export_fit<-function(ans,dd,id){
 write_t(ans$result,paste0(id,'_result'));write_t(dd,paste0(id,'_data'))
 if(!is.null(ans$model)){
  saveRDS(ans$model,file.path(O,paste0(id,'_model.rds')))
  co<-as.data.table(coef(summary(ans$model)),keep.rownames='term');write_t(co,paste0(id,'_fixed_effects'))
  vv<-as.matrix(vcov(ans$model));write_t(as.data.table(vv,keep.rownames='term'),paste0(id,'_vcov'))
  write_t(as.data.table(as.data.frame(VarCorr(ans$model))),paste0(id,'_random_variance'))
 }
}
rep<-fit_one(d,f,'CSF_reproduction','time:x','existing_result_reproduction')
ref<-z$reproduction[analysis=='CSF']
rep$result[,`:=`(reference_beta=ref$beta,reference_se=ref$se,reference_p=ref$p,pass=status=='tested'&n==ref$n&people==ref$people&abs(beta-ref$beta)<1e-7&abs(se-ref$se)<1e-7&abs(p-ref$p)<1e-8)]
export_fit(rep,d,'CSF_reproduction');stopifnot(isTRUE(rep$result$pass))
cat('SMOKE_REPRODUCTION_PASS\n')

P<-file.path(B,'results/descriptive_v1')
base<-fread(file.path(P,'subpart_baseline.tsv'),colClasses=c(PATNO='character'))
future<-fread(file.path(P,'subpart_future.tsv'),colClasses=c(PATNO='character'))
stopifnot(nrow(base)==200L,nrow(future)==1379L,uniqueN(future$PATNO)==200L,!anyDuplicated(base$PATNO),!anyDuplicated(future[,.(PATNO,EVENT_ID)]))
anc<-z$anchors$CSF[in_model==TRUE]
bi<-match(base$PATNO,anc$PATNO);stopifnot(!anyNA(bi),all(base$EVENT_ID==anc$EVENT_ID[bi]))
stopifnot(all(base$NP1RTOT+base$NP1PTOT==base$updrs1_score),all(future$NP1RTOT+future$NP1PTOT==future$updrs1_score))
stopifnot(all(base$NP1RTOT>=0&base$NP1RTOT<=24),all(base$NP1PTOT>=0&base$NP1PTOT<=28),all(future$NP1RTOT>=0&future$NP1RTOT<=24),all(future$NP1PTOT>=0&future$NP1PTOT<=28))
parts<-future[,.(PATNO,EVENT_ID,NP1RTOT,NP1PTOT,frozen_total=updrs1_score)]
common<-merge(d,parts,by=c('PATNO','EVENT_ID'),all=FALSE,sort=FALSE)
stopifnot(nrow(common)==1379L,all(common$y==common$frozen_total),all(common[,uniqueN(time),by=PATNO]$V1>=2L))
ix<-match(common$PATNO,base$PATNO);stopifnot(!anyNA(ix),all(common$y0==base$updrs1_score[ix]))
dr<-copy(common);dr[,`:=`(y=NP1RTOT,y0=base$NP1RTOT[ix])]
dp<-copy(common);dp[,`:=`(y=NP1PTOT,y0=base$NP1PTOT[ix])]
group_counts<-unique(d[,.(PATNO,genetic_group)])[,.(people=.N),by=genetic_group]
stopifnot(group_counts[genetic_group=='genetic']$people==126L,group_counts[genetic_group=='sporadic']$people==74L)
write_t(group_counts,'genetic_group_counts')
fg<-update(f,.~.+x:genetic_group+time:x:genetic_group)
jobs<-list(
 list(id='same_sample_total',dd=common,ff=f,term='time:x',role='same_sample_reference'),
 list(id='H1_rater',dd=dr,ff=f,term='time:x',role='H1'),
 list(id='H2_patient',dd=dp,ff=f,term='time:x',role='H2'),
 list(id='H3_genetic_interaction',dd=d,ff=fg,term='time:x:genetic_groupsporadic',role='H3'))
run<-function(j){ans<-fit_one(j$dd,j$ff,j$id,j$term,j$role);export_fit(ans,j$dd,j$id);ans}
answers<-parallel::mclapply(jobs,run,mc.cores=3,mc.preschedule=FALSE)
stopifnot(all(vapply(answers,is.list,logical(1))))
all<-rbindlist(lapply(answers,`[[`,'result'),fill=TRUE)
h<-all[role%in%c('H1','H2','H3')];stopifnot(nrow(h)==3L,setequal(h$role,c('H1','H2','H3')))
h[,`:=`(family_n=3L,q_batch8=p.adjust(p,'BH',n=3L))]
write_t(h,'hypotheses');write_t(all,'all_models')
if(answers[[4]]$result$status=='tested'){
 m<-answers[[4]]$model;bb<-fixef(m);vv<-as.matrix(vcov(m));rows<-list();weights<-list()
 for(g in c('genetic','sporadic')){
  L<-setNames(rep(0,length(bb)),names(bb));L['time:x']<-1
  if(g=='sporadic')L['time:x:genetic_groupsporadic']<-1
  ct<-contest1D(m,L,ddf='Satterthwaite',confint=TRUE,level=.95)
  estimate<-unname(sum(L*bb));se<-sqrt(drop(t(L)%*%vv%*%L));df<-ct[['df']];ci<-estimate+c(-1,1)*qt(.975,df)*se
  stopifnot(abs(estimate-ct[['Estimate']])<1e-10,abs(se-ct[['Std. Error']])<1e-10)
  rows[[g]]<-data.table(group=g,people=group_counts[genetic_group==g]$people,beta=estimate,se=se,df=df,lo=ci[1],hi=ci[2],inference='descriptive_linear_combination_from_one_interaction_model')
  weights[[g]]<-data.table(group=g,term=names(L),weight=unname(L))
 }
 write_t(rbindlist(rows),'genetic_group_contrasts');write_t(rbindlist(weights),'genetic_contrast_weights')
}else{
 write_t(data.table(status='not_estimable',reason=answers[[4]]$result$status),'genetic_contrast_status')
}
write_t(data.table(metric=c('new_models_fitted','existing_models_refitted','planned_hypotheses','tested_hypotheses'),value=c(4,1,3,sum(h$status=='tested'))),'execution_counts')
writeLines(c('Hypotheses fixed before execution: H1 rater, H2 patient, H3 genetic interaction. BH n=3 including unavailable tests.','Exposure scales and original anchors retained. New outcome baseline covariates use their respective subtotals.','No bootstrap, LOO, further subgroup/endpoint search or model fallback performed.','All additions are post-selection internal exploratory analyses; no external validation.'),file.path(O,'analysis_scope.txt'))
capture.output(sessionInfo(),file=file.path(O,'sessionInfo.txt'))
cat('BATCH8_MODELS_FINISHED tested=',sum(h$status=='tested'),' planned=3\n')
