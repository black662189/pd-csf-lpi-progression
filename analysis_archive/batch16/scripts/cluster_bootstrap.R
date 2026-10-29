Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
.libPaths(c(Sys.getenv('PD_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lme4);library(jsonlite)})
setDTthreads(1);options(warn=1)
P<-Sys.getenv('PD_PROJECT_ROOT');B<-file.path(P,'batch16');O<-file.path(B,'results/bootstrap_v1');stopifnot(!dir.exists(O));dir.create(O)
datasets<-list(score=fread(file.path(P,'batch9/results/classes_v1/model_data_CSF__LPI.tsv')),PI_conditional=fread(file.path(P,'batch11/results/parent_v1/full_model_data_CSF__LPI.tsv')))
all<-list();summary<-list();pilot<-list()
for(k in seq_along(datasets)){
 context<-names(datasets)[k];d<-datasets[[k]];d[,`:=`(PATNO=as.character(PATNO),sex=factor(sex),genetic_group=factor(genetic_group))]
 old<-readRDS(file.path(B,'results/diagnostics_v1',paste0(context,'__diagonal.rds')));ff<-formula(old);ids<-unique(d$PATNO);blocks<-split(d,d$PATNO)
 run<-function(i){
  set.seed(20261030+k*100000+i);sampled<-sample(ids,length(ids),replace=TRUE)
  dd<-rbindlist(lapply(seq_along(sampled),function(j){v<-copy(blocks[[sampled[j]]]);v[,PATNO:=paste0('b',j)];v}));dd[,PATNO:=factor(PATNO)]
  stopifnot(uniqueN(dd$PATNO)==length(ids),!anyDuplicated(dd[,.(PATNO,time)]))
  warnings<-character();t<-proc.time()[3]
  m<-tryCatch(withCallingHandlers(lmer(ff,data=dd,REML=TRUE,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))),warning=function(e){warnings<<-c(warnings,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(m,'error'))return(data.table(context=context,replicate=i,status='model_error',beta=NA_real_,seconds=unname(proc.time()[3]-t),detail=conditionMessage(m)))
  conv<-m@optinfo$conv$lme4$messages;drop<-attr(getME(m,'X'),'col.dropped');beta<-unname(fixef(m)['time:x'])
  status<-if(length(drop))'rank_deficient'else if(length(conv))'convergence_warning'else if(isSingular(m,tol=1e-4))'singular_fit'else if(!is.finite(beta))'nonfinite'else'tested'
  data.table(context=context,replicate=i,status=status,beta=beta,seconds=unname(proc.time()[3]-t),detail=paste(unique(c(warnings,conv)),collapse='; '))
 }
 first<-rbindlist(parallel::mclapply(1:20,run,mc.cores=8,mc.set.seed=FALSE));pilot[[context]]<-first;fwrite(rbindlist(pilot),file.path(O,'pilot.tsv'),sep='\t');stopifnot(mean(first$status=='tested')>=.95)
 rest<-rbindlist(parallel::mclapply(21:1000,run,mc.cores=8,mc.set.seed=FALSE));ans<-rbind(first,rest);all[[context]]<-ans;good<-ans[status=='tested'];qualified<-nrow(good)>=950
 ci<-if(qualified)quantile(good$beta,c(.025,.975),type=7)else c(NA_real_,NA_real_)
 summary[[context]]<-data.table(context=context,people=length(ids),observations=nrow(d),planned=1000,qualified=nrow(good),failed=1000-nrow(good),CI_qualified=qualified,original_beta=unname(fixef(old)['time:x']),bootstrap_median=median(good$beta),percentile_lo=ci[1],percentile_hi=ci[2],bootstrap_sd=sd(good$beta))
}
fwrite(rbindlist(all),file.path(O,'all_replicates.tsv'),sep='\t');fwrite(rbindlist(summary),file.path(O,'bootstrap_summary.tsv'),sep='\t');capture.output(sessionInfo(),file=file.path(O,'sessionInfo.txt'))
write_json(list(status='completed',planned=2000,workers=8),file.path(O,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE);print(rbindlist(summary))
