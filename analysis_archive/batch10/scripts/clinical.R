#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
P<-Sys.getenv('PD_PROJECT_ROOT');B<-file.path(P,'batch10')
.libPaths(c(Sys.getenv('PD_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lmerTest)})
setDTthreads(1);options(warn=1);set.seed(20261029)
mode<-commandArgs(TRUE)[1];stopifnot(mode%in%c('background_corrected','clinical'))
O<-file.path(B,'results',if(mode=='background_corrected')'background_v2' else 'clinical_v1');stopifnot(!dir.exists(O));dir.create(O,recursive=TRUE)
wt<-function(d,n)fwrite(d,file.path(O,paste0(n,'.tsv')),sep='\t')
readt<-function(p)fread(p,colClasses=c(PATNO='character',EVENT_ID='character'),na.strings=c('NA',''))
F<-y~time*(x+y0+up0+age0+sex+duration0+ledd0+genetic_group)+(1+time||PATNO)
fit<-function(d,label,formula=F){
 r<-data.table(analysis=label,people=uniqueN(d$PATNO),n=nrow(d),status='not_tested',beta=NA_real_,se=NA_real_,df=NA_real_,p=NA_real_,lo=NA_real_,hi=NA_real_,detail='')
 if(uniqueN(d$PATNO)<30||nrow(d)<60||uniqueN(d$y)<3){r$status<-'insufficient_data';return(r)}
 stopifnot(!anyDuplicated(d[,.(PATNO,time)]),all(d$time>0),all(d[,uniqueN(time),by=PATNO]$V1>=2))
 w<-character();m<-tryCatch(withCallingHandlers(lmer(formula,data=d,REML=TRUE,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(m,'error')){r$status<-'model_error';r$detail<-conditionMessage(m);return(r)}
 co<-coef(summary(m));conv<-m@optinfo$conv$lme4$messages;drop<-attr(getME(m,'X'),'col.dropped');stopifnot('time:x'%in%rownames(co))
 a<-co['time:x',];crit<-qt(.975,a['df'])
 r[,`:=`(status='tested',beta=unname(a['Estimate']),se=unname(a['Std. Error']),df=unname(a['df']),p=unname(a['Pr(>|t|)']),lo=unname(a['Estimate']-crit*a['Std. Error']),hi=unname(a['Estimate']+crit*a['Std. Error']),detail=paste(unique(c(w,conv)),collapse='; '))]
 if(length(drop))r$status<-'rank_deficient' else if(length(conv))r$status<-'convergence_warning' else if(isSingular(m,tol=1e-4))r$status<-'singular_fit'
 if(!all(is.finite(unlist(r[,.(beta,se,df,p,lo,hi)]))))r$status<-'nonfinite_inference'
 if(r$status!='tested')r$p<-NA_real_
 wt(as.data.table(co,keep.rownames='term'),paste0(label,'_all_coefficients'))
 saveRDS(m,file.path(O,paste0(label,'_fit.rds')));r
}
base<-readt(file.path(P,'batch9/results/classes_v1/model_data_CSF__LPI.tsv'))
stopifnot(nrow(base)==1389,uniqueN(base$PATNO)==200)
if(mode=='background_corrected'){
 ref<-fit(base,'frozen_reference');prior<-fread(file.path(P,'batch9/results/classes_v1/primary_score_results.tsv'))[analysis=='CSF__LPI']
 ref[,pass:=status=='tested'&abs(beta-prior$beta)<1e-7&abs(se-prior$se)<1e-7&n==prior$n&people==prior$people]
 wt(ref,'baseline_reproduction');stopifnot(isTRUE(ref$pass))
 bg<-readt(file.path(B,'input/background_scores_v2.tsv'));stopifnot(!anyDuplicated(bg$PATNO))
 d<-merge(base,bg[,.(PATNO,background)],by='PATNO');d<-d[is.finite(background)]
 wt(d,'background_model_data')
 same<-fit(d,'same_rows_reference');adj<-fit(d,'background_adjusted',update(F,.~.+background+time:background))
 wt(rbind(same,adj),'background_results')
}else{
 stopifnot(file.exists(file.path(B,'results/background_v2/COMPLETE.txt')))
 fut<-readt(file.path(B,'input/clinical_future.tsv'));anc<-readt(file.path(B,'input/clinical_anchors.tsv'))
 stopifnot(nrow(fut)==nrow(base),!anyDuplicated(fut[,.(PATNO,EVENT_ID)]),!anyDuplicated(anc$PATNO))
 ranges<-list(moca=c(0,30),gds=c(0,15),scopa=c(0,69),ess=c(0,24))
 result<-flows<-list()
 for(k in names(ranges)){
   lo<-ranges[[k]][1];hi<-ranges[[k]][2];s<-if(k=='moca')-1 else 1
   d<-merge(fut,anc[,.(PATNO,endpoint0=get(paste0(k,'0')))],by='PATNO')
   d[,`:=`(y=get(k),y0=endpoint0)]
   ok<-is.finite(d$y)&is.finite(d$y0)&d$y>=lo&d$y<=hi&d$y0>=lo&d$y0<=hi&d$y%%1==0&d$y0%%1==0
   d<-d[ok];before_n<-nrow(d);before_people<-uniqueN(d$PATNO)
   keep<-d[,.N,by=PATNO][N>=2]$PATNO;d<-d[PATNO%in%keep];setorder(d,PATNO,time)
   flows[[k]]<-data.table(endpoint=k,original_people=200,original_rows=1389,baseline_future_valid_rows=before_n,baseline_future_valid_people=before_people,final_people=uniqueN(d$PATNO),final_rows=nrow(d))
   r<-fit(d,k);r[,`:=`(endpoint=k,beta_worse=s*beta,lo_worse=if(s==1)lo else -hi,hi_worse=if(s==1)hi else -lo)];result[[k]]<-r
   wt(d[,.(PATNO,EVENT_ID,time,x,y,y0,up0,age0,sex,duration0,ledd0,genetic_group)],paste0(k,'_model_data'))
 }
 rr<-rbindlist(result);rr[,q_family4:=p.adjust(p,'BH',n=4)];wt(rr,'clinical_results');wt(rbindlist(flows),'sample_flow')
 d<-copy(fut[is.finite(ledd_visit)&ledd_visit>=0]);keep<-d[,.N,by=PATNO][N>=2]$PATNO;d<-d[PATNO%in%keep]
 d[,ledd_delta:=log1p(ledd_visit)-ledd0]
 dynamic<-fit(d,'time_varying_LEDD',update(F,.~.+ledd_delta));wt(d,'time_varying_LEDD_model_data')
 short<-base[time<=5];keep<-short[,.N,by=PATNO][N>=2]$PATNO;short<-short[PATNO%in%keep]
 window<-fit(short,'five_year');wt(short,'five_year_model_data')
 bg<-fread(file.path(B,'results/background_v2/background_results.tsv'))[analysis=='background_adjusted']
 sens<-rbindlist(list(bg,dynamic,window),fill=TRUE);sens[,q_family3:=p.adjust(p,'BH',n=3)];wt(sens,'all_sensitivity_results')
}
capture.output(sessionInfo(),file=file.path(O,'sessionInfo.txt'))
writeLines('Fixed module completed; inspect per-model status before inference.',file.path(O,'COMPLETE.txt'))
