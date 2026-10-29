#!/usr/bin/env Rscript
# Bounded post-selection diagnostics. Frozen participants, scales and observations.
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
.libPaths(c(Sys.getenv('PD_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lmerTest);library(jsonlite)})
setDTthreads(1);set.seed(20261029);options(warn=1)
P<-Sys.getenv('PD_PROJECT_ROOT'); B<-file.path(P,'batch16'); O<-file.path(B,'results/diagnostics_v1')
stopifnot(!dir.exists(O));dir.create(O)
wt<-function(d,n)fwrite(d,file.path(O,paste0(n,'.tsv')),sep='\t')
read_data<-function(f){d<-fread(file.path(P,f));d[,PATNO:=factor(PATNO)];d[,sex:=factor(sex)];d[,genetic_group:=factor(genetic_group)];d[,time_sq:=(time/5)^2];setorder(d,PATNO,time);stopifnot(!anyDuplicated(d[,.(PATNO,time)]),all(d$time>0),all(d[,uniqueN(time),by=PATNO]$V1>=2));d}
data<-list(score=read_data('batch9/results/classes_v1/model_data_CSF__LPI.tsv'),PI_conditional=read_data('batch11/results/parent_v1/full_model_data_CSF__LPI.tsv'))
stopifnot(nrow(data$score)==1389,uniqueN(data$score$PATNO)==200,nrow(data$PI_conditional)==1379,uniqueN(data$PI_conditional$PATNO)==199)
covars<-'x+y0+up0+age0+sex+duration0+ledd0+genetic_group'
formula_for<-function(context,structure){cv<-paste0(covars,if(context=='PI_conditional')'+parent' else ''); ti<-if(structure=='quadratic')'(time+time_sq)'else'time';re<-if(structure=='correlated')'(1+time|PATNO)'else'(1+time||PATNO)';as.formula(paste('y~',ti,'*(',cv,')+',re))}
fits<-list();statuses<-list();coefficients<-list();metrics<-list();resids<-list();contrasts<-list();comparisons<-list()
fit<-function(d,f,label,reml=TRUE){
 warnings<-character();t0<-proc.time()[3]
 m<-tryCatch(withCallingHandlers(lmer(f,data=d,REML=reml,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))),warning=function(e){warnings<<-c(warnings,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
 if(inherits(m,'error')){statuses[[label]]<<-data.table(model=label,status='model_error',detail=conditionMessage(m));return(NULL)}
 co<-coef(summary(m));conv<-m@optinfo$conv$lme4$messages;dropped<-attr(getME(m,'X'),'col.dropped')
 status<-if(length(dropped))'rank_deficient'else if(length(conv))'convergence_warning'else if(isSingular(m,tol=1e-4))'singular_fit'else'tested'
 statuses[[label]]<<-data.table(model=label,status=status,REML=reml,people=uniqueN(d$PATNO),observations=nrow(d),seconds=unname(proc.time()[3]-t0),formula=paste(deparse(f),collapse=' '),detail=paste(unique(c(warnings,conv)),collapse='; '))
 co<-as.data.table(co,keep.rownames='term');co[,`:=`(model=label,lo=Estimate-qt(.975,df)*`Std. Error`,hi=Estimate+qt(.975,df)*`Std. Error`)];coefficients[[label]]<<-co
 vc<-as.data.table(as.data.frame(VarCorr(m)));vf<-var(predict(m,re.form=NA));vr<-0
 for(i in seq_len(nrow(vc))){a<-vc[i];if(a$grp=='Residual')next;z1<-if(a$var1=='(Intercept)')rep(1,nrow(d))else d[[a$var1]];if(is.na(a$var2))vr<-vr+a$vcov*mean(z1^2)else{z2<-if(a$var2=='(Intercept)')rep(1,nrow(d))else d[[a$var2]];vr<-vr+2*a$vcov*mean(z1*z2)}}
 ve<-sigma(m)^2;metrics[[label]]<<-data.table(model=label,REML=reml,AIC=AIC(m),BIC=BIC(m),logLik=as.numeric(logLik(m)),parameters=attr(logLik(m),'df'),fixed_variance=vf,mean_random_variance=vr,residual_variance=ve,marginal_R2=vf/(vf+vr+ve),conditional_R2=(vf+vr)/(vf+vr+ve))
 vc[,model:=label];wt(vc,paste0(label,'_variance_components'));saveRDS(m,file.path(O,paste0(label,'.rds')));fits[[label]]<<-m
 m
}
design<-function(m,d){model.matrix(delete.response(terms(lme4::nobars(formula(m)))),data=d,contrasts.arg=attr(getME(m,'X'),'contrasts'))[,names(fixef(m)),drop=FALSE]}
get_contrast<-function(m,L,label,kind,h){a<-as.data.table(contest1D(m,L,confint=TRUE));a[,`:=`(model=label,contrast=kind,year=h)];a}
for(context in names(data)){
 d<-data[[context]]
 for(structure in c('diagonal','correlated','quadratic')){
  label<-paste(context,structure,sep='__');m<-fit(d,formula_for(context,structure),label)
  if(structure=='diagonal'){
   stopifnot(!is.null(m),statuses[[label]]$status=='tested')
   ref<-if(context=='score')fread(file.path(P,'batch9/results/classes_v1/primary_score_results.tsv'))[analysis=='CSF__LPI']else fread(file.path(P,'batch11/results/parent_v1/full_results.tsv'))[analysis=='CSF__LPI_conditional']
   stopifnot(nrow(ref)==1,abs(fixef(m)['time:x']-ref$beta)<1e-7,abs(coef(summary(m))['time:x','Std. Error']-ref$se)<1e-7)
   if(context=='score')writeLines('Frozen score model reproduced; representative fit passed.',file.path(O,'PILOT_PASS.txt'))
  }
  if(is.null(m)||statuses[[label]]$status!='tested')next
  resids[[label]]<-data.table(model=label,PATNO=as.character(d$PATNO),time=d$time,y=d$y,x=d$x,fitted_conditional=fitted(m),fitted_population=predict(m,re.form=NA),residual=residuals(m),scaled_residual=residuals(m)/sigma(m))
  pp<-d[!duplicated(PATNO)]
  for(h in c(1,3,5)){
   grid<-function(t,exposure_value=NULL){g<-copy(pp);g[,`:=`(time=t,time_sq=(t/5)^2)];if(!is.null(exposure_value))data.table::set(g,j='x',value=rep(exposure_value,nrow(g)));g}
   lx<-(colMeans(design(m,grid(h,.5)))-colMeans(design(m,grid(h,-.5))))-(colMeans(design(m,grid(0,.5)))-colMeans(design(m,grid(0,-.5))))
   lc<-colMeans(design(m,grid(h)))-colMeans(design(m,grid(0)))
   contrasts[[paste(label,h,'x')]]<-get_contrast(m,lx,label,'change_difference_per_1SD',h)
   contrasts[[paste(label,h,'cohort')]]<-get_contrast(m,lc,label,'cohort_standardized_change',h)
  }
 }
 for(structure in c('diagonal','quadratic'))fit(d,formula_for(context,structure),paste(context,structure,'ML',sep='__'),FALSE)
 ma<-fits[[paste(context,'diagonal','ML',sep='__')]];mb<-fits[[paste(context,'quadratic','ML',sep='__')]]
 if(!is.null(ma)&&!is.null(mb)){a<-as.data.table(anova(ma,mb,refit=FALSE),keep.rownames='row');a[,comparison:=paste0(context,'_quadratic_vs_linear_ML')];comparisons[[context]]<-a}
}
# Same 199 participants/1,379 visits for all incremental-fit comparisons.
d<-data$PI_conditional
clinical<-y~time*(y0+up0+age0+sex+duration0+ledd0+genetic_group)+(1+time||PATNO)
forms<-list(clinical=clinical,clinical_LPI=update(clinical,.~.+x+time:x),clinical_PI=update(clinical,.~.+parent+time:parent),clinical_PI_LPI=update(clinical,.~.+parent+time:parent+x+time:x),clinical_LPI_no_slope=update(clinical,.~.+x),clinical_PI_LPI_no_slope=update(clinical,.~.+parent+time:parent+x))
for(n in names(forms))fit(d,forms[[n]],paste0('increment__',n),FALSE)
pairs<-list(LPI_block=c('clinical','clinical_LPI'),LPI_block_given_PI=c('clinical_PI','clinical_PI_LPI'),LPI_slope=c('clinical_LPI_no_slope','clinical_LPI'),LPI_slope_given_PI=c('clinical_PI_LPI_no_slope','clinical_PI_LPI'))
for(n in names(pairs)){ids<-paste0('increment__',pairs[[n]]);if(all(vapply(ids,function(i)!is.null(fits[[i]])&&statuses[[i]]$status=='tested',logical(1)))){a<-as.data.table(anova(fits[[ids[1]]],fits[[ids[2]]],refit=FALSE),keep.rownames='row');a[,`:=`(row=ids,comparison=n)];comparisons[[n]]<-a}}
wt(rbindlist(statuses,fill=TRUE),'model_status');wt(rbindlist(coefficients,fill=TRUE),'fixed_effects');wt(rbindlist(metrics,fill=TRUE),'model_metrics');wt(rbindlist(comparisons,fill=TRUE),'likelihood_comparisons');wt(rbindlist(contrasts,fill=TRUE),'clinical_contrasts');wt(rbindlist(resids,fill=TRUE),'residuals_PRIVATE')
capture.output(sessionInfo(),file=file.path(O,'sessionInfo.txt'))
complete<-list(status='completed',models=length(statuses),qualified=sum(vapply(statuses,function(z)z$status=='tested',logical(1))),seed=20261029,threads=1,selection='Post-selection internal diagnostics; no independent validation or new discovery claim')
write_json(complete,file.path(O,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE);print(complete)
