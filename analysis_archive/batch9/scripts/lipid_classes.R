#!/usr/bin/env Rscript
# Existing-data exploratory extension; class definitions are fixed in protocol before fitting.
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2, args[1] %in% c('prepare','pilot','full'))
mode <- args[1]; O <- args[2]
ROOT <- Sys.getenv('PD_EVIDENCE_ROOT')
P <- Sys.getenv('PD_PROJECT_ROOT')
.libPaths(c(file.path(ROOT,'Rlib_r453'), .libPaths()))
suppressPackageStartupMessages({library(data.table); library(lmerTest)})
setDTthreads(1); options(warn=1); set.seed(20261029)
dir.create(O,recursive=TRUE,showWarnings=FALSE)
write_t <- function(d,n) fwrite(d,file.path(O,paste0(n,'.tsv')),sep='\t')
cache <- file.path(O,'class_inputs.rds')
classes <- c('LPI','LPC','LPE','LYSO_ALL')
F <- y ~ time*(x+y0+up0+age0+sex+duration0+ledd0+genetic_group)+(1+time||PATNO)
fit <- function(d,label){
 r <- data.table(analysis=label,people=uniqueN(d$PATNO),n=nrow(d),status='not_tested',beta=NA_real_,se=NA_real_,df=NA_real_,p=NA_real_,lo=NA_real_,hi=NA_real_,seconds=NA_real_,detail='')
 if(uniqueN(d$PATNO)<30 || nrow(d)<60 || uniqueN(d$y)<3){r$status<-'insufficient_data';return(r)}
 stopifnot(all(is.finite(d$x)),all(d$time>0),!anyDuplicated(d[,.(PATNO,time)]))
 t0<-proc.time()[['elapsed']];w<-character()
 m<-tryCatch(withCallingHandlers(lmer(F,data=d,REML=TRUE,na.action=na.fail,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=100000))),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')}),error=function(e)e)
 r$seconds<-proc.time()[['elapsed']]-t0
 if(inherits(m,'error')){r$status<-'model_error';r$detail<-conditionMessage(m);return(r)}
 co<-coef(summary(m));conv<-m@optinfo$conv$lme4$messages;drop<-attr(getME(m,'X'),'col.dropped')
 if(!'time:x'%in%rownames(co)){r$status<-'term_not_estimable';return(r)}
 a<-co['time:x',];crit<-qt(.975,a['df'])
 r[,`:=`(status='tested',beta=unname(a['Estimate']),se=unname(a['Std. Error']),df=unname(a['df']),p=unname(a['Pr(>|t|)']),lo=unname(a['Estimate']-crit*a['Std. Error']),hi=unname(a['Estimate']+crit*a['Std. Error']),detail=paste(unique(c(w,conv)),collapse='; '))]
 if(length(drop))r$status<-'rank_deficient' else if(length(conv))r$status<-'convergence_warning' else if(isSingular(m,tol=1e-4))r$status<-'singular_fit'
 if(!all(is.finite(unlist(r[,.(beta,se,df,p,lo,hi)]))))r$status<-'nonfinite_inference'
 if(r$status!='tested')r$p<-NA_real_
 r
}
score_data <- function(base, scores, label){
 a<-scores[is.finite(get(label)),.(PATNO,x=get(label))]
 d<-merge(base,a,by='PATNO');setorder(d,PATNO,time,EVENT_ID)
 stopifnot(!anyDuplicated(d[,.(PATNO,time)]),all(d[,uniqueN(time),by=PATNO]$V1>=2))
 d
}
if(mode=='prepare'){
 stopifnot(!file.exists(cache))
 paths<-c(frozen=file.path(P,'batch5/results/tg_v1/tg_frozen_inputs.rds'),lpi=file.path(P,'batch6/results/lpi_v1/lpi_inputs.rds'),raw=file.path(ROOT,'metabolomics_v2/cleaned/PPMI_named_adjusted_all_visits_cleaned_long.tsv.gz'),members=file.path(P,'batch9/protocol/feature_membership.tsv'))
 stopifnot(all(file.exists(paths)))
 write_t(data.table(role=names(paths),path=unname(paths),bytes=file.info(paths)$size,md5=unname(tools::md5sum(paths))),'input_manifest')
 old<-readRDS(paths['frozen']); lpi<-readRDS(paths['lpi']); members<-fread(paths['members'])
 feats<-members$feature;stopifnot(length(feats)==10,!anyDuplicated(feats),setequal(members$class,c('LPI','LPC','LPE')))
 stopifnot(all(vapply(feats,function(f)all(c('CSF','plasma')%in%old$slots[feature==f]$assay),logical(1))))
 raw<-fread(cmd=paste('gzip -dc',shQuote(paths['raw'])),na.strings=c('NA','','NaN'))
 raw<-raw[TESTNAME%in%feats];raw[,`:=`(PATNO=as.character(PATNO),EVENT_ID=CLINICAL_EVENT,assay=ifelse(TYPE=='Plasma','plasma','CSF'))]
 stopifnot(!anyDuplicated(raw[,.(PATNO,EVENT_ID,assay,TESTNAME)]))
 scorelist<-baselist<-data_list<-list(); scales<-avail<-featurevalues<-list()
 repro<-rbindlist(lapply(c('CSF','plasma'),function(a){r<-fit(lpi$ds[[a]],paste0('old_LPI16_',a));q<-lpi$reproduction[analysis==a];r[,pass:=status=='tested'&people==q$people&n==q$n&abs(beta-q$beta)<1e-7&abs(se-q$se)<1e-7];r}))
 write_t(repro,'original_reproduction');stopifnot(all(repro$pass))
 for(a in c('CSF','plasma')){
  anchor<-copy(lpi$anchors[[a]][,.(PATNO,EVENT_ID,visit_month)])
  anchor[,PATNO:=as.character(PATNO)];stopifnot(!anyDuplicated(anchor$PATNO))
  z<-matrix(NA_real_,nrow(anchor),length(feats),dimnames=list(anchor$PATNO,feats))
  for(f in feats){
   q<-raw[assay==a&TESTNAME==f,.(PATNO,EVENT_ID,TESTVALUE,log2_value,UNITS)]
   q<-merge(anchor,q,by=c('PATNO','EVENT_ID'),all.x=TRUE)
   stopifnot(nrow(q)==nrow(anchor))
   present<-is.finite(q$log2_value)
   stopifnot(all(q$UNITS[present]=='adjusted_area_ratio'),all(q$TESTVALUE[present]>0),all(abs(log2(q$TESTVALUE[present])-q$log2_value[present])<1e-9))
   slot<-old$slots[assay==a&feature==f]$slot_id;sc<-old$scales[slot_id==slot]
   stopifnot(length(slot)==1,nrow(sc)==1,sum(present)==sc$scale_n,abs(mean(q$log2_value[present])-sc$log2_mean)<1e-9,abs(sd(q$log2_value[present])-sc$log2_sd)<1e-9)
   q[,z_feature:=(log2_value-sc$log2_mean)/sc$log2_sd]
   z[,f]<-q$z_feature[match(anchor$PATNO,q$PATNO)]
   featurevalues[[paste(a,f)]]<-q[,.(assay=a,feature=f,PATNO,EVENT_ID,log2_value,z_feature)]
  }
  scores<-copy(anchor)
  for(cl in c(classes,'LYSO_without_LPI16')){
   fs<-if(cl=='LYSO_ALL')feats else if(cl=='LYSO_without_LPI16')setdiff(feats,'LPI(16:0)') else members[class==cl]$feature
   ok<-complete.cases(z[,fs,drop=FALSE])
   # Removing LPI16 uses the same complete subjects as the ten-lipid score.
   if(cl=='LYSO_without_LPI16')ok<-complete.cases(z)
   val<-rep(NA_real_,nrow(z));val[ok]<-rowMeans(z[ok,fs,drop=FALSE]);mu<-mean(val[ok]);ss<-sd(val[ok])
   stopifnot(sum(ok)>=30,is.finite(ss),ss>0)
   scores[,(cl):=(val-mu)/ss]
   scales[[paste(a,cl)]]<-data.table(assay=a,score=cl,members=paste(fs,collapse=';'),member_n=length(fs),complete_anchor_people=sum(ok),mean_of_member_z=mu,sd_of_member_z=ss)
  }
  ix<-old$slots[assay==a]$slot_id
  base<-unique(rbindlist(lapply(ix,function(i)copy(old$datasets[[i]][,setdiff(names(old$datasets[[i]]),'x'),with=FALSE]))))
  stopifnot(!anyDuplicated(base[,.(PATNO,time)]))
  baselist[[a]]<-base;scorelist[[a]]<-scores
  for(cl in classes){
   d<-score_data(base,scores,cl);data_list[[paste(a,cl,sep='__')]]<-d
   avail[[paste(a,cl)]]<-data.table(assay=a,score=cl,anchors=nrow(anchor),complete_anchors=sum(is.finite(scores[[cl]])),model_people=uniqueN(d$PATNO),model_rows=nrow(d))
  }
  write_t(scores,paste0('anchor_scores_',a));write_t(base,paste0('clinical_base_',a))
 }
 data_list[['CSF__LYSO_without_LPI16']]<-score_data(baselist$CSF,scorelist$CSF,'LYSO_without_LPI16')
 stopifnot(identical(data_list[['CSF__LYSO_ALL']][,.(PATNO,EVENT_ID,time,y)],data_list[['CSF__LYSO_without_LPI16']][,.(PATNO,EVENT_ID,time,y)]))
 sc<-scorelist$CSF;sp<-scorelist$plasma
 paired<-merge(sc[is.finite(LYSO_ALL)],sp[is.finite(LYSO_ALL)],by=c('PATNO','EVENT_ID','visit_month'),suffixes=c('_CSF','_plasma'))
 dc<-data_list[['CSF__LYSO_ALL']][PATNO%in%paired$PATNO];dp<-data_list[['plasma__LYSO_ALL']][PATNO%in%paired$PATNO]
 keys<-c('PATNO','EVENT_ID','time');common<-merge(dc[,..keys],dp[,..keys],by=keys)
 dc<-merge(dc,common,by=keys);dp<-merge(dp,common,by=keys)
 keep<-dc[,.(visits=.N),by=PATNO][visits>=2]$PATNO
 dc<-dc[PATNO%in%keep];dp<-dp[PATNO%in%keep];setorderv(dc,keys);setorderv(dp,keys)
 cv<-setdiff(names(dc),'x');stopifnot(identical(dc[,..cv],dp[,..cv]))
 paired<-paired[PATNO%in%keep];setorder(paired,PATNO)
 pair_scales<-pair_ds<-list()
 for(cl in classes)for(a in c('CSF','plasma')){
  v<-paired[[paste0(cl,'_',a)]];mu<-mean(v);ss<-sd(v);stopifnot(is.finite(ss),ss>0)
  d<-copy(if(a=='CSF')dc else dp);d[,x:=(paired[[paste0(cl,'_',a)]][match(PATNO,paired$PATNO)]-mu)/ss]
  pair_ds[[paste(a,cl,sep='__')]]<-d
  pair_scales[[paste(a,cl)]]<-data.table(assay=a,score=cl,people=nrow(paired),mean_original_score=mu,sd_original_score=ss)
 }
 fv<-rbindlist(featurevalues)
 cors<-list()
 for(a in c('CSF','plasma'))for(scope in c('complete_anchors','paired_participants')){
  w<-dcast(fv[assay==a],PATNO+EVENT_ID~feature,value.var='z_feature')
  if(scope=='paired_participants')w<-merge(w,paired[,.(PATNO,EVENT_ID)],by=c('PATNO','EVENT_ID'))
  w<-w[complete.cases(w[,..feats])];mat<-cor(as.matrix(w[,..feats]),method='spearman')
  cr<-as.data.table(as.table(mat));setnames(cr,c('feature1','feature2','rho'));cr[,`:=`(assay=a,scope=scope,people=nrow(w))];cors[[paste(a,scope)]]<-cr
 }
 write_t(fv,'anchor_feature_values');write_t(rbindlist(scales),'score_scaling');write_t(rbindlist(avail),'score_sample_flow')
 write_t(rbindlist(pair_scales),'paired_scaling');write_t(paired,'paired_anchor_scores');write_t(rbindlist(cors),'feature_correlations')
 for(a in c('CSF','plasma'))for(cl in classes)write_t(pair_ds[[paste(a,cl,sep='__')]],paste0('paired_model_data_',a,'_',cl))
 for(label in names(data_list))write_t(data_list[[label]],paste0('model_data_',label))
 saveRDS(list(full=data_list,paired=pair_ds,paired_people=nrow(paired),paired_rows=nrow(dc),classes=classes,members=members),cache)
 writeLines(unname(tools::md5sum(cache)),file.path(O,'cache_md5.txt'))
 cat('PREPARED nine score models; paired_people=',nrow(paired),' shared_rows=',nrow(dc),'\n')
}else{
 stopifnot(identical(readLines(file.path(O,'cache_md5.txt')),unname(tools::md5sum(cache))))
 inp<-readRDS(cache)
 if(mode=='pilot'){
  stopifnot(!file.exists(file.path(O,'primary_score_results.tsv')))
  primary<-rbindlist(lapply(names(inp$full),function(label){r<-fit(inp$full[[label]],label);r[,`:=`(assay=strsplit(label,'__',fixed=TRUE)[[1]][1],score=strsplit(label,'__',fixed=TRUE)[[1]][2])];r}))
  stopifnot(nrow(primary)==9)
  primary[,q_family9:=p.adjust(p,'BH',n=9)];write_t(primary,'primary_score_results')
  baseline<-rbindlist(lapply(names(inp$paired),function(label){r<-fit(inp$paired[[label]],label);r[,`:=`(assay=strsplit(label,'__',fixed=TRUE)[[1]][1],score=strsplit(label,'__',fixed=TRUE)[[1]][2])];r}))
  write_t(baseline,'paired_score_baseline')
  stopifnot(all(baseline$status=='tested'),inp$paired_people>=50)
 }else{
  stopifnot(file.exists(file.path(O,'pilot_pass.txt')),!file.exists(file.path(O,'paired_bootstrap.tsv')))
 }
 ids<-sort(unique(inp$paired[[1]]$PATNO));parts<-lapply(inp$paired,function(d)split(d,d$PATNO))
 set.seed(20261029);draws<-replicate(2000,sample(ids,length(ids),replace=TRUE),simplify=FALSE)
 boot<-function(j){
  draw<-draws[[j]]
  ans<-lapply(inp$classes,function(cl){
   rr<-lapply(c('CSF','plasma'),function(a){
    ps<-parts[[paste(a,cl,sep='__')]]
    d<-rbindlist(lapply(seq_along(draw),function(k){v<-copy(ps[[draw[k]]]);v[,PATNO:=paste0('b',k)];v}))
    stopifnot(uniqueN(d$PATNO)==length(ids))
    fit(d,paste(a,cl))
   })
   data.table(iter=j,score=cl,CSF_beta=rr[[1]]$beta,plasma_beta=rr[[2]]$beta,delta=rr[[1]]$beta-rr[[2]]$beta,qualified=rr[[1]]$status=='tested'&rr[[2]]$status=='tested',CSF_status=rr[[1]]$status,plasma_status=rr[[2]]$status,seconds=rr[[1]]$seconds+rr[[2]]$seconds)
  })
  rbindlist(ans)
 }
 if(mode=='pilot'){
  ans<-rbindlist(parallel::mclapply(1:20,boot,mc.cores=8,mc.preschedule=FALSE));write_t(ans,'paired_pilot')
  su<-ans[,.(planned=20,qualified=sum(qualified),failure_fraction=mean(!qualified),model_seconds=sum(seconds)),by=score];write_t(su,'paired_pilot_summary')
  if(all(su$failure_fraction<=.05))writeLines('pass',file.path(O,'pilot_pass.txt'))else stop('Bootstrap pilot failed; no full run permitted')
  cat('PILOT_COMPLETE seconds=',sum(ans$seconds),'\n')
 }else{
  pilot<-fread(file.path(O,'paired_pilot.tsv'))
  parts_dir<-file.path(O,'bootstrap_parts');dir.create(parts_dir,showWarnings=FALSE)
  rest<-parallel::mclapply(21:2000,function(j){r<-boot(j);fwrite(r,file.path(parts_dir,sprintf('rep_%04d.tsv',j)),sep='\t');r},mc.cores=8,mc.preschedule=FALSE)
  ans<-rbindlist(c(list(pilot),rest));setorder(ans,score,iter)
  stopifnot(nrow(ans)==8000,!anyDuplicated(ans[,.(iter,score)]));write_t(ans,'paired_bootstrap')
  baseline<-fread(file.path(O,'paired_score_baseline.tsv'))
  contrasts<-rbindlist(lapply(inp$classes,function(cl){
   d<-ans[score==cl];frac<-mean(!d$qualified)
   ci<-if(frac<=.05)quantile(d[qualified==TRUE]$delta,c(.025,.975,.00625,.99375))else rep(NA_real_,4)
   bc<-baseline[score==cl&assay=='CSF']$beta;bp<-baseline[score==cl&assay=='plasma']$beta
   data.table(score=cl,people=inp$paired_people,rows=inp$paired_rows,beta_CSF=bc,beta_plasma=bp,delta=bc-bp,lo95=unname(ci[1]),hi95=unname(ci[2]),lo_family98_75=unname(ci[3]),hi_family98_75=unname(ci[4]),bootstrap_planned=2000,bootstrap_qualified=sum(d$qualified),failure_fraction=frac,status=if(frac<=.05)'estimated_exploratory'else 'unqualified_failure_rate')
  }))
  contrasts[,family_supported:=status=='estimated_exploratory'&(hi_family98_75<0|lo_family98_75>0)]
  write_t(contrasts,'paired_class_contrasts');cat('FULL_COMPLETE\n')
 }
}
capture.output(sessionInfo(),file=file.path(O,paste0(mode,'_sessionInfo.txt')))
