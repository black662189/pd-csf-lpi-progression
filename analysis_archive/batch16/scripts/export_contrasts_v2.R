#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
.libPaths(c(Sys.getenv('PD_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(data.table);library(lmerTest);library(jsonlite)})
setDTthreads(1)
P<-Sys.getenv('PD_PROJECT_ROOT');B<-file.path(P,'batch16');I<-file.path(B,'results/diagnostics_v1');O<-file.path(B,'results/contrasts_v2')
stopifnot(!dir.exists(O));dir.create(O)
design<-function(m,d)model.matrix(delete.response(terms(reformulas::nobars(formula(m)))),data=d,contrasts.arg=attr(getME(m,'X'),'contrasts'))[,names(fixef(m)),drop=FALSE]
rows<-list();vectors<-list();checks<-list()
for(context in c('score','PI_conditional')){
 path<-if(context=='score')'batch9/results/classes_v1/model_data_CSF__LPI.tsv'else'batch11/results/parent_v1/full_model_data_CSF__LPI.tsv'
 d<-fread(file.path(P,path));d[,`:=`(PATNO=factor(PATNO),sex=factor(sex),genetic_group=factor(genetic_group))];setorder(d,PATNO,time);pp<-d[!duplicated(PATNO)]
 for(structure in c('diagonal','correlated','quadratic')){
  label<-paste(context,structure,sep='__');m<-readRDS(file.path(I,paste0(label,'.rds')))
  grid<-function(t,x_value=NULL){g<-copy(pp);g[,`:=`(time=t,time_sq=(t/5)^2)];if(!is.null(x_value))set(g,j='x',value=rep(x_value,nrow(g)));g}
  for(h in c(1,3,5)){
   lx<-(colMeans(design(m,grid(h,.5)))-colMeans(design(m,grid(h,-.5))))-(colMeans(design(m,grid(0,.5)))-colMeans(design(m,grid(0,-.5))))
   expected<-rep(0,length(fixef(m)));names(expected)<-names(fixef(m));expected['time:x']<-h
   if(structure=='quadratic')expected['time_sq:x']<-(h/5)^2
   err<-max(abs(lx-expected));stopifnot(err<1e-12)
   lc<-colMeans(design(m,grid(h)))-colMeans(design(m,grid(0)))
   for(kind in c('change_difference_per_1SD','cohort_standardized_change')){
    L<-if(kind=='change_difference_per_1SD')lx else lc
    a<-as.data.table(contest1D(m,L,confint=TRUE));stopifnot(all(is.finite(as.numeric(a[1]))),abs(a$Estimate-sum(L*fixef(m)))<1e-10)
    a[,`:=`(model=label,contrast=kind,year=h)];key<-paste(label,h,kind);rows[[key]]<-a
    vectors[[key]]<-data.table(model=label,contrast=kind,year=h,term=names(L),weight=as.numeric(L))
   }
   checks[[paste(label,h)]]<-data.table(model=label,year=h,contrast_vector_max_abs_error=err,pass=TRUE)
  }
 }
}
fwrite(rbindlist(rows),file.path(O,'clinical_contrasts.tsv'),sep='\t');fwrite(rbindlist(vectors),file.path(O,'contrast_vectors.tsv'),sep='\t');fwrite(rbindlist(checks),file.path(O,'contrast_checks.tsv'),sep='\t')
write_json(list(status='verified',rows=length(rows),models_refit=0,supersedes='diagnostics_v1/clinical_contrasts.tsv only',reason='Original exposure-setting helper shadowed its x parameter with the data.table x column. Corrected by explicit set() using x_value; all contrast vectors checked against analytic weights.'),file.path(O,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
print(rbindlist(rows)[year==5]);capture.output(sessionInfo(),file=file.path(O,'sessionInfo.txt'))
