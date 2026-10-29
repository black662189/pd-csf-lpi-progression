#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
args<-commandArgs(TRUE);node<-args[1];stopifnot(node%in%c('200','230'))
P<-if(node=='200')Sys.getenv('PD_PROJECT_ROOT') else Sys.getenv('PD_PROJECT_ROOT')
if(node=='200').libPaths(c(Sys.getenv('PD_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(limma);library(data.table)})
setDTthreads(1);set.seed(20261029)
B<-file.path(P,'batch10');O<-file.path(B,'results',paste0('rna_',node,'_v1'));stopifnot(!dir.exists(O));dir.create(O,recursive=TRUE)
wt<-function(d,n)fwrite(d,file.path(O,paste0(n,'.tsv')),sep='\t')
src<-if(node=='200')file.path(P,'development/rna_v1/sn__sublineage__Dopaminergic__n20__full/voom_fit.rds') else file.path(P,'validation/results_v1/primary_DA/fit_state.rds')
z<-readRDS(src);v<-z$v;design<-z$design;md<-z$metadata
stopifnot(ncol(v)==nrow(design),qr(design)$rank==ncol(design),'groupPD'%in%colnames(design),ncol(v)==if(node=='200')16L else 13L)
stopifnot(!anyDuplicated(rownames(v)),identical(colnames(v),rownames(design)))
donor<-if(node=='200')md$donor_id else md$donor_id
stopifnot(length(donor)==ncol(v),!anyDuplicated(donor))
oldpath<-file.path(dirname(src),'all_gene_effects.tsv.gz');old<-fread(cmd=paste('gzip -dc',shQuote(oldpath)))
gene_col<-if(node=='200')'gene' else 'symbol'
stopifnot(all(rownames(v)%in%old[[gene_col]]))
# Confirm cache/tested-gene alignment and exact group coefficient before gene-set testing.
lc<-z$fit$coefficients[,'groupPD'];expected<-old$logFC[match(names(lc),old[[gene_col]])]
stopifnot(length(lc)==nrow(v),max(abs(lc-expected))<1e-7)
wt(data.table(input=src,md5=unname(tools::md5sum(src)),genes=nrow(v),donors=ncol(v),PD=sum(md$group=='PD'),controls=sum(md$group=='Control'),design=paste(colnames(design),collapse=';'),coefficient_max_difference=max(abs(lc-expected))),'cache_verification')
gmt<-strsplit(readLines(file.path(B,'protocol/selected_pathways.gmt')),'\t',fixed=TRUE)
sets<-setNames(lapply(gmt,function(x)x[-c(1,2)]),vapply(gmt,`[`,character(1),1))
index<-ids2indices(sets,rownames(v));stopifnot(length(index)==3,all(lengths(index)>=10))
wt(rbindlist(lapply(names(sets),function(n)data.table(set=n,gene=sets[[n]],expressed=sets[[n]]%in%rownames(v)))),'set_membership')
# Representative software check with a synthetic two-group design; no outcome-dependent tuning.
set.seed(20261029);e<-matrix(rnorm(500*12),500,12);ds<-model.matrix(~factor(rep(c('C','P'),each=6)))
pilot<-mroast(e,index=list(pilot=1:40),design=ds,contrast=2,nrot=99,set.statistic='mean')
stopifnot(nrow(pilot)==1,is.finite(pilot$PValue),pilot$PValue>=0,pilot$PValue<=1)
wt(as.data.table(pilot,keep.rownames='set'),'synthetic_pilot')
set.seed(20261029)
rr<-mroast(v,index=index,design=design,contrast=match('groupPD',colnames(design)),nrot=9999,set.statistic='mean')
r<-as.data.table(rr,keep.rownames='set');r[,`:=`(study=if(node=='200')'GSE253462' else 'GSE243639',donors=ncol(v),PD=sum(md$group=='PD'),controls=sum(md$group=='Control'),nrot=9999,seed=20261029)]
r[,mean_logFC:=vapply(set,function(n)mean(lc[index[[n]]]),numeric(1))]
r[,median_logFC:=vapply(set,function(n)median(lc[index[[n]]]),numeric(1))]
wt(r,'roast_results_raw');capture.output(sessionInfo(),file=file.path(O,'sessionInfo.txt'))
writeLines('3 fixed sets tested. Final inference uses combined BH6, not within-cohort FDR.',file.path(O,'COMPLETE.txt'))
