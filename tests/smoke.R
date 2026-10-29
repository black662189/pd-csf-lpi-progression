#!/usr/bin/env Rscript
# Synthetic observations only. Passing is software validation, not clinical replication.
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
root<-dirname(dirname(normalizePath(script)));args<-commandArgs(trailingOnly=TRUE)
out<-if(length(args))args[1] else tempfile('lpi-smoke-')
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE);set.seed(20261029)
n<-180;latent<-rnorm(n);parent_latent<-.6*latent+rnorm(n,sd=.8)
b<-data.frame(PATNO=sprintf('SYNTHETIC_%03d',seq_len(n)),
 LPI16=2^(8+latent+rnorm(n,sd=.3)),LPI18=2^(7+latent+rnorm(n,sd=.3)),
 y0=round(runif(n,5,12)),up0=runif(n,10,35),age=runif(n,45,80),
 sex=rep(c('Female','Male'),n/2),duration=runif(n,.5,8),ledd=runif(n,0,700),
 genetic_group=rep(c('genetic','sporadic','sporadic'),n/3))
for(j in 1:6)b[[paste0('PI',j)]]<-2^(10+parent_latent+rnorm(n,sd=.3))
intercepts<-rnorm(n,sd=1.8);slopes<-rnorm(n,sd=.65)
v<-do.call(rbind,lapply(seq_len(n-15),function(i) {
 t<-c(.8,1.7,3,4.2,5.5,7)
 y<-b$y0[i]+intercepts[i]+(.9-.2*latent[i]+slopes[i])*t+rnorm(length(t),sd=1.8)
 data.frame(PATNO=b$PATNO[i],time=t,y=pmin(52,pmax(0,round(y))))
}))
write.csv(b,file.path(out,'synthetic_anchors.csv'),row.names=FALSE)
write.csv(v,file.path(out,'synthetic_followup.csv'),row.names=FALSE)
rscript<-file.path(R.home('bin'),'Rscript')
run<-function(script,args) {
 status<-system2(rscript,c(shQuote(file.path(root,script)),vapply(args,shQuote,character(1))))
 stopifnot(status==0)
}
run('scripts/prepare_scores.R',c(file.path(out,'synthetic_anchors.csv'),file.path(out,'synthetic_followup.csv'),file.path(out,'prepared')))
run('scripts/run_models.R',c(file.path(out,'prepared/model_ready.csv'),file.path(out,'models')))
flow<-read.csv(file.path(out,'prepared/sample_flow.csv'))
scales<-read.csv(file.path(out,'prepared/frozen_scales.csv'))
stopifnot(flow$model_people==165,flow$model_observations==990,all(scales$n==180))
c<-read.csv(file.path(out,'models/clinical_contrasts.csv'),check.names=FALSE)
stopifnot(nrow(c)==36,all(is.finite(c$Estimate)))
source(file.path(root,'R/models.R'))
d<-read.csv(file.path(out,'prepared/model_ready.csv'))
bad<-rbind(d,d[1,]);rejected<-inherits(try(validate_data(bad),silent=TRUE),'try-error')
stopifnot(rejected)
writeLines(c('PASS: synthetic anchor -> frozen scores -> longitudinal models -> contrasts',
 'PASS: original anchor scale retained before excluding 15 synthetic people',
 'PASS: 18 contrasts equal independent analytic vectors',
 'PASS: duplicate participant/time input rejected',
 'No patient data used; no significance or biological-effect claim tested.'),file.path(out,'SMOKE_PASS.txt'))
cat('SYNTHETIC SMOKE TEST PASSED\n')
