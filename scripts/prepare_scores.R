#!/usr/bin/env Rscript
# Input rows must already be clinically eligible original molecular anchors.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop('Usage: Rscript scripts/prepare_scores.R ANCHORS.csv FOLLOWUP.csv NEW_OUTPUT_DIRECTORY')
b <- read.csv(args[1],colClasses=c(PATNO='character'),check.names=FALSE)
v <- read.csv(args[2],colClasses=c(PATNO='character'),check.names=FALSE)
out <- args[3]; stopifnot(!dir.exists(out),!anyDuplicated(b$PATNO),!anyDuplicated(v[c('PATNO','time')]))
required <- c('PATNO','LPI16','LPI18',paste0('PI',1:6),'y0','up0','age','sex','duration','ledd','genetic_group')
stopifnot(all(required %in% names(b)),all(c('PATNO','time','y') %in% names(v)),
          all(b$age>0),all(b$duration>=0),all(b$ledd>=0))
dir.create(out,recursive=TRUE); scales <- list()
score <- function(fields,label) {
  zz <- sapply(fields,function(f) {
    x <- b[[f]]; x[!is.finite(x)|x<=0] <- NA_real_; x<-log2(x)
    mu<-mean(x,na.rm=TRUE); s<-sd(x,na.rm=TRUE); stopifnot(is.finite(mu),is.finite(s),s>0)
    scales[[length(scales)+1]] <<- data.frame(variable=f,mean=mu,SD=s,n=sum(is.finite(x)),scale='log2 species')
    (x-mu)/s
  })
  x<-rowMeans(zz); mu<-mean(x,na.rm=TRUE); s<-sd(x,na.rm=TRUE)
  stopifnot(is.finite(s),s>0)
  scales[[length(scales)+1]] <<- data.frame(variable=label,mean=mu,SD=s,n=sum(is.finite(x)),scale='mean species z')
  (x-mu)/s
}
b$x <- score(c('LPI16','LPI18'),'LPI_score')
b$parent <- score(paste0('PI',1:6),'PI_score')
b$age0 <- b$age/10; b$duration0<-b$duration/5; b$ledd0<-log1p(b$ledd)
columns <- c('PATNO','x','parent','y0','up0','age0','sex','duration0','ledd0','genetic_group')
eligible <- complete.cases(b[setdiff(columns,'parent')]) & is.finite(b$x)
v <- v[is.finite(v$time)&v$time>0&is.finite(v$y)&v$y>=0&v$y<=52,,drop=FALSE]
d <- merge(v,b[eligible,columns],by='PATNO'); counts<-table(d$PATNO)
d<-d[d$PATNO %in% names(counts[counts>=2]),,drop=FALSE]
d<-d[order(d$PATNO,d$time),,drop=FALSE]
write.csv(d,file.path(out,'model_ready.csv'),row.names=FALSE)
write.csv(do.call(rbind,scales),file.path(out,'frozen_scales.csv'),row.names=FALSE)
write.csv(data.frame(eligible_anchor_rows=nrow(b),model_people=length(unique(d$PATNO)),
                     model_observations=nrow(d)),file.path(out,'sample_flow.csv'),row.names=FALSE)
cat('Scores frozen before outcome/follow-up exclusions.\n')
