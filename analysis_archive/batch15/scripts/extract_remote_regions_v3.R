suppressPackageStartupMessages({library(Rsamtools);library(jsonlite)})
B <- file.path(Sys.getenv('PD_PROJECT_ROOT'), 'batch15')
plan <- fromJSON(file.path(B,'protocol/analysis_plan_v1.json'),simplifyVector=FALSE)
headers <- fromJSON(file.path(B,'protocol/source_headers.json'))
urls <- c(QTD000171='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000171/QTD000171.all.tsv.gz',
 QTD000206='https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000206/QTD000206.all.tsv.gz',
 GCST009325='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/harmonised/GCST009325.h.tsv.gz')
O <- file.path(B,'raw/remote_regions_v3');stopifnot(!dir.exists(O));dir.create(O)
P <- file.path(B,'raw/pilot_v3');dir.create(P)
read_query <- function(url,index,gr,label) {
 for(attempt in seq_len(3)) {
  # A fresh TabixFile for each request avoids retaining an HTTP seek state.
  z<-scanTabix(TabixFile(url,index=index),param=gr)[[1]]
  if(length(z)>0)return(z)
  cat('EMPTY_REMOTE_RESPONSE',label,'attempt',attempt,'\n');flush.console()
  if(attempt<3)Sys.sleep(3*attempt)
 }
 stop(paste('No records after three bounded HTTP/tabix attempts:',label))
}
pilot <- list();full <- list()
for (id in names(urls)) {
 u<-urls[[id]];index<-file.path(O,paste0(id,'.tbi'))
 previous<-file.path(dirname(B),'batch14/raw/pilot_v1',paste0(id,'.tbi'))
 if(!file.exists(previous))previous<-file.path(dirname(B),'batch14/raw/remote_regions_v1',paste0(id,'.tbi'))
 if(file.exists(previous))stopifnot(file.copy(previous,index)) else download.file(paste0(u,'.tbi'),index,method='libcurl',mode='wb',quiet=TRUE)
 tf<-TabixFile(u,index=index);h<-headerTabix(tf)
 for(reg in plan$regions) {
  chr<-if(reg$chromosome %in% h$seqnames)reg$chromosome else paste0('chr',reg$chromosome);stopifnot(chr %in% h$seqnames)
  key<-paste(reg$id,id,sep='__')
  cat('PILOT',key,as.character(Sys.time()),'\n');flush.console()
  gr<-GenomicRanges::GRanges(chr,IRanges::IRanges(reg$pilot_start,reg$pilot_end))
  x<-read_query(u,index,gr,key);stopifnot(length(x)>0)
  nc<-lengths(strsplit(x[seq_len(min(100,length(x)))],'\t'));stopifnot(all(nc==length(strsplit(headers[[id]],'\t')[[1]])))
  con<-gzfile(file.path(P,paste0(key,'.tsv.gz')),'wt');writeLines(c(headers[[id]],x),con);close(con)
  pilot[[key]]<-list(rows=length(x),interval=c(reg$pilot_start,reg$pilot_end));write_json(pilot,file.path(P,'progress.json'),pretty=TRUE,auto_unbox=TRUE)
  Sys.sleep(3)
 }
}
write_json(pilot,file.path(P,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
for(id in names(urls)) {
 tf<-TabixFile(urls[[id]],index=file.path(O,paste0(id,'.tbi')));h<-headerTabix(tf)
 for(reg in plan$regions) {
  chr<-if(reg$chromosome %in% h$seqnames)reg$chromosome else paste0('chr',reg$chromosome)
  key<-paste(reg$id,id,sep='__');cat('FULL',key,as.character(Sys.time()),'\n');flush.console()
  x<-read_query(urls[[id]],file.path(O,paste0(id,'.tbi')),GenomicRanges::GRanges(chr,IRanges::IRanges(reg$start,reg$end)),key);stopifnot(length(x)>0)
  stopifnot(all(lengths(strsplit(head(x,1000),'\t'))==length(strsplit(headers[[id]],'\t')[[1]])))
  f<-file.path(O,paste0(key,'.tsv.gz'));con<-gzfile(f,'wt');writeLines(c(headers[[id]],x),con);close(con)
  full[[key]]<-list(rows=length(x),region_id=reg$id,source_id=id,url=urls[[id]],interval=c(reg$start,reg$end),bytes=file.info(f)$size)
  write_json(full,file.path(O,'progress.json'),pretty=TRUE,auto_unbox=TRUE)
  cat('DONE',key,length(x),'\n');flush.console();rm(x);gc();Sys.sleep(3)
 }
}
write_json(full,file.path(O,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
