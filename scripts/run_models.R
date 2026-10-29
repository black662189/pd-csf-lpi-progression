#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop('Usage: Rscript scripts/run_models.R MODEL_READY.csv NEW_OUTPUT_DIRECTORY')
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
root <- dirname(dirname(normalizePath(script)))
source(file.path(root,'R/models.R'))
input <- normalizePath(args[1]); out <- args[2]
stopifnot(!dir.exists(out)); dir.create(out,recursive=TRUE)
d <- validate_data(read.csv(input,colClasses=c(PATNO='character'),check.names=FALSE))
writeLines(c('This run is a refit from the supplied model-ready data.',
 'No discovery-family q values are recomputed by this focused runner.'),file.path(out,'SCOPE.txt'))
coefficients <- metrics <- contrasts <- list()
for(context in c('score','PI_conditional')) {
  dd <- if(context=='score') d else droplevels(d[is.finite(d$parent),,drop=FALSE])
  stopifnot(nrow(dd)>0)
  for(structure in c('diagonal','correlated','quadratic')) {
    label <- paste(context,structure,sep='__')
    m <- fit_model(dd,context=='PI_conditional',structure)
    a <- as.data.frame(coef(summary(m))); a$term <- rownames(a); a$model <- label
    a$lo <- a$Estimate-qt(.975,a$df)*a[['Std. Error']]
    a$hi <- a$Estimate+qt(.975,a$df)*a[['Std. Error']]
    coefficients[[label]] <- a
    metrics[[label]] <- model_metrics(m,dd,label)
    contrasts[[label]] <- clinical_contrasts(m,dd,label)
    writeLines(deparse(formula(m)),file.path(out,paste0(label,'_formula.txt')))
  }
}
dd <- droplevels(d[is.finite(d$parent),,drop=FALSE]); fits <- list()
for(label in c('clinical','clinical_LPI','clinical_PI','clinical_PI_LPI')) {
  m <- fit_model(dd,parent=grepl('PI',label,fixed=TRUE) && label!='clinical_LPI',
                 exposure=grepl('LPI',label,fixed=TRUE),reml=FALSE)
  fits[[label]] <- m; metrics[[label]] <- model_metrics(m,dd,label)
}
pairs <- list(LPI_block=c('clinical','clinical_LPI'),LPI_block_given_PI=c('clinical_PI','clinical_PI_LPI'))
comparisons <- lapply(names(pairs),function(label) {
  ids <- pairs[[label]]; a <- as.data.frame(anova(fits[[ids[1]]],fits[[ids[2]]],refit=FALSE))
  a$model <- ids; a$comparison <- label; rownames(a)<-NULL; a
})
write.csv(do.call(rbind,coefficients),file.path(out,'fixed_effects.csv'),row.names=FALSE)
write.csv(do.call(rbind,metrics),file.path(out,'model_metrics.csv'),row.names=FALSE)
write.csv(do.call(rbind,contrasts),file.path(out,'clinical_contrasts.csv'),row.names=FALSE)
write.csv(do.call(rbind,comparisons),file.path(out,'likelihood_comparisons.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines(c(paste('input_md5',tools::md5sum(input)), 'models 10','clinical_contrasts 36',
             'contrast_analytic_checks 18','status completed'),file.path(out,'COMPLETE.txt'))
cat('Completed 10 qualified model fits and 36 clinical contrasts.\n')
