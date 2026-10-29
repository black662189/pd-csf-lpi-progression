.libPaths(c(Sys.getenv('PD_COLOC_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(coloc);library(data.table);library(jsonlite);library(ggplot2)})
setDTthreads(1)
B <- file.path(Sys.getenv('PD_PROJECT_ROOT'), 'batch14')
H <- file.path(B,'results/harmonized_v1');O <- file.path(B,'results/coloc_v1')
stopifnot(file.exists(file.path(H,'COMPLETE.json')),!dir.exists(O));dir.create(O)
dir.create(file.path(O,'pair_details'));dir.create(file.path(O,'figures'))
source(file.path(B,'scripts/coloc_parameter_helpers.R'))
plan <- fromJSON(file.path(B,'protocol/analysis_plan_v1.json'))
meta <- fromJSON(file.path(H,'trait_qc.json'),simplifyVector=FALSE)
names(meta) <- vapply(meta,`[[`,'','trait')
traits <- names(meta);D <- setNames(lapply(traits,function(t)fread(file.path(H,paste0(t,'.tsv.gz')))),traits)
for (t in traits) {
  d <- D[[t]];stopifnot(!anyDuplicated(d$snp),all(is.finite(d$beta)),all(d$se>0))
  if(meta[[t]]$kind!='cc') {
    sdcheck <- estimate_sdy_variable_n(d$se,d$maf,d$n)$sdY
    stopifnot(isTRUE(all.equal(sdcheck,meta[[t]]$sdY_estimated,tolerance=1e-10)))
  }
}
pairs <- list(c('cortex','LPI18'),c('cortex','LPI16'),c('nigra','LPI18'),c('nigra','LPI16'),c('cortex','PD'),c('nigra','PD'),c('LPI18','PD'),c('LPI16','PD'))
warning_records <- list();wi <- 0L
run_pair <- function(a,b,label,p12=1e-5,scale=1,expression_estimated=FALSE,snps=NULL,window=NULL,p1=1e-4,p2=1e-4,save_details=FALSE) {
  common <- intersect(D[[a]]$snp,D[[b]]$snp)
  if(!is.null(snps))common <- intersect(common,snps)
  if(!is.null(window)) common <- intersect(common,D[[a]][position>=window[1]&position<=window[2],snp])
  stopifnot(length(common)>1)
  x <- D[[a]][match(common,snp)];y <- D[[b]][match(common,snp)]
  stopifnot(identical(x$snp,y$snp),identical(x$effect_allele,y$effect_allele),identical(x$other_allele,y$other_allele),all(x$position==y$position))
  make <- function(name,d) {
    m <- meta[[name]]
    if(m$kind=='cc')return(make_pd_dataset(d))
    sd_y <- if(m$kind=='lipid') m$sdY_main*scale else if(expression_estimated) m$sdY_estimated else 1
    make_quant_dataset(d,sd_y)
  }
  ds1 <- make(a,x);ds2 <- make(b,y);messages <- character()
  invisible(capture.output(res <- withCallingHandlers(coloc.abf(ds1,ds2,p1=p1,p2=p2,p12=p12),warning=function(w){messages <<- c(messages,conditionMessage(w));invokeRestart('muffleWarning')})))
  if(length(messages)) {
    wi <<- wi+1L;warning_records[[wi]] <<- data.table(pair=paste(a,b,sep='__'),analysis=label,p12_requested=p12,scale=scale,message=unique(messages))
  }
  pp <- as.numeric(res$summary[paste0('PP.H',0:4,'.abf')]);stopifnot(length(pp)==5,all(is.finite(pp)),abs(sum(pp)-1)<1e-8)
  pr <- res$priors
  if(is.null(pr))stop('Installed coloc did not return actual priors; inspect version before reporting')
  row <- data.table(trait1=a,trait2=b,analysis=label,n_snps=length(common),p1_requested=p1,p2_requested=p2,p12_requested=p12,
       p1_used=unname(pr['p1']),p2_used=unname(pr['p2']),p12_used=unname(pr['p12']),lipid_sdY_factor=scale,
       expression_scale=if(expression_estimated)'summary_estimate'else'INT_unit',
       sdY1=if(is.null(ds1$sdY))NA_real_ else ds1$sdY,sdY2=if(is.null(ds2$sdY))NA_real_ else ds2$sdY,
       H0=pp[1],H1=pp[2],H2=pp[3],H3=pp[4],H4=pp[5],
       min_p1=min(x$p_from_beta_se),min_p2=min(y$p_from_beta_se),
       n_freq_diff_gt_0_2=sum(abs(x$eaf-y$eaf)>.2),
       interpretation=if(pp[5]>=.8)'H4_support_under_single_signal_model'else if(pp[4]>=.8)'H3_support_distinct_signals_under_model'else if(sum(pp[1:3])>=.8)'insufficient_evidence_for_both_traits'else'unresolved_H3_H4_or_mixed')
  if(save_details) {
    det <- as.data.table(res$results);det <- merge(det,x[,.(snp,position,effect_allele,other_allele,beta1=beta,se1=se,p1=p_original,eaf1=eaf)],by='snp')
    det <- merge(det,y[,.(snp,beta2=beta,se2=se,p2=p_original,eaf2=eaf)],by='snp')
    fwrite(det,file.path(O,'pair_details',paste0(a,'__',b,'.tsv.gz')),sep='\t')
    saveRDS(res,file.path(O,'pair_details',paste0(a,'__',b,'.rds')))
  }
  row
}
# Implementation pilot: real subset, not used as a scientific result.
pilot_snps <- head(intersect(D$cortex$snp,D$LPI18$snp),2000)
smoke <- run_pair('cortex','LPI18','implementation_pilot',snps=pilot_snps)
fwrite(smoke,file.path(O,'implementation_pilot.tsv'),sep='\t')
main <- rbindlist(lapply(pairs,function(p)run_pair(p[1],p[2],'primary',save_details=TRUE)))
stopifnot(nrow(main)==8);fwrite(main,file.path(O,'primary_colocalization.tsv'),sep='\t')
sens <- list();si <- 0L
for(p in pairs) {
  has_lipid <- any(vapply(meta[p],function(m)m$kind=='lipid',TRUE))
  has_expr <- any(vapply(meta[p],function(m)m$kind=='expression',TRUE))
  for(pr in plan$priors$p12_sensitivity)for(sc in if(has_lipid)c(.5,1,2)else 1)for(ex in if(has_expr)c(FALSE,TRUE)else FALSE) {
    si <- si+1L;sens[[si]] <- run_pair(p[1],p[2],'prior_and_scale',p12=pr,scale=sc,expression_estimated=ex)
  }
  si <- si+1L;sens[[si]] <- run_pair(p[1],p[2],'TSS_plusminus_500kb',window=plan$sensitivity_window)
  si <- si+1L;sens[[si]] <- run_pair(p[1],p[2],'fixed_lower_p1_p2',p1=5e-5,p2=5e-5)
}
sens <- rbindlist(sens);fwrite(sens,file.path(O,'sensitivity_colocalization.tsv'),sep='\t')
# Identical three-trait variant sets; these are pairwise compatibility checks, not joint posteriors.
triple <- list();ti <- 0L
for(e in c('cortex','nigra'))for(l in c('LPI18','LPI16')) {
  ss <- Reduce(intersect,lapply(D[c(e,l,'PD')],`[[`,'snp'))
  for(p in list(c(e,l),c(e,'PD'),c(l,'PD'))) {
    ti <- ti+1L;z <- run_pair(p[1],p[2],'three_trait_common_variants',snps=ss);z[,triple_id:=paste(e,l,'PD',sep='__')];triple[[ti]] <- z
  }
}
fwrite(rbindlist(triple),file.path(O,'three_trait_compatibility_NOT_joint_posterior.tsv'),sep='\t')
if(length(warning_records))fwrite(rbindlist(warning_records),file.path(O,'warnings.tsv'),sep='\t')

plotdat <- rbindlist(lapply(traits,function(t){x<-copy(D[[t]]);x[,trait:=t];x[,.(trait,position,p_original)]}))
plotdat[,trait:=factor(trait,levels=c('cortex','nigra','LPI18','LPI16','PD'))]
plotdat[,minus_log10_p:=-log10(pmax(p_original,1e-300))]
g <- ggplot(plotdat,aes(position/1e6,minus_log10_p))+geom_point(size=.55,alpha=.6,color='#2B6480')+
  geom_vline(xintercept=plan$gene_end/1e6,color='#B14B40',linetype=2)+facet_grid(trait~.,scales='free_y')+
  labs(title='MBOAT7 region: expression, target LPI and PD risk',subtitle='Dense regional summary statistics; GRCh38; dashed line = MBOAT7 TSS',x='Chromosome 19 position (Mb)',y='-log10(original association P)',
       caption='Different y-axis scales. Brain LPI is postmortem tissue; PD GWAS concerns risk.\nGTEx tissues are not independent cohorts. These plots do not establish mediation.')+
  theme_bw(base_size=11)+theme(strip.background=element_rect(fill='#e7edf0'),plot.caption=element_text(hjust=0),panel.grid.minor=element_blank())
ggsave(file.path(O,'figures/regional_association.pdf'),g,width=10,height=10)
ggsave(file.path(O,'figures/regional_association.png'),g,width=10,height=10,dpi=170)
post <- melt(main[,.(trait1,trait2,H0,H1,H2,H3,H4)],id.vars=c('trait1','trait2'),variable.name='hypothesis',value.name='posterior')
post[,pair:=factor(paste(trait1,trait2,sep=' / '),levels=rev(paste(main$trait1,main$trait2,sep=' / ')))]
g2 <- ggplot(post,aes(pair,posterior,fill=hypothesis))+geom_col(width=.7)+coord_flip()+
  scale_fill_manual(values=c(H0='#d8e1e8',H1='#49799b',H2='#83b5bb',H3='#d69f62',H4='#ab4a49'))+
  labs(title='Pairwise colocalization at the MBOAT7 region',x=NULL,y='Posterior probability',caption='H0-H2: at least one association not supported. H3: distinct variants. H4: shared variant.\nSingle-causal-variant model; actual priors and sensitivity results are reported separately.')+
  theme_bw(base_size=11)+theme(plot.caption=element_text(hjust=0))
ggsave(file.path(O,'figures/colocalization_posteriors.pdf'),g2,width=10,height=6.5)
ggsave(file.path(O,'figures/colocalization_posteriors.png'),g2,width=10,height=6.5,dpi=170)
sink(file.path(O,'sessionInfo.txt'));print(sessionInfo());print(.libPaths());sink()
capture.output(getFromNamespace('adjust_prior','coloc'),file=file.path(O,'installed_adjust_prior.txt'))
capture.output(getFromNamespace('approx.bf.estimates','coloc'),file=file.path(O,'installed_ABF_formula.txt'))
write_json(list(complete=TRUE,primary_pairs=nrow(main),sensitivity_rows=nrow(sens),three_trait_pair_checks=12,coloc_version=as.character(packageVersion('coloc')),
                any_primary_H4_ge_0_8=any(main$H4>=.8),LD_based_multisignal_analysis_performed=FALSE,mediation_claim=FALSE),file.path(O,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE)
print(main)
