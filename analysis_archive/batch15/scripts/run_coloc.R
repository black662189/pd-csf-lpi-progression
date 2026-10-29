.libPaths(c(Sys.getenv('PD_COLOC_R_LIB'),.libPaths()))
suppressPackageStartupMessages({library(coloc);library(data.table);library(jsonlite);library(parallel)})
setDTthreads(1)
B <- file.path(Sys.getenv('PD_PROJECT_ROOT'), 'batch15')
H<-file.path(B,'results/harmonized_v1');O<-file.path(B,'results/coloc_v1')
stopifnot(file.exists(file.path(H,'COMPLETE.json')),!dir.exists(O));dir.create(O);dir.create(file.path(O,'pair_details'))
source(file.path(B,'scripts/coloc_parameter_helpers.R'))
plan<-fromJSON(file.path(B,'protocol/analysis_plan_v1.json'),simplifyVector=FALSE)
ml<-fromJSON(file.path(H,'dataset_manifest.json'),simplifyVector=FALSE);names(ml)<-vapply(ml,`[[`,'','dataset')
regions<-plan$regions;names(regions)<-vapply(regions,`[[`,'','id')
read_d<-function(m) {d<-fread(m$path);stopifnot(!anyDuplicated(d$snp));d}
lse<-function(v) {u<-max(v);u+log(sum(exp(v-u)))}
make_ds<-function(m,d,lf=1,ee=FALSE,le=FALSE) {
 if(m$kind=='cc')return(make_pd_dataset(d))
 sd<-if(m$kind=='lipid'){if(le)m$sdY_estimated else m$sdY_main*lf}else if(ee)m$sdY_estimated else 1
 make_quant_dataset(d,sd)
}
fit<-function(a,b,x,y,label,p12=1e-5,lf=1,ee=FALSE,le=FALSE,p1=5e-5,p2=5e-5,ss=NULL,win=NULL,remove_freq=FALSE,details=FALSE) {
 common<-intersect(x$snp,y$snp)
 if(!is.null(ss))common<-intersect(common,ss)
 if(!is.null(win))common<-intersect(common,x[position>=win[1]&position<=win[2],snp])
 u<-x[match(common,snp)];v<-y[match(common,snp)]
 if(remove_freq) {keep<-is.finite(u$eaf)&is.finite(v$eaf)&abs(u$eaf-v$eaf)<=.2;u<-u[keep];v<-v[keep];common<-common[keep]}
 base<-data.table(pair_id=paste(a$dataset,b$dataset,sep='___'),region_id=a$region_id,dataset1=a$dataset,dataset2=b$dataset,gene_id=a$gene_id,phenotype=a$phenotype,
  source1=a$source_id,source2=b$source_id,analysis=label,n_snps=length(common),p1_requested=p1,p2_requested=p2,p12_requested=p12,lipid_sdY_factor=lf,expression_scale=if(ee)'summary_estimate'else'INT_unit',lipid_scale=if(le)'new_region_estimate'else'frozen_batch14')
 if(length(common)<2) {base[,`:=`(status='not_testable',reason='fewer_than_two_common_variants',H0=NA_real_,H1=NA_real_,H2=NA_real_,H3=NA_real_,H4=NA_real_)];return(list(row=base,warnings=NULL))}
 stopifnot(identical(u$snp,v$snp),identical(u$effect_allele,v$effect_allele),identical(u$other_allele,v$other_allele),all(u$position==v$position))
 ds1<-make_ds(a,u,lf,ee,le);ds2<-make_ds(b,v,lf,ee,le);ww<-character()
 invisible(capture.output(res<-withCallingHandlers(coloc.abf(ds1,ds2,p1=p1,p2=p2,p12=p12),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')})))
 pp<-as.numeric(res$summary[paste0('PP.H',0:4,'.abf')]);stopifnot(length(pp)==5,all(is.finite(pp)),abs(sum(pp)-1)<1e-10)
 pr<-res$priors;stopifnot(!is.null(pr))
 base[,`:=`(status='computed',reason='',p1_used=unname(pr['p1']),p2_used=unname(pr['p2']),p12_used=unname(pr['p12']),
  sdY1=if(is.null(ds1$sdY))NA_real_ else ds1$sdY,sdY2=if(is.null(ds2$sdY))NA_real_ else ds2$sdY,
  H0=pp[1],H1=pp[2],H2=pp[3],H3=pp[4],H4=pp[5],conditional_H4=pp[5]/sum(pp[4:5]),
  position_min=min(u$position),position_max=max(u$position),fraction_dataset1=nrow(u)/nrow(x),fraction_dataset2=nrow(v)/nrow(y),
  min_p1_original=min(u$p_original),min_p2_original=min(v$p_original),min_p1_beta_se=min(u$p_from_beta_se),min_p2_beta_se=min(v$p_from_beta_se),
  n_frequency_diff_gt_0_2=sum(abs(u$eaf-v$eaf)>.2,na.rm=TRUE))]
 if(a$kind=='lipid'||b$kind=='lipid') {
  lipid<-if(a$kind=='lipid')x else y;m<-if(a$kind=='lipid')a else b;reg<-regions[[a$region_id]]
  top<-head(lipid[order(p_original,snp),snp],50)
  z<-lipid$beta/lipid$se;V<-lipid$se^2;W<-(.15*m$sdY_main)^2;r<-W/(W+V);abf<-.5*(log1p(-r)+r*z^2)
  retained<-lipid$snp %in% common
  base[,`:=`(lipid_lead_retained=reg$lead_key %in% common,lipid_top50_retained=sum(top %in% common),lipid_top50_total=length(top),
   lipid_ABF_weight_retained=exp(lse(abf[retained])-lse(abf)),lipid_full_min_p=min(lipid$p_original))]
 }
 if(a$kind=='expression')base[,expression_lead_retained:=x$snp[which.min(x$p_original)] %in% common]
 if(details) {
  det<-as.data.table(res$results);fwrite(det,file.path(O,'pair_details',paste0(a$dataset,'___',b$dataset,'.tsv.gz')),sep='\t')
 }
 list(row=base,warnings=if(length(ww))data.table(pair_id=base$pair_id,analysis=label,p12=p12,lipid_factor=lf,message=unique(ww))else NULL)
}
jobs<-list()
for(m in ml)if(m$kind=='expression')jobs[[length(jobs)+1L]]<-list(a=m,b=ml[[paste0(m$region_id,'__lipid')]])
for(reg in regions)jobs[[length(jobs)+1L]]<-list(a=ml[[paste0(reg$id,'__lipid')]],b=ml[[paste0(reg$id,'__PD')]])
run_job<-function(job,save_details=TRUE) {
 a<-job$a;b<-job$b;x<-read_d(a);y<-read_d(b)
 for(q in list(list(m=a,d=x),list(m=b,d=y)))if(q$m$kind!='cc'&&nrow(q$d)>0)stopifnot(isTRUE(all.equal(estimate_sdy_variable_n(q$d$se,q$d$maf,q$d$n)$sdY,q$m$sdY_estimated,tolerance=1e-10)))
 main<-fit(a,b,x,y,'primary',details=save_details);sens<-list();warnings<-list(main$warnings)
 has_expr<-a$kind=='expression'||b$kind=='expression';k<-0L
 for(p in unlist(plan$sensitivity$p12))for(lf in c(.5,1,2))for(ee in if(has_expr)c(FALSE,TRUE)else FALSE) {
  k<-k+1L;res<-fit(a,b,x,y,'prior_and_scale',p12=p,lf=lf,ee=ee);sens[[k]]<-res$row;warnings[[length(warnings)+1L]]<-res$warnings
 }
 reg<-regions[[a$region_id]]
 extra<-list(fit(a,b,x,y,'lipid_full_region_estimated_scale',le=TRUE),fit(a,b,x,y,'lead_plusminus_500kb',win=c(reg$lead_position-500000,reg$lead_position+500000)),
  fit(a,b,x,y,'default_p1_p2_1e-4',p1=1e-4,p2=1e-4),fit(a,b,x,y,'exclude_EAF_difference_gt_0.2',remove_freq=TRUE))
 for(res in extra) {k<-k+1L;sens[[k]]<-res$row;warnings[[length(warnings)+1L]]<-res$warnings}
 list(main=main$row,sensitivity=rbindlist(sens,fill=TRUE),warnings=rbindlist(warnings,fill=TRUE))
}
# Real-input implementation pilot with a full first phenotype; saved separately from discovery outputs.
stopifnot(length(jobs)>2)
pilot<-run_job(jobs[[1]],save_details=FALSE)
fwrite(pilot$main,file.path(O,'implementation_pilot.tsv'),sep='\t')
cat('PILOT_OK; jobs=',length(jobs),'\n');flush.console()
ans<-mclapply(jobs,run_job,mc.cores=4,mc.preschedule=TRUE)
bad<-vapply(ans,inherits,TRUE,'try-error');if(any(bad))stop(paste(unlist(ans[bad]),collapse='\n'))
main<-rbindlist(lapply(ans,`[[`,'main'),fill=TRUE);sens<-rbindlist(lapply(ans,`[[`,'sensitivity'),fill=TRUE);warnings<-rbindlist(lapply(ans,`[[`,'warnings'),fill=TRUE)
stopifnot(nrow(main)==length(jobs));fwrite(main,file.path(O,'primary_all.tsv'),sep='\t');fwrite(sens,file.path(O,'sensitivity_all.tsv'),sep='\t')
if(nrow(warnings))fwrite(warnings,file.path(O,'warnings.tsv'),sep='\t')
ss<-sens[status=='computed',.(sensitivity_H4_min=min(H4),sensitivity_H4_max=max(H4),non_extreme_prior_H4_max=max(H4[p12_requested<5e-5]),sensitivity_n=.N),by=pair_id]
ranked<-merge(main,ss,by='pair_id',all.x=TRUE)
ranked[,classification:=fcase(status!='computed','not_testable',H4>=.8,'stronger_model_support',H4>=.5,'exploratory_model_support',sensitivity_H4_max>=.5,'parameter_sensitive_only',default='not_supported_in_tested_data')]
ranked[,coverage_label:=fifelse(lipid_lead_retained,'fixed_lipid_lead_retained','fixed_lipid_lead_missing_or_outside_cis_coverage')]
setorder(ranked,-H4);fwrite(ranked,file.path(O,'ranked_all.tsv'),sep='\t')
candidates<-ranked[gene_id!='' & (H4>=.5|sensitivity_H4_max>=.5)]
fwrite(candidates,file.path(O,'candidate_review.tsv'),sep='\t')
# Disease evidence is assessed after nomination; it is not a nomination filter.
follow<-list();triple<-list()
if(nrow(candidates))for(i in seq_len(nrow(candidates))) {
 a<-ml[[candidates$dataset1[i]]];l<-ml[[candidates$dataset2[i]]];p<-ml[[paste0(a$region_id,'__PD')]]
 x<-read_d(a);y<-read_d(l);z<-read_d(p)
 rr<-fit(a,p,x,z,'candidate_expression_PD',details=TRUE);follow[[length(follow)+1L]]<-rr$row
 common<-Reduce(intersect,list(x$snp,y$snp,z$snp))
 for(j in list(list(a=a,b=l,x=x,y=y),list(a=a,b=p,x=x,y=z),list(a=l,b=p,x=y,y=z))) {
  rr<-fit(j$a,j$b,j$x,j$y,'three_trait_common_set_pair_check',ss=common);rr$row[,candidate_expression_dataset:=a$dataset];triple[[length(triple)+1L]]<-rr$row
 }
}
if(length(follow))fwrite(rbindlist(follow,fill=TRUE),file.path(O,'candidate_expression_PD.tsv'),sep='\t')
if(length(triple))fwrite(rbindlist(triple,fill=TRUE),file.path(O,'three_trait_checks_NOT_joint_posterior.tsv'),sep='\t')
sink(file.path(O,'sessionInfo.txt'));print(sessionInfo());print(.libPaths());sink()
capture.output(getFromNamespace('approx.bf.estimates','coloc'),file=file.path(O,'installed_ABF_formula.txt'))
capture.output(getFromNamespace('adjust_prior','coloc'),file=file.path(O,'installed_adjust_prior.txt'))
summary<-list(complete=TRUE,expression_lipid_tests=sum(main$gene_id!=''),lipid_PD_tests=sum(main$gene_id==''),primary_computed=sum(main$status=='computed'),
 sensitivity_rows=nrow(sens),candidate_review_rows=nrow(candidates),main_H4_ge_0_5=sum(main$gene_id!='' & main$H4>=.5,na.rm=TRUE),
 main_H4_ge_0_8=sum(main$gene_id!='' & main$H4>=.8,na.rm=TRUE),maximum_expression_lipid_H4=max(main[gene_id!='',H4],na.rm=TRUE),
 actual_prior_adjustment_rows=sum(main$p1_requested!=main$p1_used|main$p2_requested!=main$p2_used|main$p12_requested!=main$p12_used,na.rm=TRUE),
 source='coloc 5.2.3',joint_three_trait_posterior=FALSE,mediation=FALSE,independent_replication=FALSE)
write_json(summary,file.path(O,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE);print(summary);print(head(ranked[,.(region_id,gene_id,source1,H4,sensitivity_H4_max,classification,coverage_label)],20))
