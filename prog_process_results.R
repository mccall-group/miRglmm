
#for each results file calculate beta and SE (if applicable to method)
ind_run=1 
load(file=paste0("sim results/sims_N100_m2_s1_results", ind_run, ".rda"))

#find truth based on ind_run
uniq_miRNA=unique(rowData(sims[[ind_run]]$sim_se)$miRNA)
true_logFC=rep(log(1), length(uniq_miRNA))
true_logFC[which(uniq_miRNA %in% as.character(sims[[ind_run]]$change_miRNA_up$miRNA))]=log(2)
true_logFC[which(uniq_miRNA %in% as.character(sims[[ind_run]]$change_miRNA_down$miRNA))]=log(0.5)
true_logFC_BvsA=data.frame("true_beta"=true_logFC*-1)
rownames(true_logFC_BvsA)=uniq_miRNA

#collect betas
beta_hat=get_betas(fits,  var="col_group")

#merge with truth
beta_hat=transform(merge(beta_hat, true_logFC_BvsA, by='row.names'), row.names=Row.names, Row.names=NULL)

#collect SEs
SE_hat=get_SEs(fits,  var="col_group")

#run LRT for random slope effect
LRTp=run_LRT(fits[["miRglmm"]], fits[["miRglmm_reduced"]])