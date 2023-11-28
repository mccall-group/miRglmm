library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)

#load datasets to grab truth from
load("sims_N100_m2_s1.rda")
results=list()
#for each results file calculate beta and SE (if applicable to method)
for (ind_run in seq(1,10)){ 
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

#miRglmm variance components

var_comp=get_varcomp(fits[["miRglmm"]], fits[["miRglmm_reduced"]])

#bias
bias=data.frame(beta_hat[, c("miRglmm","miRglmnb","DESeq2","edgeR","limmavoom")]-beta_hat$true_beta)
squared_error=bias^2
MSE_sim=t(data.frame(MSE=colMeans(squared_error)))
MSE_sim_by_truth=squared_error %>% group_by(true_beta) %>% summarise_all(funs(mean))

bias$true_beta=beta_hat$true_beta
squared_error$true_beta=bias$true_beta


#CI widths (Wald, 95% default)
CI_widths=getCI_widths(beta_hat, SE_hat)

#coverage indicators (95% Wald CI based)
coverage_indicators=getCoverageInd(beta_hat, SE_hat)
coverage_probability_sim=t(data.frame("coverage_prob_sim"=colSums(coverage_indicators)/dim(coverage_indicators)[1]))

coverage_indicators=transform(merge(coverage_indicators, true_logFC_BvsA, by='row.names'), row.names=Row.names, Row.names=NULL)
coverage_probability_sim_by_truth=coverage_indicators %>% group_by(true_beta) %>% summarise_all(funs(mean))


#create list of results
results[["beta_hat"]][[ind_run]]=beta_hat
results[["SE_hat"]][[ind_run]]=SE_hat
results[["LRTp"]][[ind_run]]=LRTp
results[["var_comp"]][[ind_run]]=var_comp
results[["squared_error"]][[ind_run]]=squared_error
results[["bias"]][[ind_run]]=bias
results[["CI_width"]][[ind_run]]=CI_widths
results[["coverage_indicator"]][[ind_run]]=coverage_indicators
results[["MSE_sim"]][[ind_run]]=MSE_sim
results[["MSE_sim_by_truth"]][[ind_run]]=MSE_sim_by_truth
results[["coverage_prob_sim"]][[ind_run]]=coverage_probability_sim
results[["coverage_prob_sim_by_truth"]][[ind_run]]=coverage_probability_sim_by_truth
}

miRglmm_squared_error=sapply(1:length(results[["squared_error"]]), function(i) cbind(results[["squared_error"]][[i]]$miRglmm))
rownames(miRglmm_squared_error)=rownames(results[["squared_error"]][[1]])
MSE_miRNA=data.frame("miRglmm"=rowMeans(miRglmm_squared_error))