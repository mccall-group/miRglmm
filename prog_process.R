library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)

#load datasets to grab truth from
load("bladder_testes_data_subset_filtered2.rda")
source("process_results.R")
results=list()
#for each results file calculate beta and SE (if applicable to method)
start_time=Sys.time()


load(file="bladder_testes_results/filter_neg1_results.rda")





#collect betas
beta_hat=get_betas(fits,  var="col_group")


#collect SEs
SE_hat=get_SEs(fits,  var="col_group")

#collect pvalues for coefficient of interest
pvals=get_pvals(fits,  var="col_group")

#run LRT for random slope effect
LRTp=run_LRT(fits[["miRglmm"]], fits[["miRglmm_reduced"]])

#miRglmm variance components

var_comp=get_varcomp(fits[["miRglmm"]], fits[["miRglmm_reduced"]])




#CI widths (Wald, 95% default)
CI_widths=getCI_widths(beta_hat, SE_hat)
CI_LL=getCI_LL(beta_hat, SE_hat)
CI_UL=getCI_UL(beta_hat, SE_hat)



#create list of results
results[["beta_hat"]]=beta_hat
results[["SE_hat"]]=SE_hat
results[["LRTp"]]=LRTp
results[["var_comp"]]=var_comp
results[["CI_width"]]=CI_widths
results[["CI_LL"]]=CI_LL
results[["CI_UL"]]=CI_UL
results[["pvals"]]=pvals

end_time=Sys.time()
end_time-start_time

save(results, file="bladder_testes_results/filter_neg1_processed_results.rda")
