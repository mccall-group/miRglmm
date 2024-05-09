library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)

#load datasets to grab truth from
load("study 89/study89_data_subset_filtered2.rda")
source("process_results.R")
results=list()
#for each results file calculate beta and SE (if applicable to method)
start_time=Sys.time()


load(file="study 89/results/filter_neg1_results.rda")





#collect betas
beta_hat_monocyte_vs_Blymph=get_betas(fits,  var="Monocyte")
beta_hat_NK_vs_Blymph=get_betas(fits,  var="Natural_Killer")
beta_hat_TCD4_vs_Blymph=get_betas(fits,  var="T_lymphocyte_CD4")
beta_hat_TCD8_vs_Blymph=get_betas(fits,  var="T_lymphocyte_CD8")

#collect SEs
SE_hat_monocyte_vs_Blymph=get_SEs(fits,  var="Monocyte")
SE_hat_NK_vs_Blymph=get_SEs(fits,  var="Natural_Killer")
SE_hat_TCD4_vs_Blymph=get_SEs(fits,  var="T_lymphocyte_CD4")
SE_hat_TCD8_vs_Blymph=get_SEs(fits,  var="T_lymphocyte_CD8")

#run LRT for random slope effect
LRTp=run_LRT(fits[["miRglmm"]], fits[["miRglmm_reduced"]])

#miRglmm variance components
#NEED TO ADAPT THIS CODE FOR >2 GROUP CATEGORICAL RANDOM EFFECT
#var_comp=get_varcomp(fits[["miRglmm"]], fits[["miRglmm_reduced"]])




#CI widths (Wald, 95% default)
CI_widths_monocyte=getCI_widths(beta_hat_monocyte_vs_Blymph, SE_hat_monocyte_vs_Blymph)
CI_LL_monocyte=getCI_LL(beta_hat_monocyte_vs_Blymph, SE_hat_monocyte_vs_Blymph)
CI_UL_monocyte=getCI_UL(beta_hat_monocyte_vs_Blymph, SE_hat_monocyte_vs_Blymph)

CI_widths_NK=getCI_widths(beta_hat_NK_vs_Blymph, SE_hat_NK_vs_Blymph)
CI_LL_NK=getCI_LL(beta_hat_NK_vs_Blymph, SE_hat_NK_vs_Blymph)
CI_UL_NK=getCI_UL(beta_hat_NK_vs_Blymph, SE_hat_NK_vs_Blymph)

CI_widths_TCD4=getCI_widths(beta_hat_TCD4_vs_Blymph, SE_hat_TCD4_vs_Blymph)
CI_LL_TCD4=getCI_LL(beta_hat_TCD4_vs_Blymph, SE_hat_TCD4_vs_Blymph)
CI_UL_TCD4=getCI_UL(beta_hat_TCD4_vs_Blymph, SE_hat_TCD4_vs_Blymph)

CI_widths_TCD8=getCI_widths(beta_hat_TCD8_vs_Blymph, SE_hat_TCD8_vs_Blymph)
CI_LL_TCD8=getCI_LL(beta_hat_TCD8_vs_Blymph, SE_hat_TCD8_vs_Blymph)
CI_UL_TCD8=getCI_UL(beta_hat_TCD8_vs_Blymph, SE_hat_TCD8_vs_Blymph)

beta_hat_monocyte_vs_Blymph$contrast="Monocyte vs B lymphocyte"
beta_hat_NK_vs_Blymph$contrast="Natural Killer vs B lymphocyte"
beta_hat_TCD4_vs_Blymph$contrast="CD4 T lymphocyte vs B lymphocyte"
beta_hat_TCD8_vs_Blymph$contrast="CD8 T lymphocyte vs B lymphocyte"
beta_hat=rbind(beta_hat_monocyte_vs_Blymph, beta_hat_NK_vs_Blymph, beta_hat_TCD4_vs_Blymph, beta_hat_TCD8_vs_Blymph)

SE_hat_monocyte_vs_Blymph$contrast="Monocyte vs B lymphocyte"
SE_hat_NK_vs_Blymph$contrast="Natural Killer vs B lymphocyte"
SE_hat_TCD4_vs_Blymph$contrast="CD4 T lymphocyte vs B lymphocyte"
SE_hat_TCD8_vs_Blymph$contrast="CD8 T lymphocyte vs B lymphocyte"
SE_hat=rbind(SE_hat_monocyte_vs_Blymph, SE_hat_NK_vs_Blymph, SE_hat_TCD4_vs_Blymph, SE_hat_TCD8_vs_Blymph)

CI_widths_monocyte$contrast="Monocyte vs B lymphocyte"
CI_widths_NK$contrast="Natural Killer vs B lymphocyte"
CI_widths_TCD4$contrast="CD4 T lymphocyte vs B lymphocyte"
CI_widths_TCD8$contrast="CD8 T lymphocyte vs B lymphocyte"
CI_widths=rbind(CI_widths_monocyte, CI_widths_NK, CI_widths_TCD4, CI_widths_TCD8)

CI_LL_monocyte$contrast="Monocyte vs B lymphocyte"
CI_LL_NK$contrast="Natural Killer vs B lymphocyte"
CI_LL_TCD4$contrast="CD4 T lymphocyte vs B lymphocyte"
CI_LL_TCD8$contrast="CD8 T lymphocyte vs B lymphocyte"
CI_LL=rbind(CI_LL_monocyte, CI_LL_NK, CI_LL_TCD4, CI_LL_TCD8)

CI_UL_monocyte$contrast="Monocyte vs B lymphocyte"
CI_UL_NK$contrast="Natural Killer vs B lymphocyte"
CI_UL_TCD4$contrast="CD4 T lymphocyte vs B lymphocyte"
CI_UL_TCD8$contrast="CD8 T lymphocyte vs B lymphocyte"
CI_UL=rbind(CI_UL_monocyte, CI_UL_NK, CI_UL_TCD4, CI_UL_TCD8)

#create list of results
results[["beta_hat"]]=beta_hat
results[["SE_hat"]]=SE_hat
results[["LRTp"]]=LRTp
#results[["var_comp"]]=var_comp
results[["CI_width"]]=CI_widths
results[["CI_LL"]]=CI_LL
results[["CI_UL"]]=CI_UL

end_time=Sys.time()
end_time-start_time

save(results, file="study 89/results/filter_neg1_processed_results.rda")
