library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)


source("process_results.R")
results=list()

#for each filter calculate beta and SE (if applicable to method)

load(file="ERCC/all_filters_results.rda")

filters=c("filter -1", "filter -0.5", "filter 0", "filter 0.5", "filter 1", "filter 1.5", "filter 2", "no filter")

for (ind in seq(1, length(filters))){
filter_in=filters[ind]
if (ind==1){
pull=c(filter_in, "miRglmnb", "DESeq2", "edgeR", "limvoom")
fits_in=fits[pull]
fits_in[["miRglmm"]]=fits_in[[filter_in]][["miRglmm"]]
fits_in[["miRglmm_reduced"]]=fits_in[[filter_in]][["miRglmm_reduced"]]
fits_in[[filter_in]]=NULL
#collect betas
beta_hat=get_betas(fits_in,  var="col_group")

#collect SEs
SE_hat=get_SEs(fits_in,  var="col_group")

#collect number of sequences used in model
n_seq=get_nseq(fits_in[["miRglmm"]])
colnames(n_seq)=filter_in

#run LRT for random slope effect
LRTp=run_LRT(fits_in[["miRglmm"]], fits_in[["miRglmm_reduced"]])
colnames(LRTp)=filter_in

#miRglmm variance components
var_comp=list()
var_comp[[filter_in]]=get_varcomp(fits_in[["miRglmm"]], fits_in[["miRglmm_reduced"]])




#CI widths (Wald, 95% default)
CI_widths=getCI_widths(beta_hat, SE_hat)
colnames(CI_widths)[1]=filter_in
CI_LL=getCI_LL(beta_hat, SE_hat)
colnames(CI_LL)[1]=filter_in
CI_UL=getCI_UL(beta_hat, SE_hat)
colnames(CI_UL)[1]=filter_in

colnames(beta_hat)[1]=filter_in
colnames(SE_hat)[1]=filter_in




} else {
  pull=c(filter_in, "miRglmnb", "DESeq2", "edgeR", "limvoom")
  fits_in=fits[pull]
  fits_in[["miRglmm"]]=fits_in[[filter_in]][["miRglmm"]]
  fits_in[["miRglmm_reduced"]]=fits_in[[filter_in]][["miRglmm_reduced"]]
  fits_in[[filter_in]]=NULL
  
  #collect betas and merge
  out_beta=get_betas(fits_in,  var="col_group")
  out_SE=get_SEs(fits_in,  var="col_group")
  out=out_beta["miRglmm"]
  colnames(out)=filter_in
  beta_hat=transform(merge(beta_hat, out, by='row.names', all=T), row.names=Row.names, Row.names=NULL)
  
  #collect SEs and merge
  out=out_SE["miRglmm"]
  colnames(out)=filter_in
  SE_hat=transform(merge(SE_hat, out, by='row.names', all=T), row.names=Row.names, Row.names=NULL)
  
  #collect number of sequences used in model and merge
  out=get_nseq(fits_in[["miRglmm"]])
  colnames(out)=filter_in
  n_seq=transform(merge(n_seq, out, by='row.names', all=T), row.names=Row.names, Row.names=NULL)
  
  
  #run LRT for random slope effect
  out=run_LRT(fits_in[["miRglmm"]], fits_in[["miRglmm_reduced"]])
  colnames(out)=filter_in
  LRTp=transform(merge(LRTp, out, by='row.names', all=T), row.names=Row.names, Row.names=NULL)
  
  
  #miRglmm variance components
  var_comp=list()
  var_comp[[filter_in]]=get_varcomp(fits_in[["miRglmm"]], fits_in[["miRglmm_reduced"]])
  
  
  
  
  #CI widths (Wald, 95% default)
  out=getCI_widths(out_beta, out_SE)
  out=out["miRglmm"]
  colnames(out)[1]=filter_in
  CI_widths=transform(merge(CI_widths, out, by='row.names', all=T), row.names=Row.names, Row.names=NULL)
  
  out=getCI_LL(out_beta, out_SE)
  out=out["miRglmm"]
  colnames(out)[1]=filter_in
  CI_LL=transform(merge(CI_LL, out, by='row.names', all=T), row.names=Row.names, Row.names=NULL)
  
  
  out=getCI_UL(out_beta, out_SE)
  out=out["miRglmm"]
  colnames(out)[1]=filter_in
  CI_UL=transform(merge(CI_UL, out, by='row.names', all=T), row.names=Row.names, Row.names=NULL)
  
}
}

#create list of results
results[["beta_hat"]]=beta_hat
results[["SE_hat"]]=SE_hat
results[["LRTp"]]=LRTp
results[["var_comp"]]=var_comp
results[["CI_width"]]=CI_widths
results[["CI_LL"]]=CI_LL
results[["CI_UL"]]=CI_UL
results[["n_seq"]]=n_seq


save(results, file='ERCC/allfilters_processed_results.rda')
