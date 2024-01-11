library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)


source("process_results.R")
results=list()

#for each filter calculate beta for all methods
load(file="ERCC/all_filters_results.rda")
load(file="ERCC/filter_before_agg_results.rda")

filters=c("-1", "-0.5", "0", "0.5", "1", "1.5", "2")
filters_long=c("filter -1", "filter -0.5", "filter 0", "filter 0.5", "filter 1", "filter 1.5", "filter 2")

#combine filtered agg with filtered seq models to get table of betas for each model
beta_hat=list()
for (ind in seq(1, length(filters))){
  filter_in=filters[ind]
  filter_long_in=filters_long[ind]
  fits_agg_filtered=fits_all_filters[[filter_in]]
  fits_agg_filtered[["miRglmm"]]=fits[[filter_long_in]][["miRglmm"]]
  fits_agg_filtered[["miRglmm_reduced"]]=fits[[filter_long_in]][["miRglmm_reduced"]]
  beta_hat[[filter_in]]=get_betas(fits_agg_filtered,  var="col_group")
}

#pull unfiltered results in
pull=c("no filter", "miRglmnb", "DESeq2", "edgeR", "limvoom")
fits_in=fits[pull]
fits_in[["miRglmm"]]=fits_in[["no filter"]][["miRglmm"]]
fits_in[["miRglmm_reduced"]]=fits_in[["no filter"]][["miRglmm_reduced"]]
fits_in[["no filter"]]=NULL
#collect betas
beta_hat[["no filter"]]=get_betas(fits_in,  var="col_group")
  
save(beta_hat, file="ERCC/filter_agg_betahats.rda")