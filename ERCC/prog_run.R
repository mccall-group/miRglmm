library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)
library(parallel)
library(doParallel)

load("ERCC/ERCC_filtered.rda")
ncores=8
print(ncores)

#start cluster to run miRNA in parallel
#ncores=8 #match to what was requested on BH
if (ncores>1){
  cl=makeCluster(ncores)
  registerDoParallel(cl)
}
#startTime=Sys.time()

#define groups to compare
col_group_in = colData(panel_B_filter)$Pool
adjust_var_in=colData(panel_B_filter)$Lab

source('miRglmm.R')
source('miRglmnb.R')

#fit miRglmm full and reduced models
fits[["filter -1"]] = miRglmm(panel_B_filter, col_group=col_group_in, ncores = ncores, min_med_lcpm = -1, adjust_var=adjust_var_in)
fits[["filter -0.5"]] = miRglmm(panel_B_filter, col_group=col_group_in, ncores = ncores, min_med_lcpm = -0.5, adjust_var=adjust_var_in)
fits[["filter 0"]] = miRglmm(panel_B_filter, col_group=col_group_in, ncores = ncores, min_med_lcpm = 0, adjust_var=adjust_var_in)
fits[["filter 0.5"]] = miRglmm(panel_B_filter, col_group=col_group_in, ncores = ncores, min_med_lcpm = 0.5, adjust_var=adjust_var_in)
fits[["filter 1"]] = miRglmm(panel_B_filter, col_group=col_group_in, ncores = ncores, min_med_lcpm = 1, adjust_var=adjust_var_in)
fits[["filter 1.5"]] = miRglmm(panel_B_filter, col_group=col_group_in, ncores = ncores, min_med_lcpm = 1.5, adjust_var=adjust_var_in)
fits[["filter 2"]] = miRglmm(panel_B_filter, col_group=col_group_in, ncores = ncores, min_med_lcpm = 2, adjust_var=adjust_var_in)

fits[["no filter"]] = miRglmm(panel_B_filter, col_group=col_group_in, ncores = ncores, min_med_lcpm = NULL, adjust_var=adjust_var_in)

#aggregate data to miRNAs and fit miRglmnb
miRNA_counts = t(apply(assay(panel_B_filter), 2, function(x) by(x, rowData(panel_B_filter)$miRNA, sum)))
fits[["miRglmnb"]]= miRglmnb(miRNA_counts, col_group=col_group_in, ncores = ncores)

if (ncores>1){
  stopCluster(cl)
}

#run DESeq2 
panelB_raw=SummarizedExperiment(assays=list(t(miRNA_counts)), rowData=colnames(miRNA_counts), colData=col_group_in)
names(colData(panelB_raw))=c("col_group")
ddsSE=DESeqDataSet(panelB_raw, design=~col_group)
fits[["DESeq2"]]=DESeq(ddsSE)

#run edgeR
edgeR_set=DGEList(counts=t(miRNA_counts), group=col_group_in)
design=model.matrix(~col_group_in)
edgeR_set=estimateDisp(edgeR_set, design)
et=glmQLFit(edgeR_set, design)
fits[["edgeR"]]=glmQLFTest(et)


#run limma-voom
limvoom_set=DGEList(counts=t(miRNA_counts))
y=voom(limvoom_set, design) #same design as edgeR
limvoom_fit=lmFit(y,design)
fits[["limvoom"]]=eBayes(limvoom_fit, trend=TRUE)


#will save file containing list of 1 unfiltered, 7 filtered miRglmm model fits and 4 aggregated model fits 
save(fits, file="ERCC/all_filters_results.rda")

