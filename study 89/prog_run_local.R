library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)
library(parallel)
library(doParallel)

load("study89_data_subset_filtered2.rda")
ncores=8#as.integer(commandArgs(trailingOnly = TRUE)[1])#ncores=8 #match to what was requested on BH
print(ncores)

#start cluster to run miRNA in parallel

if (ncores>1){
  cl=makeCluster(ncores)
  registerDoParallel(cl)
}
startTime=Sys.time()

#define groups to compare
col_group_in = as.factor(colData(study89_data_subset_filtered2)$cell_tissue)

source('/scratch/aberry3/McCall/miRNA/miRglmm/miRglmm.R')
source('/scratch/aberry3/McCall/miRNA/miRglmm/miRglmnb.R')


#fit miRglmm full and reduced models
fits = miRglmm(study89_data_subset_filtered2, col_group=col_group_in, ncores = ncores)

#aggregate data to miRNAs and fit miRglmnb
miRNA_counts = t(apply(assay(study89_data_subset_filtered2), 2, function(x) by(x, rowData(study89_data_subset_filtered2)$miRNA, sum)))
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
fits[["edgeR"]]=glmQLFTest(et, contrast=contrasts_in)
fits[["edgeR_et"]]=et


#run limma-voom
limvoom_set=DGEList(counts=t(miRNA_counts))
y=voom(limvoom_set, design) #same design as edgeR
limvoom_fit=lmFit(y,design)
fits[["limvoom"]]=eBayes(limvoom_fit, trend=TRUE)


#will save file containing one list of 6 model fits for each simulated dataset
save(fits, file="results/filter_neg1_results.rda")
stopTime=Sys.time()
