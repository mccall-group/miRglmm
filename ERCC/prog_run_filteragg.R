library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)
library(parallel)
library(doParallel)

load("ERCC/ERCC_filtered.rda")


#define groups to compare
col_group_in = colData(panel_B_filter)$Pool
adjust_var_in=colData(panel_B_filter)$Lab


source('miRglmnb.R')

ncores=8 #match to what was requested on BH
print(ncores)

#start cluster to run miRNA in parallel

if (ncores>1){
  cl=makeCluster(ncores)
  registerDoParallel(cl)
}

#calculate median CPM for each sequence
Y_all=as.data.frame(as.matrix(t(assays(panel_B_filter)$raw_counts)))
total_counts=rowSums(Y_all)
cpm=Y_all/(total_counts/1000000)
median_cpm=apply(cpm, 2, median)

fits_all_filters=list()

for (ind_filter in c(-1, -0.5, 0, 0.5, 1, 1.5, 2)){
  fits=list()
#filter before aggregating
idx=which(log(median_cpm)>ind_filter)
#aggregate data to miRNAs and fit miRglmnb
miRNA_counts = t(apply(assay(panel_B_filter)[idx,], 2, function(x) by(x, rowData(panel_B_filter)[idx,]$miRNA, sum)))
fits[["miRglmnb"]]= miRglmnb(miRNA_counts, col_group=col_group_in, ncores = ncores, adjust_var=adjust_var_in)



#run DESeq2 
panelB_raw=SummarizedExperiment(assays=list(t(miRNA_counts)), rowData=colnames(miRNA_counts), colData=cbind(col_group_in, adjust_var_in))
names(colData(panelB_raw))=c("col_group", "adjust_var")
ddsSE=DESeqDataSet(panelB_raw, design=~adjust_var+col_group)
fits[["DESeq2"]]=DESeq(ddsSE)

#run edgeR
edgeR_set=DGEList(counts=t(miRNA_counts), group=col_group_in)
design=model.matrix(~adjust_var_in+col_group_in)
edgeR_set=estimateDisp(edgeR_set, design)
et=glmQLFit(edgeR_set, design, coef="col_group_in")
idx=which(str_detect(colnames(design), "col_group"))
fits[["edgeR"]]=glmQLFTest(et, coef=colnames(design)[idx])


#run limma-voom
limvoom_set=DGEList(counts=t(miRNA_counts))
y=voom(limvoom_set, design) 
limvoom_fit=lmFit(y,design)
fits[["limvoom"]]=eBayes(limvoom_fit, trend=TRUE)

fits_all_filters[[as.character(ind_filter)]]=fits

}

if (ncores>1){
  stopCluster(cl)
}

#will save file containing list of 1 unfiltered, 7 filtered miRglmm model fits and 4 aggregated model fits 
save(fits_all_filters, file="ERCC/filter_before_agg_results.rda")


