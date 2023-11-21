library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)

load("sims_N100_m2_s1.rda")
ind_run=1 #can be run in parallel using array on BH


#start cluster to run miRNA in parallel
ncores_use=8 #match to what was requested on BH
if (ncores>1){
  cl=makeCluster(ncores)
  registerDoParallel(cl)
}
#startTime=Sys.time()

#define groups to compare
col_group_in = c(rep("A", 19), rep("B",20))

#fit miRglmm full and reduced models
fits = miRglmm(sims[[ind_run]]$sim_se, col_group=col_group_in, ncores = ncores_use)

#aggregate data to miRNAs and fit miRglmnb
miRNA_counts = t(apply(assay(sims[[ind_run]]$sim_se), 2, function(x) by(x, rowData(sims[[ind_run]]$sim_se)$miRNA, sum)))
fits[["miRglmnb"]]= miRglmnb(miRNA_counts, col_group=col_group_in, ncores = ncores_use)

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


#will save file containing one list of 6 model fits for each simulated dataset
save(fits, file=paste0("sim results/sims_N100_m2_s1_results", ind_run, ".rda"))

