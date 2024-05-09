
library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)


load("sims_N100_m2_s1_rtruncnorm13.rda")

#define groups to compare
col_group_in = c(rep("A", 19), rep("B",20))

source('miRglmm.R')

#find truth of each simulation
results_all=list()
for (ind_run in seq(1,length(sims))){ 
se=sims[[ind_run]]$sim_se
change_miRNA_up=sims[[ind_run]]$change_miRNA_up
change_miRNA_down=sims[[ind_run]]$change_miRNA_down
#find truth for group A effect
effect_up=sims[[ind_run]]$effect_up
iup = which(rowData(se)$miRNA %in% change_miRNA_up$miRNA)
true_FC=rep(1, dim(se)[1])
true_FC[iup]=1/effect_up
#find truth for group B effect
effect_down=sims[[ind_run]]$effect_down
idown = which(rowData(se)$miRNA %in% change_miRNA_down$miRNA)
true_FC[idown]=effect_down

true_FC=data.frame("sequence"=rowData(se)$uniqueSequence, "miRNA"=rowData(se)$miRNA, "true FC"=true_FC)
true_FC$true_logFC=log(true_FC$true.FC)

#DESeq2 on isomiRs
se$col_group=col_group_in
panelB_raw=SummarizedExperiment(assays=list(counts=as.matrix(assay(se))))
rowData(panelB_raw)=rowData(se)
panelB_raw$col_group=col_group_in
ddsSE=DESeqDataSet(panelB_raw, design=~col_group)
dds=DESeq(ddsSE)
res=results(dds)
true_FC$DESeq2_estlogFC=log(2^res$log2FoldChange)

#load sim results
load(file=paste0("sim results/sims_N100_m2_s1_rtruncnorm13_results", ind_run, ".rda"))

##### remove models where miRglmm not fit
double_warn=sapply(fits[["miRglmm"]], 'typeof')
fits[["miRglmm"]][which(double_warn=="double")]=NULL
double_warn=sapply(fits[["miRglmm_reduced"]], 'typeof')
fits[["miRglmm_reduced"]][which(double_warn=="double")]=NULL


#find isomiR-level point estimates
miRNA_vec=names(fits[["miRglmm"]])
for (ind in 1:length(miRNA_vec)){
  miRNA_in=miRNA_vec[ind]
  f1=fits[["miRglmm"]][[miRNA_in]]
  #pull out individual sequence estimates
  seq_effect=data.frame("sequence"=rownames(ranef(f1)$sequence), "est_logFC"=ranef(f1)$sequence$col_groupB+fixef(f1)[2])
  seq_effect$miRNA=miRNA_in
  if (ind==1){
  out=merge(true_FC, seq_effect, by=c("miRNA", "sequence"))
  results=out
  } else {
    out=merge(true_FC, seq_effect, by=c("miRNA", "sequence"))
    results=rbind(results, out)
  }
}

bias=data.frame(results[, c("est_logFC", "DESeq2_estlogFC")]-results$true_logFC)
squared_error=bias^2
MSE_sim=t(data.frame(MSE=colMeans(squared_error)))
squared_error$true_logFC=rep(0, dim(squared_error)[1])
squared_error$true_logFC[which(results$miRNA %in% change_miRNA_up$miRNA)]=log(0.5)
squared_error$true_logFC[which(results$miRNA %in% change_miRNA_down$miRNA)]=log(2)
MSE_sim_by_truth=squared_error %>% group_by(true_logFC) %>% summarise_all(funs(mean))


results_all[["MSE_sim"]][[ind_run]]=MSE_sim
results_all[["MSE_sim_by_truth"]][[ind_run]]=MSE_sim_by_truth
}

save(results_all, file="sim results/sims_N100_m2_s1_rtruncnorm13_seqDE_results.rda")

#MSE overall
MSE_mat=as.data.frame(matrix(unlist(results_all[["MSE_sim"]]), ncol=2, byrow=TRUE))
colnames(MSE_mat)=colnames(results_all[["MSE_sim"]][[1]])

#MSE by truth
out_miRglmm=data.frame(rbind(t(sapply(seq(1, length(sims)), function(row) results_all[["MSE_sim_by_truth"]][[row]][["est_logFC"]]))))
colnames(out_miRglmm)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
miRglmm_means=data.frame(t(colMeans(out_miRglmm)))
colnames(miRglmm_means)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
miRglmm_means$method="miRglmm"


out_DESeq2=data.frame(rbind(t(sapply(seq(1, length(sims)), function(row) results_all[["MSE_sim_by_truth"]][[row]][["DESeq2_estlogFC"]]))))
colnames(out_DESeq2)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
DESeq2_means=data.frame(t(colMeans(out_DESeq2)))
colnames(DESeq2_means)=results_all[["MSE_sim_by_truth"]][[1]][["true_logFC"]]
DESeq2_means$method="DESeq2"

MSE_by_truth=rbind(miRglmm_means, DESeq2_means)
