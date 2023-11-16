library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)

aggseq <- function(se){
  uniq_miRNA = unique(rowData(se)$miRNA)
  for (ind3 in seq(1, length(uniq_miRNA))){
  uniq_miRNA_in=uniq_miRNA[ind3]
  idx = which(rowData(se)$miRNA == uniq_miRNA[ind3])
  Y_all_sub = t(assay(se)[idx, ])
  Y_seq_labels_sub = rowData(se)$uniqueSequence[idx]
  
  
  data_wide=as.data.frame(as.matrix(Y_all_sub))
  data_wide=as.data.frame(rowSums(data_wide))
  colnames(data_wide)=uniq_miRNA_in
  if (ind3==1){
    data_miRNA=data_wide
  } else {
    data_miRNA=cbind(data_miRNA, data_wide)
  }
  }
  return(data_miRNA)
}

miRDESeq2 <- function(agg_data, col_group = c(rep("A", 18), rep("B",19))){
  library(DESeq2)
  panelB_raw=SummarizedExperiment(assays=list(counts=t(agg_data)), rowData=colnames(agg_data), colData=col_group)
  names(colData(panelB_raw))=c("Pool")
  ddsSE=DESeqDataSet(panelB_raw, design=~Pool)
  dds=DESeq(ddsSE)
  res=results(dds)
  
  deseq2_results=as.data.frame(cbind(res$log2FoldChange, res$lfcSE)) ###note this SE is on log2 scale
  rownames(deseq2_results)=rownames(res)
  colnames(deseq2_results)=c("est_log2FC", "PoolB SE")
  deseq2_results$col_groupB=log(2^deseq2_results$est_log2FC)
  return(deseq2_results)
}

miRedgeR <- function(agg_data, col_group = c(rep("A", 18), rep("B",19))){
  
  library(edgeR)
  edgeR_set=DGEList(counts=t(agg_data), group=col_group)
  design=model.matrix(~col_group)
  edgeR_set=estimateDisp(edgeR_set, design)
  et=glmQLFit(edgeR_set, design)
  et_test=glmQLFTest(et)
  edgeR_results=as.data.frame(cbind(log(2^et_test$table$logFC), et_test$table$PValue))
  rownames(edgeR_results)=rownames(et_test$table)
  colnames(edgeR_results)=c("col_groupB", "Pval")
  
  return(edgeR_results)
}

miRlimvoom <- function(agg_data, col_group = c(rep("A", 18), rep("B",19))){
  library(edgeR)
  col_data=as.data.frame(col_group)
  limvoom_set=DGEList(counts=t(agg_data))
  design=model.matrix(~col_group)
  y=voom(limvoom_set, design)
  fit=lmFit(y, design)
  eb_fit=eBayes(fit, trend=TRUE)
  coeff_out=log(2^coef(eb_fit)[,colnames(coef(eb_fit))=="col_groupB"])
  coeff_out_log2FC=coef(eb_fit)[,colnames(coef(eb_fit))=="col_groupB"]
  SE_out=sqrt(eb_fit$s2.post)*eb_fit$stdev.unscaled[,colnames(coef(eb_fit))=="col_groupB"]
  p_out=eb_fit$p.value[,colnames(eb_fit$p.value)=="col_groupB"]
  limmavoom_results=as.data.frame(cbind(coeff_out, coeff_out_log2FC, SE_out, p_out))
  rownames(limmavoom_results)=rownames(coef(eb_fit))
  colnames(limmavoom_results)=c("col_groupB", "est_log2FC", "SE", "Pval")
  
  return(limmavoom_results)
}

miRNBglm <- function(agg_data, col_group = c(rep("A", 18), rep("B",19))){
  library(MASS)
  uniq_miRNA=unique(colnames(agg_data))
  total_counts=rowSums(agg_data)
  for (ind3 in seq(1, length(uniq_miRNA))){
    uniq_miRNA_in=uniq_miRNA[ind3]
    idx = which(colnames(agg_data) == uniq_miRNA[ind3])
    Y_all_sub = data_miRNA[, idx]
    
    
    data_wide=as.data.frame(as.matrix(Y_all_sub))
    colnames(data_wide)="count"
    data_wide=cbind(col_group,total_counts, data_wide)
    
    
    tryCatch({
      f1=glm.nb(count~col_group+offset(log(total_counts)), data=data_wide, control=glm.control(maxit=1000))
      out_f1_coef=as.data.frame(t(coef(f1)))
      out_f1_SE=as.data.frame(t(summary(f1)$coefficients[,2]))
      out_f1_SE=cbind(out_f1_SE, as.data.frame(t(suppressMessages(confint(f1)[2,]))))
      colnames(out_f1_SE)=c("Int SE", "PoolB SE", "PoolB LL", "PoolB UL")
      out_f1_coef=cbind(out_f1_coef, out_f1_SE)
      out_f1_coef$miRNA=uniq_miRNA_in
      out_f1_coef$theta=summary(f1)$theta
      out_f1_coef$medianCPM=median(data_wide$cpm)
      if (ind3==1){
        miRNA_results=out_f1_coef
      } else {
        miRNA_results=rbind(miRNA_results, out_f1_coef)
      }
    }, error=function(e){cat("ERROR :", uniq_miRNA_in)})
    
  }
  
  return(miRNA_results)
}

data_miRNA=aggseq(sims[[1]]$sim_se)
tst = miRDESeq2(data_miRNA)
tst2 = miRedgeR(data_miRNA)
tst3 = miRlimvoom(data_miRNA)
tst4 = miRNBglm(data_miRNA)

uniq_miRNA = unique(rowData(sims[[1]]$sim_se)$miRNA)
true_logFC = rep(log(1), length(uniq_miRNA))
true_logFC[which(uniq_miRNA %in% as.character(sims[[1]]$change_miRNA_up$miRNA))] = log(2)
true_logFC[which(uniq_miRNA %in% as.character(sims[[1]]$change_miRNA_down$miRNA))] = log(0.5)
true_logFC_BvsA = true_logFC*-1
tst$true_logFC=true_logFC_BvsA

plot(x=true_logFC_BvsA, y=tst$col_groupB)
