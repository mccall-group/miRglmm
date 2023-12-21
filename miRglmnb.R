library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(foreach)
library(doParallel)
miRglmnb <- function(agg_data, col_group = c(rep("A", 19), rep("B",20)), ncores = 1, adjust_var=NA){
  library(MASS)
  if (is.na(adjust_var)){
    adjust_var=rep(NA, length(col_group))
  }
  uniq_miRNA=unique(colnames(agg_data))
  total_counts=rowSums(agg_data)
  if (ncores==1){
    print('running non-parallel')
    f1_list=list()
  for (ind3 in seq(1, length(uniq_miRNA))){
    uniq_miRNA_in=uniq_miRNA[ind3]
    idx = which(colnames(agg_data) == uniq_miRNA[ind3])
    Y_all_sub = agg_data[, idx]
    
    
    data_wide=as.data.frame(as.matrix(Y_all_sub))
    colnames(data_wide)="count"
    data_wide=cbind(col_group, adjust_var, total_counts, data_wide)
    if (all(is.na(adjust_var)==TRUE)){
      Formula="count ~col_group + offset(log(total_counts/1e4))"
     
    } else {
      Formula="count ~col_group + adjust_var + offset(log(total_counts/1e4))"
      
    }
    f1=0
    tryCatch({
      f1=glm.nb(Formula, data=data_wide, control=glm.control(maxit=1000))
    }, error=function(e){cat("ERROR :", uniq_miRNA_in)})
    f1_list[[uniq_miRNA[ind3]]]=f1
  }
  } else {
    print('running parallel')
    fits=foreach(ind3=1:length(uniq_miRNA), .packages=c("tidyverse", "reshape2", "MASS", "SummarizedExperiment")) %dopar% {
      uniq_miRNA_in=uniq_miRNA[ind3]
      idx = which(colnames(agg_data) == uniq_miRNA[ind3])
      Y_all_sub = agg_data[, idx]
      
      
      data_wide=as.data.frame(as.matrix(Y_all_sub))
      colnames(data_wide)="count"
      data_wide=cbind(col_group,adjust_var, total_counts, data_wide)
      if (all(is.na(adjust_var)==TRUE)){
        Formula="count ~col_group + offset(log(total_counts/1e4))"
        
      } else {
        Formula="count ~col_group + adjust_var + offset(log(total_counts/1e4))"
        
      }
      f1=0
      tryCatch({
        f1=glm.nb(Formula, data=data_wide, control=glm.control(maxit=1000))
      }, error=function(e){cat("ERROR :", uniq_miRNA_in)})
      f1_list_in=list()
      f1_list_in[[uniq_miRNA[ind3]]]=f1
      return(f1_list_in)
    }
    f1_list=list()
   f1_list=unlist(fits, recursive=FALSE)
  }
  return(f1_list)
}


