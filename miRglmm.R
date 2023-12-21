library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(foreach)
library(doParallel)
miRglmm <- function(se, col_group = c(rep("A", 19), rep("B",20)),
                    min_med_lcpm = -1, ncores = 1, adjust_var=NA){
  
  if (is.na(adjust_var)){
    adjust_var=rep(NA, length(col_group))
  }
    
  ## for each miRNA (this could be parallelized)
  uniq_miRNA = unique(rowData(se)$miRNA)
  total_counts=colSums(assay(se))
  if (ncores==1){
    print('running non-parallel')
    f1_list=list()
    f1_sub_list=list()
  for(ind3 in seq(1, length(uniq_miRNA))){
    cat(uniq_miRNA[ind3], "\n")
    
    ## subset sequences that map to the miRNA
    idx = which(rowData(se)$miRNA == uniq_miRNA[ind3])
    Y_all_sub = t(assay(se)[idx, ])
    Y_seq_labels_sub = rowData(se)$uniqueSequence[idx]
    
    ## compute median CPM for each sequence
    cpm_sub = Y_all_sub / (total_counts/1e6)
    median_cpm_sub = apply(cpm_sub, 2, median)

    ## filter sequences with median CPM less than the min_med_lcpm threshold
    run_miRNA=1
    if(!is.null(min_med_lcpm)){
      idx2 = which(log(median_cpm_sub) > min_med_lcpm)
      Y_all_sub = Y_all_sub[ ,idx2]
      Y_seq_labels_sub = Y_seq_labels_sub[idx2]
      if (length(idx2)<2){
        run_miRNA=0
      }
    }
    f1=0
    f1_sub=0
    if (run_miRNA==1){
    ## format data to fit glmer.nb models
    data_wide = as.data.frame(as.matrix(Y_all_sub))
    colnames(data_wide) = Y_seq_labels_sub
    sample_labels = rownames(Y_all_sub)
    data_wide = cbind(col_group, adjust_var, total_counts, sample_labels, data_wide)
    data_long = melt(data_wide, 
                     id.vars = c("col_group", "adjust_var", "total_counts", "sample_labels"),
                     variable.name = "sequence", 
                     value.name = "count")
    if (all(is.na(adjust_var)==TRUE)){
      Formula_full="count ~col_group + offset(log(total_counts/1e4)) + (1+col_group|sequence) + (1|sample_labels)"
      Formula_sub="count ~col_group + offset(log(total_counts/1e4)) + (1|sequence) + (1|sample_labels)"
    } else {
      Formula_full="count ~col_group + adjust_var + offset(log(total_counts/1e4)) + (1+col_group|sequence) + (1|sample_labels)"
      Formula_sub="count ~col_group + adjust_var + offset(log(total_counts/1e4)) + (1|sequence) + (1|sample_labels)"
    }
    tryCatch({
      f1 = glmer.nb(Formula_full, 
                  data=data_long, 
                  control=(glmerControl(optimizer="bobyqa", 
                                        tolPwrss = 1e-3, 
                                        optCtrl=list(maxfun=2e5))))
    }, error = function(e){cat("ERROR :", uniq_miRNA[ind3])})
    tryCatch({  
    f1_sub = glmer.nb(Formula_sub, 
                        data=data_long, 
                        control=(glmerControl(optimizer="bobyqa", 
                                              tolPwrss = 1e-3, 
                                              optCtrl=list(maxfun=2e5))))
  }, error = function(e){cat("ERROR :", uniq_miRNA[ind3])})
    f1_list[[uniq_miRNA[ind3]]]=f1
    f1_sub_list[[uniq_miRNA[ind3]]]=f1_sub
    }
  }
    fits=list()
    fits[["miRglmm"]]=f1_list
    fits[["miRglmm_reduced"]]=f1_sub_list
  } else {
    print('running parallel')
    fits_full=foreach(ind3=1:length(uniq_miRNA), .packages=c("tidyverse", "reshape2", "lme4", "SummarizedExperiment")) %dopar% {
      #cat(uniq_miRNA[ind3], "\n")
      f1_list=list()
      ## subset sequences that map to the miRNA
      idx = which(rowData(se)$miRNA == uniq_miRNA[ind3])
      Y_all_sub = t(assay(se)[idx, ])
      Y_seq_labels_sub = rowData(se)$uniqueSequence[idx]
      
      ## compute median CPM for each sequence
      cpm_sub = Y_all_sub / (total_counts/1e6)
      median_cpm_sub = apply(cpm_sub, 2, median)
      
      ## filter sequences with median CPM less than the min_med_lcpm threshold
      run_miRNA=1
      if(!is.null(min_med_lcpm)){
        idx2 = which(log(median_cpm_sub) > min_med_lcpm)
        Y_all_sub = Y_all_sub[ ,idx2]
        Y_seq_labels_sub = Y_seq_labels_sub[idx2]
        if (length(idx2)<2){
          run_miRNA=0
        }
      }
      if (run_miRNA==1){
      
      ## format data to fit glmer.nb models
      data_wide = as.data.frame(as.matrix(Y_all_sub))
      colnames(data_wide) = Y_seq_labels_sub
      sample_labels = rownames(Y_all_sub)
      data_wide = cbind(col_group, adjust_var, total_counts, sample_labels, data_wide)
      data_long = melt(data_wide, 
                       id.vars = c("col_group", "adjust_var", "total_counts", "sample_labels"),
                       variable.name = "sequence", 
                       value.name = "count")
      f1=0
      if (all(is.na(adjust_var)==TRUE)){
        Formula_full="count ~col_group + offset(log(total_counts/1e4)) + (1+col_group|sequence) + (1|sample_labels)"
         } else {
        Formula_full="count ~col_group + adjust_var + offset(log(total_counts/1e4)) + (1+col_group|sequence) + (1|sample_labels)"
       }
      tryCatch({
        f1 = glmer.nb(Formula_full, 
                      data=data_long, 
                      control=(glmerControl(optimizer="bobyqa", 
                                            tolPwrss = 1e-3, 
                                            optCtrl=list(maxfun=2e5))))
      }, error = function(e){cat("ERROR :", uniq_miRNA[ind3])})

      f1_list[[uniq_miRNA[ind3]]]=f1
      }
      return(f1_list)
    }
    fits_red=foreach(ind3=1:length(uniq_miRNA), .packages=c("tidyverse", "reshape2", "lme4", "SummarizedExperiment")) %dopar% {
      #cat(uniq_miRNA[ind3], "\n")
      f1_list=list()
      ## subset sequences that map to the miRNA
      idx = which(rowData(se)$miRNA == uniq_miRNA[ind3])
      Y_all_sub = t(assay(se)[idx, ])
      Y_seq_labels_sub = rowData(se)$uniqueSequence[idx]
      
      ## compute median CPM for each sequence
      cpm_sub = Y_all_sub / (total_counts/1e6)
      median_cpm_sub = apply(cpm_sub, 2, median)
      
      ## filter sequences with median CPM less than the min_med_lcpm threshold
      run_miRNA=1
      if(!is.null(min_med_lcpm)){
        idx2 = which(log(median_cpm_sub) > min_med_lcpm)
        Y_all_sub = Y_all_sub[ ,idx2]
        Y_seq_labels_sub = Y_seq_labels_sub[idx2]
        if (length(idx2)<2){
          run_miRNA=0
        }
      }
      if (run_miRNA==1){
      
      ## format data to fit glmer.nb models
      data_wide = as.data.frame(as.matrix(Y_all_sub))
      colnames(data_wide) = Y_seq_labels_sub
      sample_labels = rownames(Y_all_sub)
      data_wide = cbind(col_group, adjust_var, total_counts, sample_labels, data_wide)
      data_long = melt(data_wide, 
                       id.vars = c("col_group", "adjust_var", "total_counts", "sample_labels"),
                       variable.name = "sequence", 
                       value.name = "count")
      f1_sub=0
      if (all(is.na(adjust_var)==TRUE)){
        Formula_sub="count ~col_group + offset(log(total_counts/1e4)) + (1|sequence) + (1|sample_labels)"
      } else {
        Formula_sub="count ~col_group + adjust_var + offset(log(total_counts/1e4)) + (1|sequence) + (1|sample_labels)"
      }
      tryCatch({

        f1_sub = glmer.nb(Formula_sub, 
                          data=data_long, 
                          control=(glmerControl(optimizer="bobyqa", 
                                                tolPwrss = 1e-3, 
                                                optCtrl=list(maxfun=2e5))))
        
      }, error = function(e){cat("ERROR :", uniq_miRNA[ind3])})
      f1_list[[uniq_miRNA[ind3]]]=f1_sub
      }
      return(f1_list)
    }
  fits=list()
  fits[["miRglmm"]]=unlist(fits_full)
  fits[["miRglmm_reduced"]]=unlist(fits_red)
  }

  return(fits)
}


