library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
miRglmm <- function(se,  col_group = c(rep("A", 19), rep("B",20)), ncores=1,
                    min_med_lcpm = -1){
  
  ## for each miRNA (this could be parallelized)
  uniq_miRNA = unique(rowData(se)$miRNA)
  total_counts=colSums(assay(se))
  if (ncores==1){
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
    if(!is.null(min_med_lcpm)){
      idx2 = which(log(median_cpm_sub) > min_med_lcpm)
      Y_all_sub = Y_all_sub[ ,idx2]
      Y_seq_labels_sub = Y_seq_labels_sub[idx2]
    }
    
    ## format data to fit glmer.nb models
    data_wide = as.data.frame(as.matrix(Y_all_sub))
    colnames(data_wide) = Y_seq_labels_sub
    sample_labels = rownames(Y_all_sub)
    data_wide = cbind(col_group, total_counts, sample_labels, data_wide)
    data_long = melt(data_wide, 
                     id.vars = c("col_group", "total_counts", "sample_labels"),
                     variable.name = "sequence", 
                     value.name = "count")
    tryCatch({
      f1 = glmer.nb(count ~ col_group + offset(log(total_counts/1e6)) + 
                    (1+col_group|sequence) + (1|sample_labels), 
                  data=data_long, 
                  control=(glmerControl(optimizer="bobyqa", 
                                        tolPwrss = 1e-3, 
                                        optCtrl=list(maxfun=2e5))))
      ## a bunch of these produce the message:
      ## boundary (singular) fit: see help('isSingular')
      ## could we substitute f1 with f1_sub if we get this message?
      f1_sub = glmer.nb(count ~ col_group + offset(log(total_counts/1e6)) + 
                          (1|sequence) + (1|sample_labels), 
                        data=data_long, 
                        control=(glmerControl(optimizer="bobyqa", 
                                              tolPwrss = 1e-3, 
                                              optCtrl=list(maxfun=2e5))))
    }, error = function(e){cat("ERROR :", uniq_miRNA[ind3])})
    f1_list[[uniq_miRNA[ind3]]]=f1
    f1_sub_list[[uniq_miRNA[ind3]]]=f1_sub
  }
    fits=list()
    fits[["miRglmm"]]=f1_list
    fits[["miRglmm_reduced"]]=f1_sub_list
  } else {
 
    fits_full=foreach(ind3=1:length(uniq_miRNA), .packages=c("tidyverse", "reshape2", "lme4", "SummarizedExperiment")) %dopar% {
      #cat(uniq_miRNA[ind3], "\n")
      ## subset sequences that map to the miRNA
      idx = which(rowData(se)$miRNA == uniq_miRNA[ind3])
      Y_all_sub = t(assay(se)[idx, ])
      Y_seq_labels_sub = rowData(se)$uniqueSequence[idx]
      
      ## compute median CPM for each sequence
      cpm_sub = Y_all_sub / (total_counts/1e6)
      median_cpm_sub = apply(cpm_sub, 2, median)
      
      ## filter sequences with median CPM less than the min_med_lcpm threshold
      if(!is.null(min_med_lcpm)){
        idx2 = which(log(median_cpm_sub) > min_med_lcpm)
        Y_all_sub = Y_all_sub[ ,idx2]
        Y_seq_labels_sub = Y_seq_labels_sub[idx2]
      }
      
      ## format data to fit glmer.nb models
      data_wide = as.data.frame(as.matrix(Y_all_sub))
      colnames(data_wide) = Y_seq_labels_sub
      sample_labels = rownames(Y_all_sub)
      data_wide = cbind(col_group, total_counts, sample_labels, data_wide)
      data_long = melt(data_wide, 
                       id.vars = c("col_group", "total_counts", "sample_labels"),
                       variable.name = "sequence", 
                       value.name = "count")
      tryCatch({
        f1 = glmer.nb(count ~ col_group + offset(log(total_counts/1e6)) + 
                        (1+col_group|sequence) + (1|sample_labels), 
                      data=data_long, 
                      control=(glmerControl(optimizer="bobyqa", 
                                            tolPwrss = 1e-3, 
                                            optCtrl=list(maxfun=2e5))))
      }, error = function(e){cat("ERROR :", uniq_miRNA[ind3])})
      f1_list=list()
      f1_list[[uniq_miRNA[ind3]]]=f1
      return(f1_list)
    }
    fits_red=foreach(ind3=1:length(uniq_miRNA), .packages=c("tidyverse", "reshape2", "lme4", "SummarizedExperiment")) %dopar% {
      #cat(uniq_miRNA[ind3], "\n")
      ## subset sequences that map to the miRNA
      idx = which(rowData(se)$miRNA == uniq_miRNA[ind3])
      Y_all_sub = t(assay(se)[idx, ])
      Y_seq_labels_sub = rowData(se)$uniqueSequence[idx]
      
      ## compute median CPM for each sequence
      cpm_sub = Y_all_sub / (total_counts/1e6)
      median_cpm_sub = apply(cpm_sub, 2, median)
      
      ## filter sequences with median CPM less than the min_med_lcpm threshold
      if(!is.null(min_med_lcpm)){
        idx2 = which(log(median_cpm_sub) > min_med_lcpm)
        Y_all_sub = Y_all_sub[ ,idx2]
        Y_seq_labels_sub = Y_seq_labels_sub[idx2]
      }
      
      ## format data to fit glmer.nb models
      data_wide = as.data.frame(as.matrix(Y_all_sub))
      colnames(data_wide) = Y_seq_labels_sub
      sample_labels = rownames(Y_all_sub)
      data_wide = cbind(col_group, total_counts, sample_labels, data_wide)
      data_long = melt(data_wide, 
                       id.vars = c("col_group", "total_counts", "sample_labels"),
                       variable.name = "sequence", 
                       value.name = "count")
      tryCatch({
        ## a bunch of these produce the message:
        ## boundary (singular) fit: see help('isSingular')
        ## could we substitute f1 with f1_sub if we get this message?
        f1_sub = glmer.nb(count ~ col_group + offset(log(total_counts/1e6)) + 
                            (1|sequence) + (1|sample_labels), 
                          data=data_long, 
                          control=(glmerControl(optimizer="bobyqa", 
                                                tolPwrss = 1e-3, 
                                                optCtrl=list(maxfun=2e5))))
        
      }, error = function(e){cat("ERROR :", uniq_miRNA[ind3])})
      f1_list=list()
      f1_list[[uniq_miRNA[ind3]]]=f1_sub
      return(f1_list)
    }
  fits=list()
  fits[["miRglmm"]]=unlist(fits_full)
  fits[["miRglmm_reduced"]]=unlist(fits_red)
  }

  return(fits)
}

ncores=8 #match to what was requested on BH
if (ncores>1){
cl=makeCluster(ncores)
registerDoParallel(cl)
}
#startTime=Sys.time()
tst_par = miRglmm(sims[[1]]$sim_se, ncores)
if (ncores>1){
stopCluster(cl)
}
#endTime=Sys.time()
#endTime-startTime

uniq_miRNA = unique(rowData(sims[[1]]$sim_se)$miRNA)
true_logFC = rep(log(1), length(uniq_miRNA))
true_logFC[which(uniq_miRNA %in% as.character(sims[[1]]$change_miRNA_up$miRNA))] = log(2)
true_logFC[which(uniq_miRNA %in% as.character(sims[[1]]$change_miRNA_down$miRNA))] = log(0.5)
true_logFC_BvsA = true_logFC*-1
tst$true_logFC=true_logFC_BvsA

plot(x=true_logFC_BvsA, y=tst$col_groupB)
