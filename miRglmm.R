library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
miRglmm <- function(se,  col_group = c(rep("A", 19), rep("B",20)), 
                    min_med_lcpm = -1){
  
  ## for each miRNA (this could be parallelized)
  uniq_miRNA = unique(rowData(se)$miRNA)
  total_counts=colSums(assay(se))
  for(ind3 in seq(1, 10)){#} length(uniq_miRNA))){
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
      out_anova = anova(f1, f1_sub)
      out_CI = confint(f1, method="Wald")
      out_f1_coef = as.data.frame(t(fixef(f1)))
      out_f1_SE = cbind(sqrt(vcov(f1)[1,1]), sqrt(vcov(f1)[2,2]))
      out_f1_SE = cbind(out_f1_SE, t(out_CI[which(rownames(out_CI)=="col_groupB"), ]))
      colnames(out_f1_SE) = c("Int SE", "PoolB SE", "PoolB LL", "PoolB UL")
      out = VarCorr(f1)
      out_f1_varcomp = as.data.frame(cbind(out$sample_labels[1,1], 
                                           out$sequence[1,1], 
                                           out$sequence[2,2]))
      colnames(out_f1_varcomp) = c("random intercept var sample",
                                   "random intercept var sequence", 
                                   "random slope var sequence")
      out_f1_coef = cbind(out_f1_coef, out_f1_SE, out_f1_varcomp)
      out_f1_coef$miRNA = uniq_miRNA[ind3]
      #out_f1_coef$true_logFC=true_logFC_BvsA[ind3]
      #out_f1_coef$theta=f1$phis[1]
      out_f1_coef$n_seq = dim(Y_all_sub)[2]
      out_f1_coef$rand_slope_LRTp = out_anova$`Pr(>Chisq)`[2]
      if(ind3 == 1){
        f1_results = out_f1_coef
      } else{
        f1_results = rbind(f1_results, out_f1_coef)
      }
    }, error = function(e){cat("ERROR :", uniq_miRNA[ind3])})
  }
  return(f1_results)
}

startTime=Sys.time()
tst = miRglmm(sims[[1]]$sim_se)
endTime=Sys.time()
endTime-startTime


uniq_miRNA = unique(rowData(sims[[1]]$sim_se)$miRNA)
true_logFC = rep(log(1), length(uniq_miRNA))
true_logFC[which(uniq_miRNA %in% as.character(sims[[1]]$change_miRNA_up$miRNA))] = log(2)
true_logFC[which(uniq_miRNA %in% as.character(sims[[1]]$change_miRNA_down$miRNA))] = log(0.5)
true_logFC_BvsA = true_logFC*-1
tst$true_logFC=true_logFC_BvsA

plot(x=true_logFC_BvsA, y=tst$col_groupB)
