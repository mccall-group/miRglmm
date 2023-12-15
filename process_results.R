get_betas <- function(model_fits, var="col_group"){
  library(stringr)

  #miRglmm betas
  singular_warn=sapply(model_fits[["miRglmm"]], 'isSingular')
  #find full model betas
  all_coeff=sapply(model_fits[["miRglmm"]], "fixef")
  idx=which(str_detect(rownames(all_coeff), var))
  coeff_full=data.frame("full"=all_coeff[idx,], "singular_warning"=singular_warn)
  rownames(coeff_full)=colnames(all_coeff)
  
  #find reduced model betas
  all_coeff=sapply(model_fits[["miRglmm_reduced"]], "fixef")
  idx=which(str_detect(rownames(all_coeff), var))
  coeff_red=data.frame("reduced"=all_coeff[idx,])
  rownames(coeff_red)=colnames(all_coeff)
  glmm_betas=transform(merge(coeff_full, coeff_red, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  #if singularity warning choose reduced model betas
  betas=data.frame('miRglmm'=glmm_betas$full)
  rownames(betas)=rownames(glmm_betas)
  betas$miRglmm[which(glmm_betas$singular_warning==TRUE)]=glmm_betas$reduced[which(glmm_betas$singular_warning==TRUE)]
  
  
  #miRglmnb betas
  all_coeff=sapply(model_fits[["miRglmnb"]], '[[', "coefficients")
  idx=which(str_detect(rownames(all_coeff), var))
  out=data.frame('miRglmnb'=all_coeff[idx,])
  betas=transform(merge(betas, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  
  #DESeq2 betas
  library(DESeq2)
  res=results(model_fits[["DESeq2"]])
  out=data.frame('DESeq2'=log(2^res$log2FoldChange))
  rownames(out)=rownames(res)
  betas=transform(merge(betas, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  #edgeR betas
  out=data.frame('edgeR'=log(2^fits[["edgeR"]]$table$logFC))
  rownames(out)=rownames(fits[["edgeR"]]$table)
  betas=transform(merge(betas, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  #limma-voom betas
  all_coeff=coef(model_fits[["limvoom"]])
  idx=which(str_detect(colnames(all_coeff), var))
  out=data.frame('limmavoom'=log(2^all_coeff[, idx]))
  rownames(out)=rownames(all_coeff)
  betas=transform(merge(betas, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  return(betas)
  
}

get_SEs <- function(model_fits, var="col_group"){
  
  library(stringr)
  
  #miRglmm SEs
  singular_warn=sapply(model_fits[["miRglmm"]], 'isSingular')
  #find full model SEs
  all_SE=sapply(model_fits[["miRglmm"]], "vcov")
  idx1=which(str_detect(rownames(all_SE[[1]]), var)==TRUE)
  idx2=which(str_detect(colnames(all_SE[[1]]), var)==TRUE)
  SE_vec=data.frame('SE_full'=sapply(all_SE, function(x) sqrt(x[idx1,idx2])), "singular_warning"=singular_warn)
  rownames(SE_vec)=names(all_SE)
  
  
  #find reduced model SEs
  all_SE=sapply(model_fits[["miRglmm_reduced"]], "vcov")
  idx1=which(str_detect(rownames(all_SE[[1]]), var)==TRUE)
  idx2=which(str_detect(colnames(all_SE[[1]]), var)==TRUE)
  SE_vec_red=data.frame('SE_red'=sapply(all_SE, function(x) sqrt(x[idx1,idx2])))
  rownames(SE_vec_red)=names(all_SE)
  glmm_SEs=transform(merge(SE_vec, SE_vec_red, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  #if singularity warning choose reduced model betas
  SEs=data.frame('miRglmm'=glmm_SEs$SE_full)
  rownames(SEs)=rownames(glmm_SEs)
  SEs$miRglmm[which(glmm_SEs$singular_warning==TRUE)]=glmm_SEs$SE_red[which(glmm_SEs$singular_warning==TRUE)]
  
  
  #miRglmnb betas
  all_coeff=unlist(sapply(model_fits[["miRglmnb"]], function(x) summary(x)$coefficients[which(str_detect(rownames(summary(x)$coefficients), var)==TRUE), which(str_detect(colnames(summary(x)$coefficients), "Error")==TRUE)]))
  out=data.frame('miRglmnb'=all_coeff)
  rownames(out)=names(all_coeff)
  SEs=transform(merge(SEs, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  
  #DESeq2 betas
  library(DESeq2)
  res=results(model_fits[["DESeq2"]])
  #out=data.frame('beta_log2'=res$log2FoldChange, 'SE_log2'=res$lfcSE)
  #rownames(out)=rownames(res)
  #out$LL_log2=out$beta_log2-(1.96*out$SE_log2)
  #out$LL=log(2^out$LL_log2)
  #out$beta_log=log(2^out$beta_log2)
  #out$SE_log=(out$LL-out$beta_log)/-1.96 #confirmed that this provides same results as log(2^SE_log2)
  out=data.frame('DESeq2'=log(2^res$lfcSE))
  rownames(out)=rownames(res)
  SEs=transform(merge(SEs, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  
  #limma-voom betas
  SE_out=sqrt(model_fits[["limvoom"]]$s2.post)*model_fits[["limvoom"]]$stdev.unscaled[,which(str_detect(colnames(model_fits[["limvoom"]]$stdev.unscaled), var)==TRUE)]
  out=data.frame("limmavoom"=log(2^SE_out))
  rownames(out)=names(SE_out)
  SEs=transform(merge(SEs, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  return(SEs)
}


run_LRT <- function(full, reduced){
  uniq_miRNA=names(full)
out=data.frame("LRTp"=sapply(uniq_miRNA, function(row) anova(full[[row]], reduced[[row]])$`Pr(>Chisq)`[2]))
rownames(out)=uniq_miRNA
  return(out)
  
}

get_varcomp <- function(full, reduced){
  
  singular_warn=data.frame("singular"=sapply(full, 'isSingular'))
  #find full model varcomp
  uniq_miRNA=names(full)
  var_comp=data.frame('random_int_sample_var'=sapply(uniq_miRNA, function(row) VarCorr(full[[row]])$sample_labels[1,1]))
  var_comp2=data.frame('random_int_seq_var'=sapply(uniq_miRNA, function(row) VarCorr(full[[row]])$sequence[1,1]))
  var_comp=transform(merge(var_comp, var_comp2, by='row.names'), row.names=Row.names, Row.names=NULL)
  var_comp2=data.frame('random_slope_seq_var'=sapply(uniq_miRNA, function(row) VarCorr(full[[row]])$sequence[2,2]))
  var_comp=transform(merge(var_comp, var_comp2, by='row.names'), row.names=Row.names, Row.names=NULL)
  
  #find reduced model varcomp
  var_comp2=data.frame('random_int_sample_var_red'=sapply(uniq_miRNA, function(row) VarCorr(reduced[[row]])$sample_labels[1]))
  var_comp=transform(merge(var_comp, var_comp2, by='row.names'), row.names=Row.names, Row.names=NULL)
  var_comp2=data.frame('random_int_seq_var_red'=sapply(uniq_miRNA, function(row) VarCorr(reduced[[row]])$sequence[1]))
  var_comp=transform(merge(var_comp, var_comp2, by='row.names'), row.names=Row.names, Row.names=NULL)
  var_comp$`random_slope_seq_var_red`=NA
  var_comp=transform(merge(var_comp, singular_warn, by='row.names'), row.names=Row.names, Row.names=NULL)

  
  return(var_comp)

}


getCoverageInd <- function(betas, SEs, nominal_level = 0.95){
  z_alpha=qnorm(1-(1-nominal_level)/2)
  LL_mat=betas[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]-z_alpha*SEs[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]
  UL_mat=betas[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]+z_alpha*SEs[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]
  coverage_indicator=data.frame(1*((LL_mat<=betas$true_beta) & (UL_mat >= betas$true_beta)))
  return(coverage_indicator)
}

getCI_widths <- function(betas, SEs, nominal_level = 0.95){
  z_alpha=qnorm(1-(1-nominal_level)/2)
  LL_mat=betas[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]-z_alpha*SEs[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]
  UL_mat=betas[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]+z_alpha*SEs[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]
  CI_widths=UL_mat-LL_mat
  return(CI_widths)
}

getCI_LL <- function(betas, SEs, nominal_level = 0.95){
  z_alpha=qnorm(1-(1-nominal_level)/2)
  LL_mat=betas[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]-z_alpha*SEs[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]
  return(LL_mat)
}

getCI_UL <- function(betas, SEs, nominal_level = 0.95){
  z_alpha=qnorm(1-(1-nominal_level)/2)
  UL_mat=betas[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]+z_alpha*SEs[, c("miRglmm","miRglmnb","DESeq2","limmavoom")]
  return(UL_mat)
}

# #miRglmm coverage based on 95% Wald CI
# singular_warn=data.frame("singular"=sapply(model_fits[["miRglmm"]], 'isSingular'))
# #find full model CI limits
# uniq_miRNA=names(model_fits[["miRglmm"]])
# CI_lim=data.frame(t(data.frame(sapply(uniq_miRNA, function(row) confint(model_fits[["miRglmm"]][[row]], method="profile")[which(str_detect(rownames(confint(model_fits[["miRglmm"]][[row]], method="Wald")),var)==TRUE), ]))))
# rownames(CI_lim)=uniq_miRNA
# colnames(CI_lim)=c("LL", "UL")
# CI_lim=transform(merge(CI_lim, true_logFC_BvsA, by='row.names'), row.names=Row.names, Row.names=NULL)
# coverage=data.frame("miRglmm_full"=as.numeric(CI_lim$true_beta>=CI_lim$LL &  CI_lim$true_beta<=CI_lim$UL))
# rownames(coverage)=rownames(CI_lim)
# coverage=transform(merge(coverage, singular_warn, by='row.names'), row.names=Row.names, Row.names=NULL)
# 
# 
# #find reduced model betas
# uniq_miRNA=names(model_fits[["miRglmm_reduced"]])
# CI_lim=data.frame(t(data.frame(sapply(uniq_miRNA, function(row) confint(model_fits[["miRglmm_reduced"]][[row]], method="profile")[which(str_detect(rownames(confint(model_fits[["miRglmm_reduced"]][[row]], method="Wald")),var)==TRUE), ]))))
# rownames(CI_lim)=uniq_miRNA
# colnames(CI_lim)=c("LL", "UL")
# CI_lim=transform(merge(CI_lim, true_logFC_BvsA, by='row.names'), row.names=Row.names, Row.names=NULL)
# coverage_red=data.frame("miRglmm_reduced"=as.numeric(CI_lim$true_beta>=CI_lim$LL &  CI_lim$true_beta<=CI_lim$UL))
# rownames(coverage_red)=rownames(CI_lim)
# coverage=transform(merge(coverage, coverage_red, by='row.names'), row.names=Row.names, Row.names=NULL)
# 
# #if singularity warning choose reduced model betas
# coverage_out=data.frame('miRglmm'=coverage$miRglmm_full)
# rownames(coverage_out)=rownames(coverage)
# coverage_out$miRglmm[which(coverage$singular==TRUE)]=coverage$miRglmm_reduced[which(coverage$singular==TRUE)]
# coverage=coverage_out
# 
# #miRglmnb
# uniq_miRNA=names(model_fits[["miRglmnb"]])
# CI_lim=data.frame(t(data.frame(sapply(uniq_miRNA, function(row) suppressMessages(confint(model_fits[["miRglmnb"]][[row]])[which(str_detect(rownames(confint(model_fits[["miRglmnb"]][[row]], method="Wald")),var)==TRUE), ])))))
# rownames(CI_lim)=uniq_miRNA
# colnames(CI_lim)=c("LL", "UL")
# CI_lim=transform(merge(CI_lim, true_logFC_BvsA, by='row.names'), row.names=Row.names, Row.names=NULL)
# coverage_red=data.frame("miRglmm_reduced"=as.numeric(CI_lim$true_beta>=CI_lim$LL &  CI_lim$true_beta<=CI_lim$UL))