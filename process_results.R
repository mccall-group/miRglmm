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