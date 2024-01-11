library(tidyverse)
library(reshape2)
library(lme4)
library(SummarizedExperiment)
library(DESeq2)
library(edgeR)

#load truth data
ratio_pool=read.csv(file="ERCC files/FINAL_Ratiometric_SynthA_and_SynthB-1_fixlast3.csv", sep="\t")
ratio_pool$A=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.A.))
ratio_pool$B=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.B.))
ratio_pool$ratio=ratio_pool$B/ratio_pool$A
ratio_pool$true_logFC=log(ratio_pool$ratio)
out=ratio_pool["true_logFC"]
rownames(out)=ratio_pool$ratio.seqID



#load processed results
load('ERCC/filter_agg_betahats.rda')

filters=c("no filter", "-1", "-0.5", "0", "0.5", "1", "1.5", "2")


#combine filtered agg with filtered seq models to get table of betas for each model
for (ind in seq(1, length(filters))){
  filter_in=filters[ind]
  beta_hat_in=beta_hat[[filter_in]]
  beta_hat_in=transform(merge(beta_hat_in, out, by='row.names'), row.names=Row.names, Row.names=NULL)
  diff_squared=(beta_hat_in[,!(colnames(beta_hat_in)=="true_logFC")]-beta_hat_in[,"true_logFC"])^2
  
  if (ind==1){
    MSE=data.frame(colMeans(diff_squared, na.rm=TRUE))
  } else {
    MSE=cbind(MSE, data.frame(colMeans(diff_squared, na.rm=TRUE)))
  }
  
}
colnames(MSE)=filters

MSE_long=stack(MSE)
colnames(MSE_long)=c("MSE", "filter")
MSE_long$method=rep(rownames(MSE), times=ncol(MSE))
ggplot(MSE_long, aes(x=filter, y=MSE, group=method, color=method))+
  geom_line()

ggsave("figures/figure11.tif", plot=last_plot(), device="tiff", width=200, height=150, units="mm", dpi=320, bg="white")
