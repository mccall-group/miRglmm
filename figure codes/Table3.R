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
load('ERCC/allfilters_processed_results_miRglmmreducedonly.rda')


beta_hat=results[["beta_hat"]]
beta_hat=transform(merge(beta_hat, out, by='row.names'), row.names=Row.names, Row.names=NULL)
beta_hat=beta_hat[,c("filter..1","filter..0.5","filter.0","filter.0.5","filter.1" ,"filter.1.5","filter.2",
            "no.filter","miRglmnb","DESeq2","edgeR","limmavoom", "true_logFC")]
colnames(beta_hat)=c("filter -1","filter -0.5","filter 0","filter 0.5","filter 1" ,"filter 1.5","filter 2",
                         "no filter","miRglmnb","DESeq2","edgeR","limma-voom", "True LogFC")


diff_squared=(beta_hat[,!(colnames(beta_hat)=="True LogFC")]-beta_hat[,"True LogFC"])^2
MSE=data.frame(MSE=colMeans(diff_squared, na.rm=TRUE))
n_miRNA=data.frame(n_miRNA=dim(beta_hat)[1]-colSums(is.na(beta_hat[,!(colnames(beta_hat)=="True LogFC")])))
n_seq=results[["n_seq"]]
n_seq_used=data.frame(n_seq_min=apply(n_seq,2,min, na.rm=TRUE), n_seq_max=apply(n_seq, 2, max, na.rm=TRUE))
rownames(n_seq_used)=c("filter -1","filter -0.5","filter 0","filter 0.5","filter 1" ,"filter 1.5","filter 2",
                       "no filter")


CI_LL=results[["CI_LL"]]
CI_LL=transform(merge(CI_LL, out, by='row.names'), row.names=Row.names, Row.names=NULL)
CI_LL=CI_LL[,c("filter..1","filter..0.5","filter.0","filter.0.5","filter.1" ,"filter.1.5","filter.2",
                     "no.filter","miRglmnb","DESeq2","limmavoom", "true_logFC")]
colnames(CI_LL)=c("filter -1","filter -0.5","filter 0","filter 0.5","filter 1" ,"filter 1.5","filter 2",
                     "no filter","miRglmnb","DESeq2","limma-voom", "True LogFC")
gr_LL=CI_LL[,!(colnames(CI_LL)=="True LogFC")]<=CI_LL[,"True LogFC"]


CI_UL=results[["CI_UL"]]
CI_UL=transform(merge(CI_UL, out, by='row.names'), row.names=Row.names, Row.names=NULL)
CI_UL=CI_UL[,c("filter..1","filter..0.5","filter.0","filter.0.5","filter.1" ,"filter.1.5","filter.2",
               "no.filter","miRglmnb","DESeq2","limmavoom", "true_logFC")]
colnames(CI_UL)=c("filter -1","filter -0.5","filter 0","filter 0.5","filter 1" ,"filter 1.5","filter 2",
                  "no filter","miRglmnb","DESeq2","limma-voom", "True LogFC")
lt_UL=CI_UL[,!(colnames(CI_UL)=="True LogFC")]>=CI_UL[,"True LogFC"]

cov_prob=data.frame("coverage probability"=colMeans(gr_LL& lt_UL, na.rm=TRUE))

table_out=cbind(MSE, "n_miRNA"=n_miRNA[match(rownames(MSE), rownames(n_miRNA)),], 
                n_seq_used[match(rownames(MSE), rownames(n_seq_used)),],
                "coverage probability"=cov_prob[match(rownames(MSE), rownames(cov_prob)),1])
write.csv(table_out, file="figures/Table3.csv")


