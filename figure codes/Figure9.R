library(ggplot2)


#load truth data
ratio_pool=read.csv(file="ERCC files/FINAL_Ratiometric_SynthA_and_SynthB-1_fixlast3.csv", sep="\t")
ratio_pool$A=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.A.))
ratio_pool$B=as.numeric(gsub("x","",ratio_pool$X.SAMPLE.B.))
ratio_pool$ratio=ratio_pool$B/ratio_pool$A
ratio_pool$true_logFC=log(ratio_pool$ratio)
out=ratio_pool["true_logFC"]
rownames(out)=ratio_pool$ratio.seqID



#load processed results
load('ERCC/allfilters_processed_results.rda')


beta_hat=results[["beta_hat"]]
beta_hat=transform(merge(beta_hat, out, by='row.names'), row.names=Row.names, Row.names=NULL)
beta_hat=beta_hat[,c("filter..1", "no.filter","miRglmnb","DESeq2","edgeR","limmavoom", "true_logFC")]
colnames(beta_hat)=c("filter -1","no filter","miRglmnb","DESeq2","edgeR","limma-voom", "True LogFC")

diff=beta_hat[,!(colnames(beta_hat)=="True LogFC")]-beta_hat[,"True LogFC"]
diff_long=stack(diff)
colnames(diff_long)=c("diff", "method")
diff_long$`True LogFC`=as.factor(round(as.numeric(as.character(rep(beta_hat[,"True LogFC"], times=ncol(diff)))), digits=3))
ggplot(diff_long, aes(x=`True LogFC`, y=diff, fill=method))+geom_boxplot()+ylab('Pool LogFC Estimate-Truth')+xlab('True Pool LogFC')

ggsave("figures/figure9.tif", plot=last_plot(), device="tiff", width=300, height=150, units="mm", dpi=320, bg="white")
