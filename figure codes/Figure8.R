####### Figure 8A- ERCC summary histograms
library(ggplot2)
library(DESeq2)
library(ggpubr)
library(SummarizedExperiment)
## read in data produced by extract_bladder_testes_data_subset.R
load("ERCC/ERCC_filtered.rda")
total_counts=colSums(assay(panel_B_filter))

## calculate miRNA-level summaries
nseq_miRNA = table(rowData(panel_B_filter)$miRNA)
miRNA_counts = apply(assay(panel_B_filter), 2, function(x) by(x, rowData(panel_B_filter)$miRNA, sum))
miRNA_cpms = t(t(miRNA_counts) / (total_counts / 1e6))
identical(names(nseq_miRNA), rownames(miRNA_cpms))
n_seq_out = data.frame(nseq_miRNA = nseq_miRNA, 
                       mean_cpm = rowMeans(miRNA_cpms), 
                       median_cpm = apply(miRNA_cpms, 1, median))
names(n_seq_out)[1] = "miRNA"
names(n_seq_out)[2] = "nseq_miRNA"




p1=ggplot(data.frame("log_total_counts"=log(total_counts)), aes(x=log_total_counts))+geom_histogram(color="black", fill="gray", bins=10)+xlab('log(total counts)')+ylab('number of subjects')+ggtitle("Total counts")+theme(plot.title=element_text(hjust=0.5))
p2=ggplot(data.frame("log_medCPM"=log(n_seq_out$median_cpm)), aes(x=log_medCPM))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of miRNA')+ggtitle("miRNA expression")+theme(plot.title=element_text(hjust=0.5))
cpm_mat=t(assay(panel_B_filter))/(total_counts/1e6)
med_cpm_mat=apply(cpm_mat, 2,median)
log_med_cpm_mat=log(med_cpm_mat)
#scaleFUN=function(x) format(x, scientific=TRUE, digits=1)  #sprintf("%.1f", x)
p4=ggplot(data.frame(log_med_cpm_mat), aes(x=log_med_cpm_mat))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of sequences')+ggtitle("sequence expression")+theme(plot.title=element_text(hjust=0.5))
p3=ggplot(data.frame("n_seq"=n_seq_out$nseq_miRNA), aes(x=n_seq))+geom_histogram(color="black", fill="gray", bins=100)+xlab('number of sequences')+ylab('number of miRNA')+ggtitle("number of sequences per miRNA")+theme(plot.title=element_text(hjust=0.5))
ggarrange(p1, p2, p3, p4, ncol=1, nrow=4)
ggsave("figures/figure8A.tif", plot=last_plot(), device="tiff", width=165, height=250, units="mm", dpi=320, bg="white")

#####Figure 8 Background vs miRNA expression

panelB=readRDS('/scratch/mmccall2_lab/ERCC_UMI/panel_B_SE.rds')
ddsColl=collapseReplicates(panelB, panelB$GEO_Sample)
raw_counts=assay(ddsColl)
total_counts=colSums(raw_counts)

cpm_sub=t(raw_counts)/(total_counts/1000000)
median_cpm_sub=apply(cpm_sub, 2, median)

source=rep("miRNA", length(median_cpm_sub))
idx=which(rowData(panelB)$miRNA=="-")
source[idx]="background (non-miRNA)"


df=data.frame(log(median_cpm_sub), source)
names(df)=c("log(median CPM)", "source")
ggplot(df, aes(x=`log(median CPM)`, fill=source))+geom_density(alpha=0.25)+ggtitle("Expression in sequences mapping to miRNA vs background expression")
ggsave("figures/figure8B.tif", plot=last_plot(), device="tiff", width=200, height=100, units="mm", dpi=320, bg="white")




