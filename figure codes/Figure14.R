## produce figure 14 plot

####### Figure 14A- ERCC summary histograms
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
ggsave("figures/figure14A.tif", plot=last_plot(), device="tiff", width=165, height=250, units="mm", dpi=320, bg="white")




################################################# Panel B
## read in data produced by extract_bladder_testes_data_subset.R and keep exact sequence matches only
load("bladder_testes_data_subset_filtered.rda")
bladder_testes_data_subset_filtered=subset(bladder_testes_data_subset_filtered, rowData(bladder_testes_data_subset_filtered)$match=="exact miRNA")
#load("bladder_testes_subset_total_counts.rda")
total_counts=colSums(assay(bladder_testes_data_subset_filtered))

## calculate miRNA-level summaries
nseq_miRNA = table(rowData(bladder_testes_data_subset_filtered)$miRNA)
miRNA_counts = apply(assay(bladder_testes_data_subset_filtered), 2, function(x) by(x, rowData(bladder_testes_data_subset_filtered)$miRNA, sum))
miRNA_cpms = t(t(miRNA_counts) / (total_counts / 1e6))
identical(names(nseq_miRNA), rownames(miRNA_cpms))
n_seq_out = data.frame(nseq_miRNA = nseq_miRNA, 
                       mean_cpm = rowMeans(miRNA_cpms), 
                       median_cpm = apply(miRNA_cpms, 1, median))
names(n_seq_out)[1] = "miRNA"
names(n_seq_out)[2] = "nseq_miRNA"


## keep only miRNA with a median CPM > 5
## this it the min of miRNA level median CPM in ERCC data
keep_miRNA = rownames(n_seq_out)[log(n_seq_out$median_cpm) > 5] 
bladder_testes_data_subset_filtered2 <- bladder_testes_data_subset_filtered[rowData(bladder_testes_data_subset_filtered)$miRNA %in% keep_miRNA, ]

total_counts_new=colSums(assay(bladder_testes_data_subset_filtered2))

## calculate miRNA-level summaries after filtering low expressing miRNA
nseq_miRNA = table(rowData(bladder_testes_data_subset_filtered2)$miRNA)
miRNA_counts = apply(assay(bladder_testes_data_subset_filtered2), 2, function(x) by(x, rowData(bladder_testes_data_subset_filtered2)$miRNA, sum))
miRNA_cpms = t(t(miRNA_counts) / (total_counts_new / 1e6))
identical(names(nseq_miRNA), rownames(miRNA_cpms))
n_seq_out = data.frame(nseq_miRNA = nseq_miRNA, 
                       mean_cpm = rowMeans(miRNA_cpms), 
                       median_cpm = apply(miRNA_cpms, 1, median))
names(n_seq_out)[1] = "miRNA"
names(n_seq_out)[2] = "nseq_miRNA"




p5=ggplot(data.frame("log_total_counts"=log(total_counts_new)), aes(x=log_total_counts))+geom_histogram(color="black", fill="gray", bins=10)+xlab('log(total counts)')+ylab('number of subjects')+ggtitle("Total counts")+theme(plot.title=element_text(hjust=0.5))
p6=ggplot(data.frame("log_medCPM"=log(n_seq_out$median_cpm)), aes(x=log_medCPM))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of miRNA')+ggtitle("miRNA expression")+theme(plot.title=element_text(hjust=0.5))
cpm_mat=t(assay(bladder_testes_data_subset_filtered2))/(total_counts_new/1e6)
med_cpm_mat=apply(cpm_mat, 2,median)
log_med_cpm_mat=log(med_cpm_mat)
p8=ggplot(data.frame(log_med_cpm_mat), aes(x=log_med_cpm_mat))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of sequences')+ggtitle("sequence expression")+theme(plot.title=element_text(hjust=0.5))

p7=ggplot(data.frame("n_seq"=n_seq_out$nseq_miRNA), aes(x=n_seq))+geom_histogram(color="black", fill="gray", bins=100)+xlab('number of sequences')+ylab('number of miRNA')+ggtitle("number of sequences per miRNA")+theme(plot.title=element_text(hjust=0.5))
ggarrange(p5, p6, p7, p8, ncol=1, nrow=4)
ggsave("figures/figure14B.tif", plot=last_plot(), device="tiff", width=165, height=250, units="mm", dpi=320, bg="white")


####################################### monocyte exact match only log(median CPM)>5 miRNA filtered Fig14C
load(file = "monocyte_exact_subset_filtered2.rda")
total_counts=colSums(assay(exact_subset_filtered2))
## calculate miRNA-level summaries
nseq_miRNA = table(rowData(exact_subset_filtered2)$miRNA)
miRNA_counts = apply(assay(exact_subset_filtered2), 2, function(x) by(x, rowData(exact_subset_filtered2)$miRNA, sum))
miRNA_cpms = t(t(miRNA_counts) / (total_counts / 1e6))
identical(names(nseq_miRNA), rownames(miRNA_cpms))
n_seq_out = data.frame(nseq_miRNA = nseq_miRNA, 
                       mean_cpm = rowMeans(miRNA_cpms), 
                       median_cpm = apply(miRNA_cpms, 1, median))
names(n_seq_out)[1] = "miRNA"
names(n_seq_out)[2] = "nseq_miRNA"


p1=ggplot(data.frame("log_total_counts"=log(total_counts)), aes(x=log_total_counts))+geom_histogram(color="black", fill="gray", bins=10)+xlab('log(total counts)')+ylab('number of subjects')+ggtitle("Total counts")+theme(plot.title=element_text(hjust=0.5))
p2=ggplot(data.frame("log_medCPM"=log(n_seq_out[,3])), aes(x=log_medCPM))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of miRNA')+ggtitle("miRNA expression")+theme(plot.title=element_text(hjust=0.5))
cpm_mat=t(assay(exact_subset_filtered2))/(total_counts/1000000)
med_cpm_mat=apply(cpm_mat, 2,median)
log_med_cpm_mat=log(med_cpm_mat)
p4=ggplot(data.frame(log_med_cpm_mat), aes(x=log_med_cpm_mat))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of sequences')+ggtitle("sequence expression")+theme(plot.title=element_text(hjust=0.5))

p3=ggplot(data.frame("n_seq"=n_seq_out$nseq_miRNA), aes(x=n_seq))+geom_histogram(color="black", fill="gray", bins=100)+xlab('number of sequences')+ylab('number of miRNA')+ggtitle("number of sequences per miRNA")+theme(plot.title=element_text(hjust=0.5))
ggarrange(p1, p2, p3, p4, ncol=1, nrow=4)
ggsave("figures/figure14C.tif", plot=last_plot(), device="tiff", width=165, height=250, units="mm", dpi=320, bg="white")
