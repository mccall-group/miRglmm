## produce figure 2 plot
library(ggpubr)
## read in data produced by extract_monocyte_data_subset.R
load("bladder_testes_data_subset_filtered.rda")
load("bladder_testes_subset_total_counts.rda")

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




p1=ggplot(data.frame("log_total_counts"=log(total_counts)), aes(x=log_total_counts))+geom_histogram(color="black", fill="gray", bins=10)+xlab('log(total counts)')+ylab('number of subjects')+ggtitle("Total counts")+theme(plot.title=element_text(hjust=0.5))
p2=ggplot(data.frame("log_medCPM"=log(n_seq_out$median_cpm)), aes(x=log_medCPM))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of miRNA')+ggtitle("miRNA expression")+theme(plot.title=element_text(hjust=0.5))
cpm_mat=t(assay(bladder_testes_data_subset_filtered))/(total_counts/1e6)
med_cpm_mat=apply(cpm_mat, 2,median)
log_med_cpm_mat=log(med_cpm_mat+1)
scaleFUN=function(x) format(x, scientific=TRUE, digits=1)  #sprintf("%.1f", x)
p4=ggplot(data.frame(log_med_cpm_mat), aes(x=log_med_cpm_mat))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM+1)')+ylab('number of sequences')+ggtitle("sequence expression")+theme(plot.title=element_text(hjust=0.5))+scale_y_continuous(trans='log', labels=scaleFUN)
p3=ggplot(data.frame("n_seq"=n_seq_out$nseq_miRNA), aes(x=n_seq))+geom_histogram(color="black", fill="gray", bins=100)+xlab('number of sequences')+ylab('number of miRNA')+ggtitle("number of sequences per miRNA")+theme(plot.title=element_text(hjust=0.5))+scale_x_continuous(trans='log')
#ggarrange(p1, p2, p3, p4, ncol=2, nrow=4)

#############################################################

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
p6=ggplot(data.frame("log_medCPM"=log(n_seq_out$median_cpm)), aes(x=log_medCPM))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM+1)')+ylab('number of miRNA')+ggtitle("miRNA expression")+theme(plot.title=element_text(hjust=0.5))
cpm_mat=t(assay(bladder_testes_data_subset_filtered2))/(total_counts_new/1e6)
med_cpm_mat=apply(cpm_mat, 2,median)
log_med_cpm_mat=log(med_cpm_mat)
log_med_cpm_mat_plot=log(med_cpm_mat+1)
p8=ggplot(data.frame(log_med_cpm_mat_plot), aes(x=log_med_cpm_mat_plot))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM+1)')+ylab('number of sequences')+ggtitle("sequence expression")+theme(plot.title=element_text(hjust=0.5))+scale_y_continuous(trans='log', labels=scaleFUN)

p7=ggplot(data.frame("n_seq"=n_seq_out$nseq_miRNA), aes(x=n_seq))+geom_histogram(color="black", fill="gray", bins=100)+xlab('number of sequences')+ylab('number of miRNA')+ggtitle("number of sequences per miRNA")+theme(plot.title=element_text(hjust=0.5))+scale_x_continuous(trans='log')
#ggarrange(p1, p5, p2, p6, p3, p7, p4, p8, ncol=2, nrow=4)




##############################################################
#filter low expressing sequences
idx2 = which(log_med_cpm_mat > -1)
bladder_testes_data_subset_filtered3 <- bladder_testes_data_subset_filtered2[idx2, ]
log_med_cpm_mat_plot=log_med_cpm_mat_plot[idx2]
## calculate miRNA-level summaries after filtering low expressing sequences
nseq_miRNA = data.frame(table(rowData(bladder_testes_data_subset_filtered3)$miRNA))

#create empty plots for total counts and miRNA expression after filtering--these are not impacted by the filtering
p9=ggplot()+theme_void()
p10=ggplot()+theme_void()

p12=ggplot(data.frame(log_med_cpm_mat_plot), aes(x=log_med_cpm_mat_plot))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM+1)')+ylab('number of sequences')+ggtitle("sequence expression")+theme(plot.title=element_text(hjust=0.5))

p11=ggplot(nseq_miRNA, aes(x=Freq))+geom_histogram(color="black", fill="gray", bins=100)+xlab('number of sequences')+ylab('number of miRNA')+ggtitle("number of sequences per miRNA")+theme(plot.title=element_text(hjust=0.5))+scale_x_continuous(trans='log')

ggarrange(p1, p5, p9, p2, p6, p10, p3, p7, p11, p4, p8, p12, ncol=3, nrow=4)
ggsave("figures/figure2.tif", plot=last_plot(), device="tiff", width=500, height=250, units="mm", dpi=320, bg="white")
