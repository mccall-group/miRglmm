## filter miRNA outliers

## read in data produced by extract_monocyte_data_subset.R
load("monocyte_data_subset.rda")

## remove sequences containing N
library(stringr)
containsN = str_detect(rowData(monocyte_data_subset)$uniqueSequence, "N")
monocyte_data_subset = monocyte_data_subset[!containsN, ]

## compute total counts for each sample
total_counts = colSums(assays(monocyte_data_subset)$counts)

## select sequences that exactly match a miRNA
iexact <- which(rowData(monocyte_data_subset)$match == "exact miRNA")
exact_subset = monocyte_data_subset[iexact, ]

## why filter these two specific samples? 
## they have less than 1e6 total count but so do 4 other samples
## also this should probably be in extract_monocyte_data_subset.R
## where the other selection of samples is done
irm = which(exact_subset$sample_id %in% c("SRR5755891", "SRR5755870"))
exact_subset = exact_subset[ ,-irm]
total_counts = total_counts[-irm]
save(total_counts, file = "monocyte_subset_total_counts.rda")

## only keep sequences with counts than zero in some sample
exact_subset_filtered = exact_subset[rowSums(assay(exact_subset)) > 0, ]
save(exact_subset_filtered, file = "exact_subset_filtered.rda")

## calculate miRNA-level summaries
nseq_miRNA = table(rowData(exact_subset_filtered)$miRNA)
miRNA_counts = apply(assay(exact_subset_filtered), 2, function(x) by(x, rowData(exact_subset_filtered)$miRNA, sum))
miRNA_cpms = t(t(miRNA_counts) / (total_counts / 1e6))
identical(names(nseq_miRNA), rownames(miRNA_cpms))
n_seq_out = data.frame(nseq_miRNA = nseq_miRNA, 
                       mean_cpm = rowMeans(miRNA_cpms), 
                       median_cpm = apply(miRNA_cpms, 1, median))
names(n_seq_out)[1] = "miRNA"
names(n_seq_out)[2] = "nseq_miRNA"

########################################################################
## move these plots to files to remake the specific figures in the paper
hist(n_seq_out[,1])
hist(log(n_seq_out[,2]),120)
hist(log(n_seq_out[,3]), 100)

library(ggpubr)
p1=ggplot(data.frame("log_total_counts"=log(total_counts)), aes(x=log_total_counts))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(total counts)')+ylab('number of subjects')+ggtitle("Total counts")+theme(plot.title=element_text(hjust=0.5))
p2=ggplot(data.frame("log_medCPM"=log(n_seq_out[,3])), aes(x=log_medCPM))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of miRNA')+ggtitle("miRNA expression")+theme(plot.title=element_text(hjust=0.5))
cpm_mat=Y_all/(total_counts/1000000)
med_cpm_mat=apply(cpm_mat, 2,median)
log_med_cpm_mat=log(med_cpm_mat)
p4=ggplot(data.frame(log_med_cpm_mat), aes(x=log_med_cpm_mat))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of sequences')+ggtitle("sequence expression")+theme(plot.title=element_text(hjust=0.5))

p3=ggplot(data.frame("n_seq"=n_seq_out[,1]), aes(x=n_seq))+geom_histogram(color="black", fill="gray", bins=100)+xlab('number of sequences')+ylab('number of miRNA')+ggtitle("number of sequences per miRNA")+theme(plot.title=element_text(hjust=0.5))
ggarrange(p1, p2, p3, p4, ncol=1, nrow=4)
########################################################################

## keep only miRNA with a median CPM > 5
## this it the min of miRNA level median CPM in ERCC data
keep_miRNA = rownames(n_seq_out)[log(n_seq_out[,3]) > 5] 
exact_subset_filtered2 <- exact_subset_filtered[rowData(exact_subset_filtered)$miRNA %in% keep_miRNA, ]
save(exact_subset_filtered2, file = "exact_subset_filtered2.rda")
n_seq_out = n_seq_out[keep_miRNA, ]
save(n_seq_out, file = "n_seq_out.rda")

########################################################################
## move these plots to 
#make plots after filtering
for (ind2 in seq(1, length(uniq_miRNA))){
  uniq_miRNA_in=uniq_miRNA[ind2]
  idx=which(Y_miRNA_labels_all==uniq_miRNA_in)
  if (length(idx)>1){
    miRNA_sub=rowSums(Y_all[,idx])
  } else {
    miRNA_sub=Y_all[,idx]
  }
  cpm_miRNA=miRNA_sub/(total_counts/1000000)
  
  if (ind2==1){
    n_seq_out=c(length(idx), mean(cpm_miRNA), median(cpm_miRNA))
  } else {
    n_seq_out=rbind(n_seq_out, c(length(idx), mean(cpm_miRNA), median(cpm_miRNA)))
  }
}

total_counts_new=rowSums(Y_all)
library(ggpubr)
p1=ggplot(data.frame("log_total_counts"=log(total_counts_new)), aes(x=log_total_counts))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(total counts)')+ylab('number of subjects')+ggtitle("Total counts")+theme(plot.title=element_text(hjust=0.5))
p2=ggplot(data.frame("log_medCPM"=log(n_seq_out[,3])), aes(x=log_medCPM))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of miRNA')+ggtitle("miRNA expression")+theme(plot.title=element_text(hjust=0.5))
cpm_mat=Y_all/(total_counts/1000000)
med_cpm_mat=apply(cpm_mat, 2,median)
log_med_cpm_mat=log(med_cpm_mat)
p4=ggplot(data.frame(log_med_cpm_mat), aes(x=log_med_cpm_mat))+geom_histogram(color="black", fill="gray", bins=100)+xlab('log(median CPM)')+ylab('number of sequences')+ggtitle("sequence expression")+theme(plot.title=element_text(hjust=0.5))

p3=ggplot(data.frame("n_seq"=n_seq_out[,1]), aes(x=n_seq))+geom_histogram(color="black", fill="gray", bins=100)+xlab('number of sequences')+ylab('number of miRNA')+ggtitle("number of sequences per miRNA")+theme(plot.title=element_text(hjust=0.5))
ggarrange(p1, p2, p3, p4, ncol=1, nrow=4)




