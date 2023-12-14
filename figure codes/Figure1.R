library(ggplot2)

load("bladder_testes_data_subset_filtered2.rda")
miRNA_counts = apply(assay(bladder_testes_data_subset_filtered2), 2, function(x) by(x, rowData(bladder_testes_data_subset_filtered2)$miRNA, sum))
total_counts=colSums(assay(bladder_testes_data_subset_filtered2))

#find canonical sequences from ERCC file
ratio_pool=read.csv(file="ERCC files/FINAL_Ratiometric_SynthA_and_SynthB-1_fixlast3.csv", sep="\t")



########################################Figure 1A
miRNA_plot="hsa-let-7g-5p"

#find aggregated counts
agg_choose=data.frame(count=miRNA_counts[rownames(miRNA_counts)==miRNA_plot,])
agg_choose$total_counts=total_counts
agg_choose$cpm=agg_choose$count/(agg_choose$total_counts/1000000)
agg_choose$sequence="miRNA"
agg_choose$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue

# calculate sequence median CPM expression
idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA==miRNA_plot)
  Y_all_sub=t(assay(bladder_testes_data_subset_filtered2)[idx,])
  Y_seq_labels_sub=rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]
  cpm_sub=Y_all_sub/(total_counts/1000000)
  median_cpm_sub=apply(cpm_sub, 2, median)
  logmedCPM=data.frame(sequence=Y_seq_labels_sub, logmedCPM=log(median_cpm_sub))
  
  #find canonical sequence
  ratio_pool_in=subset(ratio_pool, ratio.seqID==miRNA_plot)
  data_canon=data.frame(cpm=cpm_sub[, which(Y_seq_labels_sub==ratio_pool_in$sequence)])
  data_canon$sequence=ratio_pool_in$sequence
  data_canon$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  
  

# from remaining sequences find the top two expressing sequences
  data_left=data.frame(as.matrix(cpm_sub[, which(Y_seq_labels_sub!=ratio_pool_in$sequence)]))
  colnames(data_left)=Y_seq_labels_sub[which(Y_seq_labels_sub!=ratio_pool_in$sequence)]
  logmedCPM_left=subset(logmedCPM, sequence %in% colnames(data_left))
  next_two_seq=logmedCPM_left$sequence[order(logmedCPM_left$logmedCPM, decreasing=TRUE)[1:2]]
  data_next_one=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[1])])
  data_next_one$sequence=next_two_seq[1]
  data_next_one$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  
  data_next_two=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[2])])
  data_next_two$sequence=next_two_seq[2]
  data_next_two$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  data_next_two=rbind(data_next_one, data_next_two)
  
  #boxplots
  data_comb=rbind(agg_choose[, c("cpm", "sequence", "col_group")], data_canon[,c("cpm", "sequence", "col_group")],
                  data_next_two[, c("cpm", "sequence", "col_group")])
  ggplot(data_comb, aes(x=col_group, y=cpm))+geom_boxplot()+facet_wrap(~factor(sequence, levels=unique(data_comb$sequence)), scales="free")+xlab("tissue")
  ggsave(paste0("figures/figure1A_", miRNA_plot, ".tif"), plot=last_plot(), device="tiff", width=354, height=190, units="mm", dpi=320, bg="white")
  
  

########################################Figure 1B
  miRNA_plot="hsa-miR-26a-5p"
  
  #find aggregated counts
  agg_choose=data.frame(count=miRNA_counts[rownames(miRNA_counts)==miRNA_plot,])
  agg_choose$total_counts=total_counts
  agg_choose$cpm=agg_choose$count/(agg_choose$total_counts/1000000)
  agg_choose$sequence="miRNA"
  agg_choose$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  
  # calculate sequence median CPM expression
  idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA==miRNA_plot)
  Y_all_sub=t(assay(bladder_testes_data_subset_filtered2)[idx,])
  Y_seq_labels_sub=rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]
  cpm_sub=Y_all_sub/(total_counts/1000000)
  median_cpm_sub=apply(cpm_sub, 2, median)
  logmedCPM=data.frame(sequence=Y_seq_labels_sub, logmedCPM=log(median_cpm_sub))
  
  #find canonical sequence
  ratio_pool_in=subset(ratio_pool, ratio.seqID==miRNA_plot)
  data_canon=data.frame(cpm=cpm_sub[, which(Y_seq_labels_sub==ratio_pool_in$sequence)])
  data_canon$sequence=ratio_pool_in$sequence
  data_canon$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  
  
  
  # from remaining sequences find the top two expressing sequences
  data_left=data.frame(as.matrix(cpm_sub[, which(Y_seq_labels_sub!=ratio_pool_in$sequence)]))
  colnames(data_left)=Y_seq_labels_sub[which(Y_seq_labels_sub!=ratio_pool_in$sequence)]
  logmedCPM_left=subset(logmedCPM, sequence %in% colnames(data_left))
  next_two_seq=logmedCPM_left$sequence[order(logmedCPM_left$logmedCPM, decreasing=TRUE)[1:2]]
  data_next_one=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[1])])
  data_next_one$sequence=next_two_seq[1]
  data_next_one$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  
  data_next_two=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[2])])
  data_next_two$sequence=next_two_seq[2]
  data_next_two$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  data_next_two=rbind(data_next_one, data_next_two)
  
  #boxplots
  data_comb=rbind(agg_choose[, c("cpm", "sequence", "col_group")], data_canon[,c("cpm", "sequence", "col_group")],
                  data_next_two[, c("cpm", "sequence", "col_group")])
  ggplot(data_comb, aes(x=col_group, y=cpm))+geom_boxplot()+facet_wrap(~factor(sequence, levels=unique(data_comb$sequence)), scales="free")+xlab("tissue")
  ggsave(paste0("figures/figure1B_", miRNA_plot, ".tif"), plot=last_plot(), device="tiff", width=354, height=190, units="mm", dpi=320, bg="white")
  
  
  ########################################Figure 1C
  miRNA_plot="hsa-let-7a-5p"
  
  #find aggregated counts
  agg_choose=data.frame(count=miRNA_counts[rownames(miRNA_counts)==miRNA_plot,])
  agg_choose$total_counts=total_counts
  agg_choose$cpm=agg_choose$count/(agg_choose$total_counts/1000000)
  agg_choose$sequence="miRNA"
  agg_choose$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  
  # calculate sequence median CPM expression
  idx=which(rowData(bladder_testes_data_subset_filtered2)$miRNA==miRNA_plot)
  Y_all_sub=t(assay(bladder_testes_data_subset_filtered2)[idx,])
  Y_seq_labels_sub=rowData(bladder_testes_data_subset_filtered2)$uniqueSequence[idx]
  cpm_sub=Y_all_sub/(total_counts/1000000)
  median_cpm_sub=apply(cpm_sub, 2, median)
  logmedCPM=data.frame(sequence=Y_seq_labels_sub, logmedCPM=log(median_cpm_sub))
  
  #find canonical sequence
  ratio_pool_in=subset(ratio_pool, ratio.seqID==miRNA_plot)
  data_canon=data.frame(cpm=cpm_sub[, which(Y_seq_labels_sub==ratio_pool_in$sequence)])
  data_canon$sequence=ratio_pool_in$sequence
  data_canon$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  
  
  
  # from remaining sequences find the top two expressing sequences
  data_left=data.frame(as.matrix(cpm_sub[, which(Y_seq_labels_sub!=ratio_pool_in$sequence)]))
  colnames(data_left)=Y_seq_labels_sub[which(Y_seq_labels_sub!=ratio_pool_in$sequence)]
  logmedCPM_left=subset(logmedCPM, sequence %in% colnames(data_left))
  next_two_seq=logmedCPM_left$sequence[order(logmedCPM_left$logmedCPM, decreasing=TRUE)[1:2]]
  data_next_one=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[1])])
  data_next_one$sequence=next_two_seq[1]
  data_next_one$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  
  data_next_two=data.frame(cpm=data_left[, which(colnames(data_left)==next_two_seq[2])])
  data_next_two$sequence=next_two_seq[2]
  data_next_two$col_group=colData(bladder_testes_data_subset_filtered2)$cell_tissue
  data_next_two=rbind(data_next_one, data_next_two)
  
  #boxplots
  data_comb=rbind(agg_choose[, c("cpm", "sequence", "col_group")], data_canon[,c("cpm", "sequence", "col_group")],
                  data_next_two[, c("cpm", "sequence", "col_group")])
  ggplot(data_comb, aes(x=col_group, y=cpm))+geom_boxplot()+facet_wrap(~factor(sequence, levels=unique(data_comb$sequence)), scales="free")+xlab("tissue")
  ggsave(paste0("figures/figure1C_", miRNA_plot, ".tif"), plot=last_plot(), device="tiff", width=354, height=190, units="mm", dpi=320, bg="white")
  
  
