
library(SummarizedExperiment)
library(DESeq2)
# load ERCC object (there is a step of processing done before this---not sure where we want this code base to start from)
panelB=readRDS('/scratch/mmccall2_lab/ERCC_UMI/panel_B_SE.rds')
ddsColl=collapseReplicates(panelB, panelB$GEO_Sample)


#keep needed colData
coldata_in=data.frame(sample_id=colData(ddsColl)$Run, Exp=colData(ddsColl)$GEO_Sample, Pool=colData(ddsColl)$Pool, Lab=colData(ddsColl)$Lab, lib=colData(ddsColl)$lib.method.detail, lib_simple=colData(ddsColl)$lib.method.simple, instrument=colData(ddsColl)$Instrument.y)

#remove sequences that do not match to a miRNA
idx=which(rowData(panelB)$miRNA!="-")
panel_B_filter=SummarizedExperiment(assays=list(raw_counts=assay(ddsColl)[idx,]), rowData=rowData(panelB)[idx,], colData=coldata_in)

#remove miRNA with length>45
ratio_pool=read.csv(file="ERCC files/FINAL_Ratiometric_SynthA_and_SynthB-1_fixlast3.csv", sep="\t")
miRNA_keep=ratio_pool$ratio.seqID[ratio_pool$Length<45]
panel_B_filter=subset(panel_B_filter, rowData(panel_B_filter)$miRNA %in% miRNA_keep)


#remove sequences that match to multiple miRNA
uniq_seq=unique(rowData(panel_B_filter)$Template)
n_uniq_seq=length(uniq_seq)
Y_seq_labels_all=rowData(panel_B_filter)$Template
Y_miRNA_labels_all=rowData(panel_B_filter)$miRNA


dup_miRNA=c()
for (ind in seq(1, n_uniq_seq)){
  seq_in=uniq_seq[ind]
  idx=which(Y_seq_labels_all==seq_in)
  if (ind==1){
    n_miRNA=length(idx)
  } else {
    n_miRNA=rbind(n_miRNA, length(idx))
  }
  if (length(idx)>1){
    dup_miRNA=c(dup_miRNA, Y_miRNA_labels_all[idx])
  }
}

keep_seq=uniq_seq[which(n_miRNA==1)]
panel_B_filter=subset(panel_B_filter, rowData(panel_B_filter)$Template %in% keep_seq)


save(panel_B_filter, file='ERCC/ERCC_filtered.rda')

