## filter miRNA outliers

## read in data produced by extract_study89_data_subset.R
load("study 89/study89_data_subset.rda")

## remove sequences containing N
library(stringr)
containsN = str_detect(rowData(study89_data_subset)$uniqueSequence, "N")
study89_data_subset = study89_data_subset[!containsN, ]

## compute total counts for each sample
total_counts = colSums(assays(study89_data_subset)$counts)

## keep exact and iso
#iexact <- which(rowData(study89_data_subset)$match == "exact miRNA")
#exact_subset = study89_data_subset[iexact, ]


save(total_counts, file = "study 89/study89_subset_total_counts.rda")

## only keep sequences with counts than zero in some sample
study89_data_subset_filtered = study89_data_subset[rowSums(assay(study89_data_subset)) > 0, ]
save(study89_data_subset_filtered, file = "study 89/study89_data_subset_filtered.rda")

## calculate miRNA-level summaries
nseq_miRNA = table(rowData(study89_data_subset_filtered)$miRNA)
miRNA_counts = apply(assay(study89_data_subset_filtered), 2, function(x) by(x, rowData(study89_data_subset_filtered)$miRNA, sum))
miRNA_cpms = t(t(miRNA_counts) / (total_counts / 1e6))
identical(names(nseq_miRNA), rownames(miRNA_cpms))
n_seq_out = data.frame(nseq_miRNA = nseq_miRNA, 
                       mean_cpm = rowMeans(miRNA_cpms), 
                       median_cpm = apply(miRNA_cpms, 1, median))
names(n_seq_out)[1] = "miRNA"
names(n_seq_out)[2] = "nseq_miRNA"


########################################################################

## keep only miRNA with a median CPM > 5
## this it the min of miRNA level median CPM in ERCC data
keep_miRNA = rownames(n_seq_out)[log(n_seq_out$median_cpm) > 5] 
study89_data_subset_filtered2 <- study89_data_subset_filtered[rowData(study89_data_subset_filtered)$miRNA %in% keep_miRNA, ]
save(study89_data_subset_filtered2, file = "study 89/study89_data_subset_filtered2.rda")
n_seq_out = n_seq_out[keep_miRNA, ]
save(n_seq_out, file = "study 89/study89_n_seq_out.rda")

