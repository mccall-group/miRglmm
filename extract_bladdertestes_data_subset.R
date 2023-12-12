## extract the miRNA data for bladder/testes samples from the microRNAome

## read in sequence level data
## this file is 9.8GB; be patient
load(file="/scratch/mmccall2_lab/microRNAome/seqlevel_microRNAome.rda")

## identify the samples from study 46 that were not flagged as poor quality
ind_col = which(seqlevel_microRNAome$study_id==46 &
                  !seqlevel_microRNAome$flagged)

## identify sequences that align to miRNAs
idx_iso = which(rowData(seqlevel_microRNAome)$isomiR.miRNA != "")
idx_exact = which(rowData(seqlevel_microRNAome)$exact.miRNA != "")

## create a variable (called match) that denotes whether a sequence aligns to a miRNA or isomiR
rowData(seqlevel_microRNAome)$match = NA
rowData(seqlevel_microRNAome)$match[idx_iso] = "isomiR"
rowData(seqlevel_microRNAome)$match[idx_exact] = "exact miRNA"

## create a variable (called label) that denotes the miRNA a sequence aligns to
rowData(seqlevel_microRNAome)$miRNA_label = NA
rowData(seqlevel_microRNAome)$miRNA_label[idx_iso] = as.character(rowData(seqlevel_microRNAome)$isomiR.miRNA)[idx_iso]
rowData(seqlevel_microRNAome)$miRNA_label[idx_exact] = as.character(rowData(seqlevel_microRNAome)$exact.miRNA)[idx_exact]

## create dataset of just exact matches or isomiRs
ind_row = union(idx_exact, idx_iso)
bladder_testes_data_subset = seqlevel_microRNAome[ind_row, ind_col]

## convert unique sequence from factor to character 
## this make the object much smaller
rowData(bladder_testes_data_subset)$uniqueSequence = as.character(rowData(bladder_testes_data_subset)$uniqueSequence)

## remove original miRge variables from rowData
## this gets the object under 1GB in memory
rowData(bladder_testes_data_subset) = rowData(bladder_testes_data_subset)[ ,-c(1,3:12)]

## rename miRNA_label as miRNA
names(rowData(bladder_testes_data_subset))[3] = "miRNA"

## convert match variable to factor
rowData(bladder_testes_data_subset)$match = factor(rowData(bladder_testes_data_subset)$match,
                                             levels = c("exact miRNA", "isomiR"))

## save data object
save(bladder_testes_data_subset, file="bladder_testes_data_subset.rda")
