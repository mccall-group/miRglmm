## simulation function
sim_random_signal <- function(se_orig, bins, mean_effect=2, sd_effect=1){
  ## permute remove any column data effect
  new_sample_order = sample(ncol(se_orig))
  se = se_orig[ ,new_sample_order]
  sample_labels = colnames(se)

  ## randomly select miRNA stratified by bin to add signal
  set.seed(NULL) # does this need to be here?
  sample_change = bins %>% group_by(rank) %>% sample_n(size=2)
  change_miRNA_up = sample_change %>% group_by(rank) %>% sample_n(size=1)
  change_miRNA_down = sample_change %>%
    filter(!miRNA %in% change_miRNA_up$miRNA)

  ## add signal to chosen miRNA
  ## up is increased in samples 1:19
  iup = which(rowData(se)$miRNA %in% change_miRNA_up$miRNA)
  effect_up = rnorm(length(iup), mean_effect, sd_effect)
  effect_up[which(effect_up < 0)] = 0 # truncate effect at zero below
  assay(se)[iup, 1:19] = round(assay(se)[iup, 1:19] * effect_up)
  
  ## down is increased in samples 20:39
  idown = which(rowData(se)$miRNA %in% change_miRNA_down$miRNA)
  effect_down = rnorm(length(idown), mean_effect, sd_effect)
  effect_down[which(effect_down < 0)] = 0 # truncate effect at zero below
  assay(se)[idown, 20:39] = round(assay(se)[idown, 20:39] * effect_down)
  
  return(list(sim_se=se, 
              change_miRNA_up=change_miRNA_up, 
              change_miRNA_down=change_miRNA_down, 
              effect_up=effect_up, 
              effect_down=effect_down))
}

## generate N simulated data sets
N <- 100

## load filtered monocyte data subset
load(file = "monocyte_exact_subset_filtered2.rda")
load(file = "monocyte_n_seq_out.rda")

## set up bins to add signal to miRNA
## do not allow highest expressing miRNA to be chosen for added effect
library(tidyverse)
bins = n_seq_out %>% filter(median_cpm < 59000) %>% 
  mutate(rank = ntile(median_cpm, 20))
bins$miRNA <- as.character(bins$miRNA)

## run the function with mean 2 and sd 1
sims <- list()
for(k in 1:N) sims[[k]] <- sim_random_signal(exact_subset_filtered2, bins,
                                             mean_effect=2, sd_effect=1) 
save(sims, file = "sims_N100_m2_s1.rda")

## run the function with mean 2 and sd 0.5
sims <- list()
for(k in 1:N) sims[[k]] <- sim_random_signal(exact_subset_filtered2, bins,
                                             mean_effect=2, sd_effect=0.5) 
save(sims, file = "sims_N100_m2_s0.5.rda")

## run the function with mean 2 and sd 0.1
sims <- list()
for(k in 1:N) sims[[k]] <- sim_random_signal(exact_subset_filtered2, bins,
                                             mean_effect=2, sd_effect=0.1) 
save(sims, file = "sims_N100_m2_s0.1.rda")
