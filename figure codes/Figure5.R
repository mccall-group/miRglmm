load(file='bladder_testes_results/filter_neg1_processed_results.rda')
LRTp=data.frame(results[["LRTp"]])
ggplot(LRTp, aes(x=LRTp))+geom_histogram(color="black", fill="gray", bins=50)+xlab('LRT p')+ylab('number of miRNA')
ggsave("figures/figure5.tif", plot=last_plot(), device="tiff", width=200, height=125, units="mm", dpi=320, bg="white")
