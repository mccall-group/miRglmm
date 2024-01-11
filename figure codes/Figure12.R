library(ggplot2)
load(file='ERCC/allfilters_processed_results.rda')
LRTp=data.frame(LRTp=results[["LRTp"]][["filter..1"]])
ggplot(LRTp, aes(x=LRTp))+geom_histogram(color="black", fill="gray", bins=50)+xlab('LRT p')+ylab('number of miRNA')+xlim(c(0,NA))



ggsave("figures/figure12.tif", plot=last_plot(), device="tiff", width=200, height=125, units="mm", dpi=320, bg="white")
