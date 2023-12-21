library(vegan)
library(ggplot2)

load("ERCC/ERCC_filtered.rda")


Bray_curtis_NMDS=metaMDS(as.matrix(t(assay(panel_B_filter))), k=2)
out=Bray_curtis_NMDS$points
col_data=as.data.frame(colData(panel_B_filter))
col_data$NMDS1=out[,1]
col_data$NMDS2=out[,2]


ggplot(col_data, aes(x=NMDS1, y=NMDS2, color=Lab))+geom_point(aes(shape=Pool, size=4))+
  scale_size(guide='none')+guides(color=guide_legend(override.aes=list(size=3)))+
  guides(shape=guide_legend(override.aes=list(size=3)))+theme(legend.text=element_text(size=12))
ggsave("figures/figure7.tif", plot=last_plot(), device="tiff", width=250, height=150, units="mm", dpi=320, bg="white")

