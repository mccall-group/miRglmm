library(UpSetR)

load(file='bladder_testes_results/filter_neg1_processed_results.rda')


sig_df=data.frame(sign(data.frame(results[["CI_LL"]]))==sign(data.frame(results[["CI_UL"]])))
colnames(sig_df)=colnames(data.frame(results[["CI_LL"]]))
rownames(sig_df)=rownames(data.frame(results[["CI_LL"]]))

#find edgeR significance from fit object
load(file='bladder_testes_results/filter_neg1_results.rda')
sig_df_edgeR=data.frame('edgeR'=fits[["edgeR"]][["table"]][["PValue"]]<0.05)
rownames(sig_df_edgeR)=rownames(fits[["edgeR"]][["table"]])
sig_df=transform(merge(sig_df, sig_df_edgeR, by='row.names'), row.names=Row.names, Row.names=NULL)




list_Input=list(miRglmm=rownames(sig_df)[sig_df$miRglmm==TRUE], DESeq2=rownames(sig_df)[sig_df$DESeq2==TRUE],`NB GLM`=rownames(sig_df)[sig_df$miRglmnb==TRUE], 
                edgeR=rownames(sig_df)[sig_df$edgeR==TRUE], `limma-voom`=rownames(sig_df)[sig_df$limmavoom==TRUE])
myplot=upset(fromList(list_Input), order.by="freq")
tiff(file="figures/figure5.tif", units="mm", height=100, width=200, res=320)
print(myplot)
dev.off()
