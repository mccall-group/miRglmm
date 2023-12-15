library(irr)
library(vcd)

load(file='bladder_testes_results/filter_neg1_processed_results.rda')

beta_df=data.frame(results[["beta_hat"]])
sig_df=data.frame(sign(data.frame(results[["CI_LL"]]))==sign(data.frame(results[["CI_UL"]])))
colnames(sig_df)=colnames(data.frame(results[["CI_LL"]]))
rownames(sig_df)=rownames(data.frame(results[["CI_LL"]]))

#find edgeR significance from fit object
load(file='bladder_testes_results/filter_neg1_results.rda')
sig_df_edgeR=data.frame('edgeR'=fits[["edgeR"]][["table"]][["PValue"]]<0.05)
rownames(sig_df_edgeR)=rownames(fits[["edgeR"]][["table"]])
sig_df=transform(merge(sig_df, sig_df_edgeR, by='row.names'), row.names=Row.names, Row.names=NULL)

#create interaction significant groups of all aggregated methods vs miRglmm
sig_groups_v_miRglmm=data.frame("miRglmnb"=as.character(interaction(sig_df$miRglmm, sig_df$miRglmnb)), "DESeq2"=as.character(interaction(sig_df$miRglmm, sig_df$DESeq2)),
                      "edgeR"=as.character(interaction(sig_df$miRglmm, sig_df$miRglmnb)), "limmavoom"=as.character(interaction(sig_df$miRglmm, sig_df$limmavoom)))
rownames(sig_groups_v_miRglmm)=rownames(sig_df)

sig_groups_v_miRglmm[sig_groups_v_miRglmm=="TRUE.TRUE"]="significant for both"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="TRUE.FALSE"]="significant for miRglmm only"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="FALSE.TRUE"]="significant for aggregated model only"
sig_groups_v_miRglmm[sig_groups_v_miRglmm=="FALSE.FALSE"]="significant for neither"

#combine with the betas to make sure miRNAs align with correct beta
comb_df=transform(merge(beta_df, sig_groups_v_miRglmm, by='row.names'), row.names=Row.names, Row.names=NULL)


#calculate ICC and Kappa to include in plotting

icc_out=data.frame(icc=t(sapply(seq(2,5), function(x) icc(comb_df[, c(1,x)], model="oneway", type="agreement", unit="single")$value)))
colnames(icc_out)=colnames(comb_df)[2:5]
kappa_out=data.frame(kappa=t(sapply(seq(2,5), function(x) Kappa(table(sig_df[,1], sig_df[,x]))[["Unweighted"]][["value"]])))
colnames(kappa_out)=colnames(sig_df)[2:5]

#scatter plots
ggplot(comb_df, aes(x=miRglmnb.x, y=miRglmm, color=miRglmnb.y))+geom_point()+geom_abline()+xlab('NB GLM estimate')+ylab('miRglmm estimate')+scale_colour_discrete(drop=FALSE)+labs(title="NB GLM vs miRglmm", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$miRglmnb.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$miRglmnb,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)
ggsave("figures/figure3A.tif", plot=last_plot(), device="tiff", width=200, height=150, units="mm", dpi=320, bg="white")

ggplot(comb_df, aes(x=DESeq2.x, y=miRglmm, color=miRglmnb.y))+geom_point()+geom_abline()+xlab('DESeq2 estimate')+ylab('miRglmm estimate')+scale_colour_discrete(drop=FALSE)+labs(title="DESeq2 vs miRglmm", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$DESeq2.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$DESeq2,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)
ggsave("figures/figure3B.tif", plot=last_plot(), device="tiff", width=200, height=150, units="mm", dpi=320, bg="white")

ggplot(comb_df, aes(x=edgeR.x, y=miRglmm, color=miRglmnb.y))+geom_point()+geom_abline()+xlab('edgeR estimate')+ylab('miRglmm estimate')+scale_colour_discrete(drop=FALSE)+labs(title="edgeR vs miRglmm", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$edgeR.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$edgeR,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)
ggsave("figures/figure3C.tif", plot=last_plot(), device="tiff", width=200, height=150, units="mm", dpi=320, bg="white")

ggplot(comb_df, aes(x=limmavoom.x, y=miRglmm, color=miRglmnb.y))+geom_point()+geom_abline()+xlab('limma-voom estimate')+ylab('miRglmm estimate')+scale_colour_discrete(drop=FALSE)+labs(title="limma-voom vs miRglmm", color="Significance")+
  annotate("text", x=Inf, y=-Inf, label=paste0(paste("ICC =", format(round(icc_out$limmavoom.x,2), nsmall=2)), '\n', paste("Kappa =", format(round(kappa_out$limmavoom,2), nsmall=2))), hjust=1.05, vjust=-0.3, fontface=2, size=5)
ggsave("figures/figure3D.tif", plot=last_plot(), device="tiff", width=200, height=150, units="mm", dpi=320, bg="white")


library(UpSetR)
list_Input=list(sequence=f1_results_neg10$miRNA[f1_results_neg10$significant==1], DESeq2=deseq2_results$miRNA[deseq2_results$significant==1],`NB GLM`=miRNA_results$miRNA[miRNA_results$significant==1], 
                edgeR=edgeR_results$miRNA[edgeR_results$significant==1], `limma-voom`=limmavoom_results$miRNA[limmavoom_results$significant==1])
upset(fromList(list_Input), order.by="freq")

#find difference in significant miRNA
table(df_merge$significant.x, df_merge$significant.y)
out_disp=subset(df_merge, significant.x!=significant.y)

#find significant miRNA with different direction of sign between two methods
test_sub=subset(df_merge, (sign(col_grouptestes.x)!=sign(col_grouptestes.y))&significant.x==1&significant.y==1)
