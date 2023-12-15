library(ggplot2)
library(lme4)


## read in results
load(file='bladder_testes_results/filter_neg1_results.rda')

#############################panel A
miRNA_plot="hsa-let-7g-5p"

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

#pull out individual sequence estimates
groupA=fixef(f1)[1]
groupB=fixef(f1)[1]+fixef(f1)[2]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_overall=rbind(groupA, groupB)

miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
groupA=coef(miRNA_sub)[1]
groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_miRNA=rbind(groupA, groupB)
ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm random sequence estimates"))+geom_line()+geom_point()+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm fixed effect estimates"))+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="NB GLM estimates"))+xlab("tissue")+
  scale_color_manual(name="", values=c("miRglmm random sequence estimates"="black", "miRglmm fixed effect estimates"="red", "NB GLM estimates"="blue"))
ggsave(paste0("figures/figure6A_", miRNA_plot, ".tif"), plot=last_plot(), device="tiff", width=250, height=150, units="mm", dpi=320, bg="white")

#############################panel B
miRNA_plot="hsa-miR-26a-5p"

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

#pull out individual sequence estimates
groupA=fixef(f1)[1]
groupB=fixef(f1)[1]+fixef(f1)[2]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_overall=rbind(groupA, groupB)

miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
groupA=coef(miRNA_sub)[1]
groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_miRNA=rbind(groupA, groupB)
ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm random sequence estimates"))+geom_line()+geom_point()+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm fixed effect estimates"))+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="NB GLM estimates"))+xlab("tissue")+
  scale_color_manual(name="", values=c("miRglmm random sequence estimates"="black", "miRglmm fixed effect estimates"="red", "NB GLM estimates"="blue"))
ggsave(paste0("figures/figure6B_", miRNA_plot, ".tif"), plot=last_plot(), device="tiff", width=250, height=150, units="mm", dpi=320, bg="white")

#############################panel C
miRNA_plot="hsa-let-7a-5p"

f1=fits[["miRglmm"]][[miRNA_plot]]

#pull out individual sequence estimates
groupA=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`
groupB=fixef(f1)[1]+ranef(f1)$sequence$`(Intercept)`+ranef(f1)$sequence$col_grouptestes+fixef(f1)[2]
groupA=data.frame("estimate"=groupA, "sequence"=rownames(ranef(f1)$sequence))
groupB=data.frame("estimate"=groupB, "sequence"=rownames(ranef(f1)$sequence))
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df=rbind(groupA, groupB)

#pull out individual sequence estimates
groupA=fixef(f1)[1]
groupB=fixef(f1)[1]+fixef(f1)[2]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_overall=rbind(groupA, groupB)

miRNA_sub=fits[["miRglmnb"]][[miRNA_plot]]
groupA=coef(miRNA_sub)[1]
groupB=coef(miRNA_sub)[1]+coef(miRNA_sub)[2]
groupA=data.frame("estimate"=groupA)
groupB=data.frame("estimate"=groupB)
groupA$Pool="Bladder"
groupB$Pool="Testes"
comb_df_miRNA=rbind(groupA, groupB)
ggplot(comb_df, aes(x=Pool, y=estimate, group=sequence, color="miRglmm random sequence estimates"))+geom_line()+geom_point()+
  geom_line(data=comb_df_overall, aes(x=Pool, y=estimate, group=1, color="miRglmm fixed effect estimates"))+
  geom_line(data=comb_df_miRNA, aes(x=Pool, y=estimate, group=1, color="NB GLM estimates"))+xlab("tissue")+
  scale_color_manual(name="", values=c("miRglmm random sequence estimates"="black", "miRglmm fixed effect estimates"="red", "NB GLM estimates"="blue"))
ggsave(paste0("figures/figure6C_", miRNA_plot, ".tif"), plot=last_plot(), device="tiff", width=250, height=150, units="mm", dpi=320, bg="white")

