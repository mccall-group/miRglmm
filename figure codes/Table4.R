### summarize simulation results

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

MSE_mat=as.data.frame(matrix(unlist(results[["MSE_sim"]]), ncol=5, byrow=TRUE))
colnames(MSE_mat)=colnames(results[["MSE_sim"]][[1]])

mean_MSE=colMeans(MSE_mat)
std_MSE=apply(MSE_mat,2 ,sd)
median_MSE=apply(MSE_mat, 2, median)
min_MSE=apply(MSE_mat,2,min)
max_MSE=apply(MSE_mat,2,max)
n_min_MSE=apply(MSE_mat,1, which.min)
n_min_MSE=factor(n_min_MSE, levels=c("1", "2", "3", "4", "5"))
out=data.frame(table(n_min_MSE))

MSE_table=data.frame("mean MSE"=mean_MSE, "SD MSE"=std_MSE, "median MSE"=median_MSE, "minimum MSE"=min_MSE, "maximum MSE"=max_MSE, "number of times minimizing MSE"=out[,2])

write.csv(MSE_table, file="figures/Table4.csv")


mean(MSE_mat$miRglmm[which(n_min_MSE==1)]-MSE_mat$DESeq2[which(n_min_MSE==1)])
mean(MSE_mat$DESeq2[which(n_min_MSE==3)]-MSE_mat$miRglmm[which(n_min_MSE==3)])
median(MSE_mat$miRglmm[which(n_min_MSE==1)]-MSE_mat$DESeq2[which(n_min_MSE==1)])
median(MSE_mat$DESeq2[which(n_min_MSE==3)]-MSE_mat$miRglmm[which(n_min_MSE==3)])