### summarize simulation results

load(file="sim results/sims_N100_m2_s1_rtruncnorm13_combinedresults.rda")

cov_mat=as.data.frame(matrix(unlist(results[["coverage_prob_sim"]]), ncol=4, byrow=TRUE))
colnames(cov_mat)=colnames(results[["coverage_prob_sim"]][[1]])

mean_cov=colMeans(cov_mat)
std_cov=apply(cov_mat,2 ,sd)
median_cov=apply(cov_mat, 2, median)
min_cov=apply(cov_mat,2,min)
max_cov=apply(cov_mat,2,max)


cov_table=data.frame("mean cov"=mean_cov, "SD cov"=std_cov, "median cov"=median_cov, "minimum cov"=min_cov, "maximum cov"=max_cov)

write.csv(cov_table, file="figures/Table5.csv")


