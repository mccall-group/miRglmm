test_miRglmm <- function(){
    data(se_example)
    checkEquals(class(oncogene2013)[1],"SummarizedExperiment")
    fit = miRglmm(se_example, col_group=rep(c("A","B"), c(3,3)))
    checkEquals(class(fit), "list")
    checkEquals(names(fit), c("miRglmm","miRglmm_reduced"))
    checkEquals(names(fit$miRglmm), c("hsa-miR-19b-3p","hsa-miR-101-3p","hsa-let-7f-5p"))
    checkEqualsNumeric(coef(fit$miRglmm$`hsa-miR-19b-3p`)$sequence[5,2], 
                       0.9218401, tolerance=1.0e-4)
}
