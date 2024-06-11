\name{miRglmm}
\alias{miRglmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to fit a generalized linear mixed model (GLMM) on isomiR-level miRNA-seq data. 
}
\description{
Fits a negative binomial GLMM to each miRNA using isomiR-level counts from miRNA-seq data and returns the fitted model objects.
}
\usage{
miRglmm(se, 
        col_group,
        min_med_lcpm = -1, 
        ncores = 1, 
        adjust_var=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{se}{a Summarized Experiment object. Assay(se) must have column names that indicate unique sample IDs. rowData(se) must contain a column called "miRNA" that indicates unique miRNA IDs and a column called "uniqueSequence" that indicates unique isomiR ids}
  \item{col_group}{a vector with the same length as the number of columns of se denoting the fixed effect of interest.}
  \item{min_med_lcpm}{threshold for isomiR level expression based on the minimum median log counts per million. isomiRs with expression below this threshold are filtered prior to model fitting.}
  \item{ncores}{number of cores to use. defaults to non-parallel evaluation.}
  \item{adjust_var}{a vector with the same length as the number of columns of se denoting a variable to adjust for in the model as a fixed effect. This could be a variable from colData(se).}
}
\details{
Model fitting is performed used the glmer.nb function from lme4. The model includes a fixed effect for col_group, and adjust_var, if specified. The model includes a random intercept for sample (defined by colnames(assay(se))), and a random intercept and a random col_group slope for uniqueSequence. A reduced model without the random col_gropu slope term for uniqueSequence is also fit, to allow calculating the likelihood ratio test for the random slope parameter.
}
\value{
Two lists of glmerMod model fit objects called "miRglmm" and "miRglmm_reduced", each with length equal to the number of miRNA modelled. The elements of the lists are named by miRNA (based on unique entries in rowData(se)$miRNA). The "miRglmm" list contains results from the full miRglmm model (including random col_group slope for uniqueSequence). The "miRglmm_reduced" list contains results from the reduced miRglmm model (excluding the random col_group slope for uniqueSequence).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Andrea Baran <andrea_baran@urmc.rochester.edu>, Matthew N. McCall <mccallm@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("se_example")
fit = miRglmm(se_example, col_group=rep(c("A","B"), c(3,3)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
