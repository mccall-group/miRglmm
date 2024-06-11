\name{se_example}
\alias{se_example}
\docType{data}
\title{
Example data set to demonstrate miRglmm functionality.
}
\description{
A SummarizedExperiment object containing a subset of miRNA-seq data from 6 monocyte samples. 
}
\usage{data("se_example")}
\format{
  A Summarized experiment with isomiR-level read counts for 3 miRNAs and 6 monocyte samples. 
}
\details{
The column names of assay(se) indicate unique sample identifiers. rowData(se) contains columns called "miRNA" and "uniqueSequence" as identifiers of miRNA and isomiR for each read count. 
}
\source{
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100467
}
\references{
Juzenas S, Venkatesh G, Hübenthal M, Hoeppner MP et al. A comprehensive, cell specific microRNA catalogue of human peripheral blood. Nucleic Acids Res 2017 Sep 19;45(16):9290-9301. PMID: 28934507
}
\examples{
data(se_example)
show(se_example)
}
\keyword{datasets}
