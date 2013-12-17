\name{geneExprTest}
\alias{geneExprTest}
\title{Prepare differential expression testing results for spike-in analysis}
\usage{
geneExprTest(expDat)
}
\arguments{
  \item{expDat}{list, contains input data and stores
  analysis results}
}
\description{
Prepare differential expression testing results for
spike-in analysis
}
\examples{
# After running initDat and est_r_m functions provide the resulting expDat
# list as input to the geneExprTest function

expDat <- geneExprTest(expDat)
}

