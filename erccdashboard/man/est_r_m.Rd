\name{est_r_m}
\alias{est_r_m}
\title{Estimate the mRNA fraction differences for the pair of samples using
replicate data}
\usage{
est_r_m(expDat)
}
\arguments{
  \item{expDat}{list, contains input data and stores
  analysis results}
}
\description{
Estimate the mRNA fraction differences for the pair of
samples using replicate data
}
\examples{
# After running initDat to initialize the expDat list run the following
# command to estimate r_m
expDat <- est_r_m(expDat)

# View the mRNA fraction results
expDat$r_m.res

# View the per ERCC r_m estimates
expDat$Figures$r_mPlot
}

