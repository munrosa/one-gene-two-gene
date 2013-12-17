\name{erccROC}
\alias{erccROC}
\title{Produce Receiver Operator Characteristic (ROC) Curves and AUC statistics}
\usage{
erccROC(expDat, pValDat = NULL)
}
\arguments{
  \item{expDat}{list, contains input data and stores
  analysis results}

  \item{pValDat}{string, Pval csv file name, default is
  NULL}
}
\description{
Produce Receiver Operator Characteristic (ROC) Curves and
AUC statistics
}
\examples{
# Requires P-value data in a csv file as an input with the naming format
# "filenameRoot ERCC Pvals.csv", if the pValDat input argument is not given
# the function will search the working directory for the csv file

expDat <- erccROC(expDat)
}

