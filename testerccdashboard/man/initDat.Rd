\name{initDat}
\alias{initDat}
\title{Combine technical replicates}
\usage{
  initDat(countTable, totalReads, designFactors,
    sample1Name, sample2Name, erccMixes, ERCCdilution = 1,
    spikeVol = 1, totalRNAmass = 1, printPDF = T,
    filenameRoot, DEtest = T, choseFDR = 0.05,
    totalSeqReads = T, libeSizeNorm = T,
    myYLimMA = c(-3.5, 3.5), myXLim = c(-10, 15),
    myYLim = NULL)
}
\arguments{
  \item{type}{Either "sum" or "mean"}

  \item{libeSizeNorm}{If TRUE (the default), library size
  normalize the data}

  \item{totalSeqReads}{use the total number of sequence
  reads in FASTA/FASTQ file for library size normalization}

  \item{totalReads}{use the total number of reads}

  \item{expressDat}{the expression data}

  \item{designMat}{the design matrix}

  \item{sampleNameList}{the sample names, character list}

  \item{libeList}{the library list, sample characters}
}
\description{
  Combine technical replicates
}
\examples{
mtcars[with(mtcars, order(cyl, disp)), ]
arrange(mtcars, cyl, disp)
arrange(mtcars, cyl, desc(disp))
}
\keyword{manip}

