\name{initDat}
\alias{initDat}
\title{Initialize the expDat list}
\usage{
initDat(countTable, totalReads, filenameRoot, sample1Name, sample2Name,
  ERCCMixes, ERCCdilution, spikeVol, totalRNAmass, choseFDR)
}
\arguments{
  \item{countTable}{data frame, the first column contains
  names of genes or transcripts (Feature) and the remaining
  columns are counts for sample replicates spiked with ERCC
  controls}

  \item{totalReads}{vector of totalReads to use for library
  size normalization}

  \item{filenameRoot}{string root name for output files}

  \item{sample1Name}{string name for sample 1 in the gene
  expression experiment}

  \item{sample2Name}{string name for sample 2 in the gene
  expression experiment}

  \item{ERCCMixes}{Name of ERCC mixture design,
  "Ambion4plexPair" is default, the other option is
  "AmbionSingle"}

  \item{ERCCdilution}{unitless dilution factor used in
  dilution of the Ambion ERCC spike-in mixture solutions}

  \item{spikeVol}{volume in microliters of diluted ERCC mix
  spiked into the total RNA samples}

  \item{totalRNAmass}{mass in micrograms of total RNA
  spiked with diluted ERCC mixtures}

  \item{choseFDR}{False Discovery Rate for differential
  expression testing , default is 0.05}
}
\description{
Initialize the expDat list
}
\examples{
load(file = system.file("data/SEQC.RatTox.Example.RData",
     package = "erccdashboard"))

expDat = initDat(countTable = COH.RatTox.ILM.MET.CTL.countTable, totalReads =
                 COH.RatTox.ILM.MET.CTL.totalReads, filenameRoot = "COH.ILM",
                 sample1Name = "MET", sample2Name = "CTL",
                 ERCCMixes = "Ambion4PlexPair", ERCCdilution = 1/100,
                 spikeVol = 1, totalRNAmass = 0.500)

summary(expDat)
}

