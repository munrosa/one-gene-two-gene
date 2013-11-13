#' Combine technical replicates
#'
#' @param type Either "sum" or "mean"
#' @param libeSizeNorm If TRUE (the default), library size normalize the data
#' @param totalSeqReads use the total number of sequence reads in FASTA/FASTQ file for library size normalization
#' @param totalReads use the total number of reads 
#' @param expressDat the expression data
#' @param designMat the design matrix
#' @param sampleNameList the sample names, character list
#' @param libeList the library list, sample characters
#' 
#' @keywords manip
#' @export
#' @examples
#' mtcars[with(mtcars, order(cyl, disp)), ]
#' arrange(mtcars, cyl, disp)
#' arrange(mtcars, cyl, desc(disp))

initDat <- function(countTable, totalReads, filenameRoot, sample1Name, 
                    sample2Name, ERCCdilution, spikeVol, totalRNAmass, choseFDR,
                    printPDF = T){
  ## These variables may become options for user later, for now they will be defined as
  ## internal default values
  erccMixes = "Ambion4plexPair"
  totalSeqReads = T
  libeSizeNorm = T
  myYLimMA = c(-3.5,3.5)
  myXLim = c(-10,15)
  myYLim = NULL
  ##############################
  
  library(ggplot2)
  library(reshape2)
  library(plyr)
  library(scales)
  library(edgeR)
  library(locfit)
  library(QuasiSeq)
  library(grid)
  library(gridExtra)
  library(stringr)

  sampleInfo = list(sample1Name = sample1Name,
                    sample2Name = sample2Name, choseFDR = choseFDR,
                    ERCCdilution = ERCCdilution, spikeVol = spikeVol,
                    totalRNAmass = totalRNAmass, printPDF = printPDF, 
                    totalSeqReads = totalSeqReads, 
                    libeSizeNorm = libeSizeNorm, myYLimMA = myYLimMA,
                    myXLim = myXLim, myYLim = myYLim)
  
  expDat = list(sampleInfo = sampleInfo)
  
  if (!is.null(filenameRoot)){
    expDat <- dashboardPDF(expDat,filenameRoot = filenameRoot)  
  }else{
    stop("The filenameRoot character string has not been defined!")
  }
  
  
  ###############################################################################
  # Run loadERCCInfo function to obtain ERCC information
  expDat <- loadERCCInfo(expDat, erccMixes="Ambion4plexPair")
  
  ###############################################################################
  # Assume user has created data frame countTable and totalReads vector
  # process those data files to add to expDat structure
  expDat <- loadExpMeas(expDat, countTable, totalReads)
  
  ###############################################################################
  # library size normalize the data
  expDat <- libeSizeNorm(expDat)
  
  ###############################################################################
  # length normalize the ERCC concentrations
  expDat <- prepERCCDat(expDat)
  
  # Estimate the mean library size factor for the data to use to estimate
  # corresponding concentrations for LODR
  expDat <- estMnLibeFactor(expDat, cnt = expDat$Transcripts)
  
  return(expDat)
  
}