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

initDat <- function(countTable, totalReads, designFactors, 
                    sample1Name, sample2Name, erccMixes, 
                    ERCCdilution = 1, spikeVol = 1, totalRNAmass = 1, 
                    printPDF = T, filenameRoot, DEtest = T, 
                    choseFDR= 0.05, totalSeqReads = T, libeSizeNorm = T, 
                    myYLimMA = c(-3.5,3.5), myXLim = c(-10,15),myYLim = NULL){
  library(ggplot2)
  library(reshape2)
  library(plyr)
  library(scales)
  library(edgeR)
  library(locfit)
  library(QuasiSeq)
  library(grid)
  
  sampleInfo = list(designFactors = designFactors, sample1Name = sample1Name,
                    sample2Name = sample2Name, choseFDR = choseFDR,
                    ERCCdilution = ERCCdilution, spikeVol = spikeVol,
                    totalRNAmass = totalRNAmass, printPDF = printPDF, 
                    DEtest = DEtest, totalSeqReads = totalSeqReads, 
                    libeSizeNorm = libeSizeNorm, myYLimMA = myYLimMA,
                    myXLim = myXLim, myYLim = myYLim)
  
  expDat = list(sampleInfo = sampleInfo)
  
  if (!is.null(filenameRoot)){
    expDat <- dashboardPDF(expDat,filenameRoot = filenameRoot)  
  }else{
    expDat <- dashboardPDF(expDat)
  }
  
  
  ###############################################################################
  # Run loadERCCInfo function to obtain ERCC information
  expDat <- loadERCCInfo(expDat, erccMixes="Ambion4plexPair")
  
  ###############################################################################
  # Assume user has created data frame countTable and totalReads vector
  # process those data files to add to expDat structure
  expDat <- loadExpMeas(expDat, countTable, designFactors, totalReads)
  
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