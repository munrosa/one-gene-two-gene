#' Initialize the expDat list
#'
#' @param countTable    data frame, the first column contains names of 
#'                      genes or transcripts (Feature) and the remaining columns
#'                      are counts for sample replicates spiked with ERCC 
#'                      controls
#' @param totalReads    vector of totalReads to use for library size
#'                      normalization
#' @param filenameRoot  string root name for output files
#' @param sample1Name   string name for sample 1 in the gene expression 
#'                      experiment
#' @param sample2Name   string name for sample 2 in the gene expression
#'                      experiment
#' @param ERCCMixes     Name of ERCC mixture design, "Ambion4plexPair" is 
#'                      default, the other option is "AmbionSingle"
#' @param ERCCdilution  unitless dilution factor used in dilution of the Ambion 
#'                      ERCC spike-in mixture solutions 
#' @param spikeVol      volume in microliters of diluted ERCC mix spiked into
#'                      the total RNA samples
#' @param totalRNAmass  mass in micrograms of total RNA spiked with diluted ERCC
#'                      mixtures 
#' @param choseFDR      False Discovery Rate for differential expression testing
#'                      , default is 0.05
#' 
#' @export
#' @examples
#' load(file = system.file("data/SEQC.RatTox.Example.RData",
#'      package = "erccdashboard"))
#'      
#' expDat = initDat(countTable = COH.RatTox.ILM.MET.CTL.countTable, totalReads = 
#'                  COH.RatTox.ILM.MET.CTL.totalReads, filenameRoot = "COH.ILM",
#'                  sample1Name = "MET", sample2Name = "CTL", 
#'                  ERCCMixes = "Ambion4PlexPair", ERCCdilution = 1/100, 
#'                  spikeVol = 1, totalRNAmass = 0.500)
#'                  
#' summary(expDat)

initDat <- function(countTable, totalReads, filenameRoot, sample1Name, 
                    sample2Name, ERCCMixes, ERCCdilution, spikeVol, 
                    totalRNAmass, choseFDR){
  ## These variables may become options for user later, for now they will be defined as
  ## internal default values
  totalSeqReads = T
  libeSizeNorm = T
  myYLimMA = c(-3.5,3.5)
  myXLim = c(-10,15)
  myYLim = NULL
  if (is.null(choseFDR)){
    choseFDR = 0.05
    cat(paste("Default choseFDR = ", choseFDR))
  }
  ##############################
  
#library(ggplot2)
#   library(reshape2)
#   library(plyr)
#   library(scales)
#   library(edgeR)
#   library(locfit)
#   library(QuasiSeq)
#   library(grid)
#   library(gridExtra)
#   library(stringr)

  sampleInfo = list(sample1Name = sample1Name,
                    sample2Name = sample2Name, choseFDR = choseFDR,
                    ERCCdilution = ERCCdilution, spikeVol = spikeVol,
                    totalRNAmass = totalRNAmass, 
                    totalSeqReads = totalSeqReads, 
                    libeSizeNorm = libeSizeNorm, myYLimMA = myYLimMA,
                    myXLim = myXLim, myYLim = myYLim)
  
  expDat = list(sampleInfo = sampleInfo)
  
  if (!is.null(filenameRoot)){
    expDat <- dashboardFile(expDat,filenameRoot = filenameRoot)  
  }else{
    stop("The filenameRoot character string has not been defined!")
  }
  
  
  ###############################################################################
  # Run loadERCCInfo function to obtain ERCC information
  expDat <- loadERCCInfo(expDat, erccMixes=ERCCMixes)
  
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