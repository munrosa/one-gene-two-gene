#' Prepare differential expression testing results for spike-in analysis
#'
#' @param expDat    list, contains input data and stores analysis results
#' 
#' @export
#' @examples
#' # After running initDat and est_r_m functions provide the resulting expDat 
#' # list as input to the geneExprTest function
#' 
#' expDat <- geneExprTest(expDat)
#' 
geneExprTest <- function(expDat){
  
  cnt <- expDat$Transcripts
  designMat <- expDat$designMat
  sampleInfo <- expDat$sampleInfo
  info <- designMat
  choseFDR <- sampleInfo$choseFDR
  
  # set initial p.thresh
  p.thresh<-.1
  # First 3 columns of qvalFile must contain Feature, pvals, and qvals
  
  qvalFile <- paste(sampleInfo$filenameRoot,"quasiSeq.res.csv",sep=".")
  pvalERCC <- paste(sampleInfo$filenameRoot, "ERCC Pvals.csv",sep=" ")
  #dispPlotFile = paste(sampleInfo$filenameRoot, ".DispPlot.RData",sep = "")
  
    # Decide to reuse results or run testDE
  if (file.exists(qvalFile) == TRUE){
    deRes <- read.csv(qvalFile)
    if (names(deRes)[3] != "qvals"){
      stop("qvals column is missing in qvalFile")
    }
    if (file.exists(pvalERCC) == TRUE){
    cat(paste("\n Differential expression test results exist, will use   \n",
              "existing P-values and Q-values for analysis.\n",
              "Delete", qvalFile, "if you want to repeat differential \n",
              "differential expression testing or view dispersion plots.    \n"))
      if(any(deRes$qvals<choseFDR)){
        p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
      }
    }else{
      cat("No P-values exist for ERCCs, starting differential expression tests")
      expDat <- suppressWarnings(testDE(sampleInfo, expDat, cnt = cnt, 
                                        info = info ))  
      deRes <- read.csv(qvalFile)
      if(any(deRes$qvals<choseFDR)){
        p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
      }
    }
    
  }
  if(file.exists(qvalFile) == FALSE){
    pvalERCC = paste(sampleInfo$filenameRoot, "ERCC Pvals.csv",sep=" ")
    if(file.exists(pvalERCC)){
      cat(paste("\n P-values exist for ERCCs, but Q-values are missing\n"))
      # First three columns must be Feature, MnCnt, Pval
      # Need to get Threshold p-value from user
      getPThresh<- function(){
        cat("\nTo continue with P-values for LODR estimation\n")
        readline("Enter the threshold P-value: ")
      }
      p.thresh <- as.numeric(getPThresh())
      
    }else{
      cat("No P-values exist for ERCCs, starting differential expression tests")
      expDat <- suppressWarnings(testDE(sampleInfo, expDat, cnt = cnt, 
                                        info = info ))  
      deRes <- read.csv(qvalFile)
      if(any(deRes$qvals<choseFDR)){
        p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
      } 
    }
    
  }


  if(is.null(expDat$Figures$dispPlot)){
      cat(paste("\nDE testing results supplied without companion dispersion\n",
              "plot. Dispersion plot is unavailable to print.\n"))
  }  
  cat("\nThreshold P-value\n")
  cat(p.thresh,"\n")
  
  if (p.thresh > .1){
    cat(paste("Threshold P-value is high for the chosen FDR of ", 
                as.character(choseFDR)))
    cat(paste("The sample comparison indicates a large amount of \n",
              "differential expression in the measured transcript \n",
              "populations\n"))
  }
  
  
  expDat$p.thresh <- p.thresh
  return(expDat)
  
}
