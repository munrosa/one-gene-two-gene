geneExprTest <- function(expDat, cnt, designMat ){
  sampleInfo <- expDat$sampleInfo
  info <- designMat
  choseFDR <- sampleInfo$choseFDR
  qvalFile <- paste(sampleInfo$filenameRoot,"quasiSeq.res.csv",sep=".")
  # First 3 columns of qvalFile must contain Feature, pvals, and qvals
  # Decide to reuse results or run testDE
  if (file.exists(qvalFile) == TRUE){
    cat(paste("\n Differential expression test results exist, will use   \n",
              "existing results for analysis. No dispersion plots will\n",
              "be produced. Delete", qvalFile, "if  \n",
              "you want to repeat differential expression testing.    \n"))
  }else{
    suppressWarnings(testDE(sampleInfo, expDat, cnt = cnt, info = info ))
  }
  deRes <- read.csv(file = paste(sampleInfo$filenameRoot,"quasiSeq.res.csv",
                                 sep="."))
  
  p.thresh<-.1
  
  if(any(deRes$qvals<choseFDR)) p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
  
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
