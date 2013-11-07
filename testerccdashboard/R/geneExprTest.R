geneExprTest <- function(expDat, cnt, designMat ){
  sampleInfo <- expDat$sampleInfo
  info <- designMat
  DEtest <- sampleInfo$DEtest
  choseFDR <- sampleInfo$choseFDR
  if (DEtest == TRUE){
    suppressWarnings(testDE(sampleInfo, expDat, cnt = cnt, info = info ))
  }
  if (DEtest == FALSE){
    print(paste("Not testing for differential expression, no dispersion plots",
                "will be produced and pre-existing DE test result files will", 
                "be used for additional analysis"))
  }
  
  deRes <- read.csv(file = paste(sampleInfo$filenameRoot,"quasiSeq.res.csv",sep="."))
  
  p.thresh<-.1
  
  if(any(deRes$qvals<choseFDR)) p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
  
  print("Threshold P-value")
  print(p.thresh)
  
  if (p.thresh > .1){
    print(paste("Threshold P-value is high for the chosen FDR of", 
                as.character(choseFDR)))
    print(paste("The sample comparison indicates a large amount of differential",
                "expression in the measured transcript populations"))
  }
  expDat$p.thresh <- p.thresh
  return(expDat)
}
