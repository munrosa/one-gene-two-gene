annotLODR <- function(expDat){
  
  LODR.annot.ERCC <- printLODRres(expDat)
  
  expDat <- dynRangePlotLODR(dynRangeRes = expDat$Figures$plotdynRange,
                   LODR.annot.ERCC = LODR.annot.ERCC)
  
  expDat<- maConcPlot(expDat, LODR.annot.ERCC, alphaPoint = 0.8, r_mAdjust = T, 
                        replicate = T)
  
  return(expDat)
}