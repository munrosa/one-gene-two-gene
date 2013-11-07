saveResults <- function(expDat, erccROC.res, maPlotAB, lodr.ERCC){
  
  filenameRoot <- expDat$sampleInfo$filenameRoot
  
  # Name and consolidate metrics for the interlaboratory comparison
  #nam <- paste(filenameRoot, "expDat",sep = ".")
  #assign(nam,expDat)
  
  nam <- paste(filenameRoot, "AUC",sep = ".")
  assign(nam,erccROC.res$AUCdat)
  
  nam <- paste(filenameRoot, "lodr.ERCC",sep = ".")
  assign(nam,lodr.ERCC)
  
  nam <- paste(filenameRoot, "sdGlobal",sep = ".")
  assign(nam,maPlotAB$sdGlobal)
  
  nam <- paste(filenameRoot, "modRatVar",sep = ".")
  assign(nam,maPlotAB$modRatVar)
  
  nam <- paste(filenameRoot, "r_m",sep = ".")
  assign(nam,expDat$r_m.res$r_m.mn)
  
  nam <- paste(filenameRoot, "r_m_lower",sep = ".")
  assign(nam,expDat$r_m.res$r_m.lower)
  
  nam <- paste(filenameRoot, "r_m_upper",sep = ".")
  assign(nam,expDat$r_m.res$r_m.upper)
  
  nam <- paste(filenameRoot, "p.thresh",sep = ".")
  assign(nam,expDat$p.thresh)
  
  to.save <- ls()
  
  save(list = to.save[grepl(pattern = filenameRoot,x=to.save)],
       file=paste(filenameRoot,"Results","RData",sep = "."))
  
  if (expDat$sampleInfo$printPDF == T){
    dev.off()  
  }
  
}