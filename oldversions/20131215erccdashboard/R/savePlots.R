savePlots<-function(expDat,plotsPerPg = "manuscript", plotlist = NULL){
  #Options are either the default of printing the plots as shown in publication
  # plotsPerPg = "manuscript" and plotlist is NULL or plotsPerPg = "single" and
  # any combination of the plots can be printed, one per page
  # Open PDF file to write results
  filenameUse <- expDat$sampleInfo$filenameRoot 
  if (plotsPerPg == "manuscript"){
    cols = 2
    pwidth = 7*cols
    pheight = 7*6/cols
    pdf(file = paste(filenameUse,"pdf",sep="."),width=pwidth,height = pheight)
    multiplot(expDat$Figures$ROCplot,expDat$Figures$plotdynRangeAnnot, 
              expDat$Figures$plotLODR.ERCC,expDat$Figures$plotERCCeffects, 
              expDat$Figures$dispPlot,expDat$Figures$plotRatioAnnot,cols=2)
    dev.off()
  } 
  if (plotsPerPg == "single"){
    if (is.null(plotlist)){
      plotlist = expDat$Figures
    } 
    pdf(file = paste(filenameUse,"pdf",sep="."),onefile=T,width=7,height = 7)
    print(plotlist)
    dev.off()
  }
  
}
