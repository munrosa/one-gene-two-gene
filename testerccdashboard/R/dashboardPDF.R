dashboardPDF <-function(expDat, filenameRoot){
  sampleInfo <- expDat$sampleInfo
 
  filenameUse = with(sampleInfo, paste(filenameRoot,sample1Name,
                                          sample2Name,sep = "."))  
  print(filenameUse)
  theme_set(theme_bw(base_size=12))
  # Open PDF file to write results
  if (sampleInfo$printPDF == T){
    pdf(file = paste(filenameUse,"pdf",sep="."),onefile=T,width=7,height = 7)  
  }
  
  expDat$sampleInfo$filenameRoot <- filenameUse
  
  return(expDat)
}