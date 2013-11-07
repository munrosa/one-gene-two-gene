estMnLibeFactor <- function( expDat, cnt = expDat$TranscriptsAB){
  sampleInfo = expDat$sampleInfo
  if(sampleInfo$totalSeqReads == F){
    sampleLibeSums = colSums(cnt[-c(1)],na.rm =T)
    mnLibeFactor = (mean(as.vector(sampleLibeSums)))/10^6
    print("Using total mapped reads mean library size = ")
    print(mnLibeFactor)
  }else{
    mnLibeFactor = (mean(as.vector(expDat$totalReads)))/10^6
    print("Using total sequencing reads mean library size = ")
    print(mnLibeFactor)
    sampleLibeSums = expDat$totalReads
  }
  expDat$mnLibeFactor <- mnLibeFactor
  expDat$sampleLibeSums <- sampleLibeSums
  return(expDat)
}
  
  