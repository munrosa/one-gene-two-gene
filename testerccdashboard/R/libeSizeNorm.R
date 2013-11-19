libeSizeNorm <- function(expDat){
  #attach(expDat)
  sampleInfo = expDat$sampleInfo
  expressDat = expDat$Transcripts
  # Library size normalize the data
  
  if (sampleInfo$totalSeqReads == F){
    TranscriptsAll = expressDat
    TranscriptMappedReadSums = colSums(TranscriptsAll[-c(1)],na.rm = T)
    libeSize = TranscriptMappedReadSums
    datCols = expressDat[-c(1)]
    libeSize = libeSize/(10^6) #per million mapped reads
    #Library size normalize the data  
    libAdjust = sweep(datCols, 2, libeSize,"/")
    expressDat = cbind(expressDat[c(1)], libAdjust)
  }
  if (sampleInfo$totalSeqReads == T){
    TranscriptsAll = expressDat
    libeSize = expDat$totalReads
    datCols = expressDat[-c(1)]
    libeSize = libeSize/(10^6) #per million mapped reads
    #Library size normalize the data  
    libAdjust = sweep(datCols, 2, libeSize,"/")
    expressDat = cbind(expressDat[c(1)], libAdjust)
  }
  
  cat("Library sizes:\n")
  cat(libeSize)
  
  # get just ERCC data in the expression data frame
  expressDat = expressDat[c(grep("ERCC-0", expressDat$Feature)),]
  
  #detach(expDat)
  expDat$expressDat <- expressDat
  expDat$libeSize <- libeSize
  return(expDat)
  
}