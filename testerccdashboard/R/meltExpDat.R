meltExpDat <- function(expDat, cnt, designMat){
  sampleInfo <- expDat$sampleInfo
  libeSize <- expDat$sampleLibeSums
  datNames <- colnames(designMat)[-1]
  sample1 <- expDat$sample1
  sample2 <- expDat$sample2
  
  datCols = cnt[-c(1)]
  libAdjust = sweep(datCols, 2, libeSize,"/")
  sampleLibeDataNorm = cbind(cnt[c(1)],libAdjust)
  myDataERCC = sampleLibeDataNorm
  expressDat = myDataERCC[-c(1)] 
  sampleNameList = c(sample1,sample2)
  
  expressDatAandBERCCs = merge(expDat$idColsAdj[c(1,4)],myDataERCC)
  
  expressDatAandBERCCs$Feature <- as.factor(as.character(
    expressDatAandBERCCs$Feature))
  expressDatAandBERCCs$Ratio <- as.factor(as.character(
    expressDatAandBERCCs$Ratio))
  
  expressDat_l <- melt(expressDatAandBERCCs)
  
  colAdd <- colsplit(expressDat_l$variable,pattern="_",names=datNames)
  colAdd = as.data.frame(lapply(colAdd,as.factor))
  expressDat_l<- expressDat_l[,-c(3)]
  colnames(expressDat_l)[3]<-"NormCounts"
  expressDat_l <- cbind(expressDat_l, colAdd)
  
  #expressDat_l$Sample <- as.factor(as.character(expressDat_l$Sample))
  #expressDat_l$Replicate <- as.factor(as.character(expressDat_l$Replicate))
  expDat$expressDat_l <- expressDat_l
  return(expDat)
}