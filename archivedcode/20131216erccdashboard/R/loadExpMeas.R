loadExpMeas<- function(expDat, countTable, totalReads = NULL){
  
  designFactors <- c("Sample","Rep")
  
  # Check for countable input errors
  for (i in 2:length(colnames(countTable))){
      if (str_count(colnames(countTable[c(i)]),"_") != 1){
        stop("Check countTable column names for use of underscore (_)")
      }
  }
  if (anyDuplicated(colnames(countTable)) > 0){
    stop("Column names must be unique sample IDs")  
  }
  sampleInfo = expDat$sampleInfo
  # Pull idCols out of expDat
  idCols <- sampleInfo$idCols 

  if (sampleInfo$totalSeqReads == F){
    expDat$totalReads = NULL
  }else{
    expDat$totalReads = totalReads
  }
  
  Transcripts = countTable
  
  # Import data based on analysis type from SEQC main project
    
   # force names to be ERCC- and first column name to Feature
   names(Transcripts)[1] = "Feature"
   Transcripts$Feature = gsub("ERCC_","ERCC-",Transcripts$Feature)
   #Transcripts$Feature = gsub(".","-",Transcripts$Feature)
   Transcripts$Feature = gsub(":","",Transcripts$Feature)
   
   # get data frames with just the ERCCs and just the human genes
   TranscriptsERCCOnly = Transcripts[c(grep("ERCC-0", Transcripts$Feature)),]
   TranscriptsHumanOnly = Transcripts[-c(grep("ERCC-0", Transcripts$Feature)),]
   
   # Remove ERCCs in the definition file that are not in the count data file
   idCols = idCols[match(TranscriptsERCCOnly$Feature,idCols$Feature),]
   
   # Remove ERCCs without a Ratio
   idCols = idCols[which(is.finite(idCols$Ratio)),]
   
   # Remove ERCCs from count data and idCols that are absent from the experiment
   TranscriptsERCCOnly = TranscriptsERCCOnly[match(idCols$Feature,
                                                  TranscriptsERCCOnly$Feature),]
   Transcripts = rbind(TranscriptsERCCOnly, TranscriptsHumanOnly)
   
   #############################################################################
   sample1 <- sampleInfo$sample1Name 
   sample2 <- sampleInfo$sample2Name

#   # idxsample <- which((rowMeans(Transcripts[-c(1)])>1)&(rowSums(Transcripts[-c(1)]!=0)>=2))
#    
#   zerofilterIdx = NULL
#   for (i in 1:length(Transcripts$Feature)){
#     idxAdd  = which(Transcripts[i,-c(1)] == 0)
#     if (length(idxAdd) >= 2){
#       zerofilterIdx = append(zerofilterIdx, i) 
#     }
#     
#   }
#   print(head(zerofilterIdx))
#   
#   idxsample = zerofilterIdx
#    Transcripts <- Transcripts[-idxsample,]
#    
#    Transcripts$Feature <- as.factor(as.character(Transcripts$Feature))
#    
#   lowfilterIdx = which(rowSums(Transcripts[-c(1)])<ncol(Transcripts[-c(1)]))
#   print(lowfilterIdx)
# #   lowfilterIdx = NULL
# #   for (i in 1:length(Transcripts$Feature)){
# #     idxAdd  = which(rowsum(Transcripts[i,-c(1)] == 0)
# #     if (length(idxAdd) >= 2){
# #       filterIdx = append(filterIdx, i) 
# #     }
# #     
# #   }
# #   print(head(filterIdx))
# #   
# #   idxsample = filterIdx
# #   Transcripts <- Transcripts[-idxsample,]
#   
#   Transcripts$Feature <- as.factor(as.character(Transcripts$Feature))
#   
#   
#   
#   
#   
#   
#    measERCCs <- Transcripts$Feature[grep("ERCC-0", Transcripts$Feature)]
#    
#    insuffDat <- setdiff(idCols$Feature, measERCCs)
#    
#    print(paste("Transcripts were removed with a mean count < 1 or more than 2",
#                "replicates with 0 counts.")) 
#    print(paste("A total of",length(insuffDat),"out of",length(idCols$Feature),
#                "ERCC controls were filtered"))
#    print("The excluded ERCCs are:")
#    print(insuffDat)
#    print(paste("The remaining",length(measERCCs),"ERCC controls were analyzed"))   
#   
  designMat <- getDesignMat(expressionData = Transcripts,
                               factorList = designFactors,
                               patternSplit = '_')
  # write Transcript csv file to directory
  #write.csv(Transcripts, paste(sampleInfo$filenameRoot,"Transcripts.csv",sep="."),
  #          row.names = F)
  # collect everything to add to expDat
  expDat = append(expDat, list(Transcripts = Transcripts,
                               designMat = designMat,
                               sampleNames = c(sample1,sample2),
                               idCols = idCols,
                               totalReads = totalReads))
  return(expDat)

}