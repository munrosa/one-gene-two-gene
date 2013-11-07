siteList = c("AGR.NCTR.ILM.A.B.Transcripts.csv","BGI.NCTR.ILM.A.B.Transcripts.csv","COH.NCTR.ILM.A.B.Transcripts.csv","CNL.NCTR.ILM.A.B.Transcripts.csv","MAY.NCTR.ILM.A.B.Transcripts.csv","NVS.NCTR.ILM.A.B.Transcripts.csv","NWU.NIST.LIF.A.B.Transcripts.csv","PSU.NIST.LIF.A.B.Transcripts.csv","SQW.NIST.LIF.A.B.Transcripts.csv")
  
fileList = list()
for (site in 1:length(siteList)){
    fileRead <- read.csv(siteList[site],header = T)
    #print(grep("ERCC-0",fileRead[,c(1)]))
    fileERCC <- as.character(fileRead[grep("ERCC-0",fileRead[,c(1)]),c(1)])
    fileList[[site]] <- fileERCC
} 
print(str(fileList))
#uniqERCC <- unique(fileList)

intersectERCC <- Reduce(intersect, fileList)
#rm(list=setdiff(ls(), "intersectERCC"))