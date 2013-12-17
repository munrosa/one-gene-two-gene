##############################################################################
# Analysis of per site results from SEQC from a given analysis pipeline
##############################################################################

### Run this commented code line here and then source the rest of the script
# siteName = "COH";platform = "Array";analysis = "COH";libeSizeNorm = F;printPDF = F; DEtest = T; totalSeqReads = T; sample1 = "A";sample2 = "B";  FCcode = data.frame(Ratio = c("a","b","c","d"),FC =  c(4,1,.667,.5)); myYLimMA = c(-3,3); myXLim = c(-7,18); choseFDR = 0.05; legendLabels = c("4:1","1:1","1:1.5","1:2")


#####
# Input parameter adjustments
# For platform = "ILM" use analysis = "NCTR", siteName may be equal to "MAY","BGI","CNL","NVS","AGR",or "COH"
# For platform = "LIF" use analysis = "NIST", siteName may be equal to "SQW","NWU",or "PSU"

###############################################################################
# Examples of alternative SEQC ratios to compare, substitute for FCcode above if needed
# A:B ratios are default FCcode matrix
# A:C 
# FCcode = data.frame(Ratio = c("a","b","c","d"),FC = c(1.23,1,.89,.80))
# B:D 
# FCcode = data.frame(Ratio = c("a","b","c","d"),FC = c(0.57,1,1.09,1.14))


###############################################################################
# Set all warnings to report as they appear
options(warn = 1)

# Get required libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(edgeR)
library(DESeq)
library(QuasiSeq)
library(grid)

# Source package scripts
source("R/combineTechReps.R") # combine technical replicates
source("R/compareReps.R") # compare technical replicates for a sample 
source("R/designMat.R") # build factor matrix from sample names
source("R/dynRangePlot.R") # dose response curve for dynamic range
source("R/est_LODR.R") # estimate LODR
source("R/est_r_m.R") # estimate r_m
source("R/loadERCCInfo.R") # load ERCC definition files
source("R/loadExpMeas.R") # load expression measures from file
source("R/maConcPlot.R") # MA plot coded with LODR
source("R/multiplot.R") # plot multiple plots on the same page
source("R/rocPlot.R") # ROC Plot
source("R/testDE.R") # test for DE, obtain p and q values

###############################################################################################################

# Run loadExpMeas function to process count data file
Transcripts <- loadExpMeas(siteName = siteName, platform = platform, analysis = analysis)

###############################################################################################################

# Run loadERCCInfo function to obtain ERCC information
erccInfo <- loadERCCInfo()
idCols = erccInfo$idCols
MixDef = erccInfo$MixDef

###############################################################################################################

# get data frames with just the ERCCs and just the human genes
TranscriptsERCCOnly = Transcripts[c(grep("ERCC-0", Transcripts$Feature)),]
TranscriptsHumanOnly = Transcripts[-c(grep("ERCC-0", Transcripts$Feature)),]

# Get Column sums for all Transcripts, ERCCs + endogenous transcripts
TranscriptMappedReadSums = colSums(Transcripts[-c(1)], na.rm=T)
# Get Column sums for only the ERCCs 
ERCCMappedReadSums = colSums(TranscriptsERCCOnly[-c(1)],na.rm = T)


# Check for and remove ERCCs in the definition file that are not in the count data file
idCols = idCols[match(TranscriptsERCCOnly$Feature,idCols$Feature),]

# Remove ERCCs without a Ratio
idCols = idCols[which(is.finite(idCols$Ratio)),]

# Remove ERCCs from count data and idCols that are not present in the experiment
TranscriptsERCCOnly = TranscriptsERCCOnly[match(idCols$Feature,TranscriptsERCCOnly$Feature),]
Transcripts = rbind(TranscriptsERCCOnly, TranscriptsHumanOnly)



###############################################################################################################

# Generate the Design Matrix for the table example name is SEQC_ILM_BGI_A_1_L01_ATCACG_AC0AYTACXX
if (analysis == "NCTR"){
  designMat = generateDesignMat(TranscriptsERCCOnly, factorList = c("Study","Platform","Site","Sample","Library","Lane","Barcode","Flowcell"), patternSplit = '_')
}
if(analysis == "NIST"){
 designMat = generateDesignMat(TranscriptsERCCOnly, factorList = c("Study","Platform","Site","Sample","Library","Lane","Barcode","Flowcell","SiteInfo1","SiteInfo2"), patternSplit = '_')
}

###############################################################################################################
# Get total Reads (from the original fastq files)
if(totalSeqReads == T){
mainReads = read.delim(paste("Data/MainSEQC/all_qc_results/",platform,"_",siteName,"_qc_results.txt", sep=""))
mainReads = mainReads[c(1:9)]

}

###############################################################################################################
filenameRoot = paste(siteName,analysis,platform,sample1,sample2,sep = ".")
# Open PDF file to write results
if (printPDF == T){
  pdf(file = paste(filenameRoot,"%03d","pdf",sep="."),onefile=T,width=7,height = 7)  
}

###############################################################################################################

# Subset to get just the A and B samples for the Count Matrix and for the designMatrix
select = subset(designMat, (Sample == sample1)|(Sample == sample2))
select <- as.data.frame(lapply(select,as.character))
select <- as.data.frame(lapply(select,as.factor))
designMatAB = select
dataAB = Transcripts[-c(1)]
TranscriptsAB = cbind(Transcripts[c(1)],dataAB[c(match(select$countSet, names(dataAB)))])

# Change the sample names to UHRR and HBRR first
colnames(TranscriptsAB) <- gsub(pattern = "_A_",replacement="_UHRR_",x=colnames(TranscriptsAB))
colnames(TranscriptsAB) <- gsub(pattern = "_B_",replacement="_HBRR_",x=colnames(TranscriptsAB))
sample1 = "UHRR";sample2 = "HBRR";

if (analysis == "NCTR"){
  designMatAB = generateDesignMat(TranscriptsAB, factorList = c("Study","Platform","Site","Sample","Library","Lane","Barcode","Flowcell"), patternSplit = '_')
}
if(analysis == "NIST"){
  designMatAB = generateDesignMat(TranscriptsAB, factorList = c("Study","Platform","Site","Sample","Library","Lane","Barcode","Flowcell","SiteInfo1","SiteInfo2"), patternSplit = '_')
}

# write csv file to directory
write.csv(TranscriptsAB, paste(filenameRoot,"Transcripts.csv",sep="."),row.names = F)


###############################################################################################################

# Library size normalize the data if libeSizeNorm = T
expressDat = TranscriptsAB

if ((libeSizeNorm == T)&(totalSeqReads == F)){
  TranscriptsAll = expressDat
  TranscriptMappedReadSums = colSums(TranscriptsAll[-c(1)],na.rm = T)
  libeSize = TranscriptMappedReadSums
  datCols = expressDat[-c(1)]
  libeSize = libeSize/(10^6) #per million mapped reads
  #Library size normalize the data  
  libAdjust = sweep(datCols, 2, libeSize,"/")
  expressDat = cbind(expressDat[c(1)], libAdjust)
  print("ERCC % of total mapped reads (endo + ERCC) :")
  ERCCofTotal = ERCCMappedReadSums/TranscriptMappedReadSums
  print(summary(ERCCofTotal*100))
  print("Total Mapped Reads Library size summary:")
  print(summary(libeSize))
  }
if ((libeSizeNorm == T)&(totalSeqReads == T)){
  TranscriptsAll = expressDat
  
  totalReads = mainReads$numReads[c(match(select$countSet,mainReads$sampleName))]
  libeSize = totalReads
  datCols = expressDat[-c(1)]
  libeSize = libeSize/(10^6) #per million mapped reads
  #Library size normalize the data  
  libAdjust = sweep(datCols, 2, libeSize,"/")
  expressDat = cbind(expressDat[c(1)], libAdjust)
  print("ERCC % of total sequenced reads:")
  ERCCOnly = expressDat[c(grep("ERCC-0", expressDat$Feature)),]
  ERCCMappedReadSums = colSums(ERCCOnly[-c(1)],na.rm = T)
  
  ERCCofTotal = ERCCMappedReadSums/totalReads
  print(summary(ERCCofTotal*100))
  print("Total Sequenced Reads Library size summary:")
  print(summary(libeSize))
  
}


  # get just ERCC data in the expression data frame
  expressDat = expressDat[c(grep("ERCC-0", expressDat$Feature)),]

  # Length normalize the Expected ERCC concentrations
  lengthFactor = (idCols$Length)/(1000)
  
  # If length normalization of the expected concentrations is desired (default)
  idCols$Conc1 = (idCols$Conc1*lengthFactor)
  idCols$Conc2 = (idCols$Conc2*lengthFactor)  

###############################################################################################################
# Dynamic Range Plots of Sample-Library Combinations
###############################################################################################################

print("Signal-Abundance Plots for dynamic range estimation...")

# Revert to raw counts to sum fluidic replicates and then sum by library replicate
# exclude library 5

sampleLibeData = combineTechReps(type = "sum", libeSizeNorm = T, totalSeqReads= totalSeqReads, totalReads = totalReads, expressDat = TranscriptsAB, designMat = designMatAB, sampleNameList = c(sample1,sample2), libeList = levels(designMatAB$Library)[-5])

sampleLibeDataNorm = sampleLibeData[[1]]

sampleLibeDataNoNorm = sampleLibeData[[2]]

sampleLibeSums = colSums(sampleLibeDataNoNorm[-c(1)],na.rm =T)

if(totalSeqReads==T){
  sampleLibeSums = sampleLibeData[[3]]
  totalReads = sampleLibeSums
} 

# create Data frame of the Sample library pairs to return to the workspace if desired
write.csv(sampleLibeDataNoNorm,file=paste(siteName,analysis,platform,sample1,sample2,"samplelibedataNoLibeSize.csv",sep = "."))

# create Data frame of the Sample library pairs to return to the workspace if desired
write.csv(sampleLibeDataNorm,file=paste(siteName,analysis,platform,sample1,sample2,"samplelibedataNorm.csv",sep = "."))

designMatNorm = generateDesignMat(sampleLibeDataNorm, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')

expressDatSampleMean = combineTechReps(type = "mean", libeSizeNorm = T, totalSeqReads = totalSeqReads, totalReads = totalReads, expressDat = sampleLibeDataNorm, designMat = designMatNorm, sampleNameList = c(sample1,sample2), libeList = levels(designMatAB$Library)[-5])

expressDatSampleMean = expressDatSampleMean[[1]]
expressDatMeanERCC = expressDatSampleMean[c(grep("ERCC-0", expressDatSampleMean$Feature)),]

write.csv(expressDatMeanERCC,file=paste(filenameRoot,"expressDatMeanERCC.csv",sep = "."))


###############################################################################################################
  ### Estimate r_m for the sample pair using a negative binomial glm
  r_m.res = glm_est_r_m(site = filenameRoot, cnt = sampleLibeDataNoNorm, libeSizeNorm = libeSizeNorm,totalSeqReads = totalSeqReads,totalReads = totalReads, idCols = idCols, FCcode = FCcode,sample1 = sample1, sample2 = sample2,legendLabels=legendLabels)
  
#Return the mean estimate of r_m
r_m.mn = r_m.res$r_m.mn
print(r_m.mn)
#Return the upper bound of the 95% confidence interval for r_m
r_m.upper = r_m.res$r_m.upper

#Return the lower bound of the 95% confidence interval for r_m
r_m.lower = r_m.res$r_m.lower

# Return the r_m adjusted concentrations in idCols
idColsAdj = r_m.res$idCols

#expressDatAandBERCCs = merge(idCols[c(1,4)],expressDatSampleAB)
r_m.mnlog = exp(r_m.mn)



###############################################################################################################

### Generate new design matrix
  designMatSum <- generateDesignMat(sampleLibeDataNorm, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')

###############################################################################################################

### Evaluate dynamic range
  dynRangeDat = dynRangePlot(expressDat = sampleLibeDataNorm, designMat = designMatSum,idCols = idColsAdj, sampleNames = c(sample1,sample2),myXLim=myXLim,myYLim=NULL, legendLabels=legendLabels,FCcode=FCcode)
  
  print(dynRangeDat$fit.coeff)
  fit.coeff = dynRangeDat$fit.coeff

###############################################################################################################

### Estimate average totalSeqReads to use for mnLibeFactor estimation
  if(totalSeqReads == F){
  sampleLibeSums = colSums(TranscriptsAB[-c(1)],na.rm =T)
  mnLibeFactor = (mean(as.vector(sampleLibeSums)))/10^6
  print("Using total mapped reads")
  print(mnLibeFactor)
  }else{
  mnLibeFactor = (mean(as.vector(totalReads)))/10^6
  print("Using total sequencing reads")
  print(mnLibeFactor)
  sampleLibeSums = totalReads
  }
###############################################################################################################

  #####
  # After observing that Library 5 is an outlier it is excluded from further analysis
  # The four remaining library replicates for each sample are collected into a data frame with the header names:
  # Feature Ratio Sample1 Sample2 Library
  #####
  myDataERCC = sampleLibeDataNorm
  expressDat = myDataERCC[-c(1)] 
  sampleNameList = c(sample1,sample2)
  libenum = c("1","2","3","4") # Library 5 is not included
  #if(siteName == "NWU")libenum = c("1","2","3")
  expressDatAandBERCCs = merge(idColsAdj[c(1,4)],myDataERCC)
  expressDatAandBERCCs$Feature <- as.factor(as.character(expressDatAandBERCCs$Feature))
  expressDat = expressDatAandBERCCs[-c(1:2)]
  expressDatSampleAB = expressDatAandBERCCs[c(1:2)]
  
  # Create long data frame of counts for each library in a for loop
  for (libe in 1:length(libenum)){
    if(libe == 1){
      select = subset(designMatSum, (Library == libenum[libe]))
      select <- as.data.frame(lapply(select,as.character))
      select <- as.data.frame(lapply(select,as.factor))
      expressDatForSum = expressDat[c(match(select$countSet, names(expressDat)))]
      expressDatSampleAdd = cbind(expressDatSampleAB, expressDatForSum)
      expressDatSampleAdd$Library = libe
      names(expressDatSampleAdd)[3:4]= c(sample1,sample2)
    }else{
      select = subset(designMatSum, (Library == libenum[libe]))
      select <- as.data.frame(lapply(select,as.character))
      select <- as.data.frame(lapply(select,as.factor))
      expressDatForSum = expressDat[c(match(select$countSet, names(expressDat)))]
      #meanCounts = data.frame(rowMeans(expressDatForSum))
      addRows = cbind(expressDatSampleAB[c(1:2)],expressDatForSum)
      addRows$Library = libe
      #print(addRows)
      names(addRows)[3:4]= c(sample1,sample2)
      expressDatSampleAdd = rbind(expressDatSampleAdd,addRows )  
    } 
  }
  expressDatSampleAdd$Library = as.factor(expressDatSampleAdd$Library)
  expressDatSampleAdd$Feature = as.factor(expressDatSampleAdd$Feature)
  
###############################################################################################################

  # Test for differential expression if DEtest == T
  # else use existing pvalues of the data in csv files
  if (DEtest == T){
    suppressWarnings(testDE(filenameRoot = filenameRoot, cnt = sampleLibeDataNoNorm, info = designMatAB, idCols = idColsAdj, FCcode = FCcode, r_m.mn = r_m.mn, totalSeqReads = totalSeqReads, totalReads = totalReads,legendLabels=legendLabels))
  } 
  
  deRes <- read.csv(file = paste(filenameRoot,"quasiSeq.res.csv",sep="."))
  p.thresh<-.1
  if(any(deRes$qvals<choseFDR)) p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
  print("Threshold P-value")
  print(p.thresh)
#   if (p.thresh > .1){
#     print(paste("Threshold P-value is high for the chosen FDR of ", as.character(choseFDR)))
#     print("The sample comparison indicates a large amount of differential expression in the measured transcript populations")
#     print("Resetting Threshold to 0.1")
#     p.thresh<- 0.1
#   }
###############################################################################################################

  # Find LODR estimates using the ERCC data pvalues
  lodr.ERCC = findLODR(pval.cutoff=p.thresh, prob=.9,filenameRoot=filenameRoot, kind="ERCC", FCcode = FCcode, legendLabels=legendLabels)
  
  # Find LODR estimates using data pvalues simulated from endogenous transcripts
  lodr.Sim = findLODR(pval.cutoff=p.thresh, prob=.9,filenameRoot=filenameRoot,kind="Sim", FCcode = FCcode, legendLabels=legendLabels)  
  
  print(lodr.ERCC)
  lodr.ERCC = data.frame(lodr.ERCC)
  Fold = lodr.ERCC[c(1)]
  Count = lodr.ERCC$Estimate
  
  print(lodr.ERCC$Estimate)
  
  #ConcEst = (log2(lodr.ERCC$Estimate))/fit.coeff[2,1]
  print("Log2 Counts (read depth normalized) LODR")
  # Convert LODR count estimate to library size normalized and 0 count adjusted (+1) count value
  logCount = log2((lodr.ERCC$Estimate/(mnLibeFactor))+1)
  print(logCount)
  
  print(fit.coeff[1])
  print(fit.coeff[2])
  
  print("Log2 Minimum ERCC Design Concentration:")
  ConcEst = (log2(lodr.ERCC$Estimate/(mnLibeFactor))-fit.coeff[1])/fit.coeff[2]
  print(ConcEst)
  
  print("Minimum ERCC Design Concentration:")
  print(2^(ConcEst))
  
  LODR.print.res = data.frame(Fold, Count, logCount, ConcEst)
  
  names(LODR.print.res)<- c("Fold","Count","Log2Count_normalized","Log2Conc")
  print(LODR.print.res)

###############################################################################################################

  ### Generate MA plots of erccs coded by concentrations from LODR
    
  maPlotAB = maConcPlot(idCols=idColsAdj,siteName = siteName, analysis = analysis, countPair = expressDatSampleAdd, r_m.mn = r_m.mnlog, sample1 = sample1,sample2 = sample2, myYLim = myYLimMA, myXLim=myXLim, replicate = T, cutoffs = ConcEst[which(!(is.na(ConcEst)))], FCcode = FCcode, legendLabels=legendLabels)

  print("Plotting MA plot with coding")
  print(maPlotAB[[1]])
  sdGlobal = maPlotAB$sdGlobal  
###############################################################################################################
# Generate ROC curves for the selected differential ratios
  erccROC.res = erccROC(filenameRoot = filenameRoot, kind = "ERCC", folds = FCcode, legendLabels = legendLabels)
  print(erccROC.res$ROCplot)
  AUCdat = erccROC.res$AUCdat

###############################################################################################################
# name and consolidate metrics for the interlaboratory comparison
  nam <- paste(filenameRoot, "AUC",sep = ".")
  assign(nam,AUCdat)

  nam <- paste(filenameRoot, "lodr.ERCC",sep = ".")
  assign(nam,lodr.ERCC)

  nam <- paste(filenameRoot, "sdGlobal",sep = ".")
  assign(nam,sdGlobal)

  nam <- paste(filenameRoot, "r_m",sep = ".")
  assign(nam,r_m.mn)

  nam <- paste(filenameRoot, "r_m_lower",sep = ".")
  assign(nam,r_m.lower)

nam <- paste(filenameRoot, "r_m_upper",sep = ".")
assign(nam,r_m.upper)

nam <- paste(filenameRoot, "p.thresh",sep = ".")
assign(nam,p.thresh)

###############################################################################################################
if (printPDF == T){
  dev.off()  
}

to.save <- ls()

save(list = to.save[grepl(pattern = filenameRoot,x=to.save)],file=paste(filenameRoot,"RData",sep = "."))

#rm(list = ls())

###########  END OF SCRIPT  ###############