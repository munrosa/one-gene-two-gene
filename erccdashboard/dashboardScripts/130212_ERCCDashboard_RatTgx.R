##############################################################################
# Analysis of Rat Toxicogenomics data from SEQC
##############################################################################
# Run the following at the command line with necessary adjustments and then 
# source this script for the ERCC dashboard diagnostics

# siteName = "COH";platform = "ILM";analysis = "RatTox";libeSizeNorm = T;printPDF = T; DEtest = T; totalSeqReads = T; sample1 = "NIT";sample2 = "CTL";  FCcode = data.frame(Ratio = c("a","b","c","d"),FC =  c(4,1,.667,.5)); myYLimMA = c(-3,3); myXLim = c(-10,15); choseFDR = 0.05; deliveryMet = "NN_IP"; legendLabels = c("4:1","1:1","1:1.5","1:2"); ERCCdilution = 1/100; spikeVol = 1; totalRNAmass = 0.500;ERCCxlabel = expression(paste("Log2 ERCC Spike Amount(attomol nt/",mu,"g total RNA)",sep = ""))
#ERCCxlabel  = expression(paste("Log2 ERCC Spike Amount(attomol nt/",mu,"g)",sep = ""))

######
### Explanation of input parameters that can be adjusted for this script
# totalSeqReads - if TRUE, use the number of sequenced reads in fasta file for library size normalization, if FALSE use number of mapped reads for library size normalization
# sample1 - can be adjusted for "NIT","THI","3ME","MET","NAP"
# deliveryMet - can be "NN_IP" or "NN_OG" according to sample1 choice
# ERCC dilution corresponds to dilution of the Ambion ERCC mixes prior to spiking, if no dilution was performed then ERÇCdilution = 1
# spikeVol = the amount of diluted ERCC mix spiked into the total RNA sample, units are in microliters
# totalRNAmass = mass of total RNA spiked with the diluted ERCC mix, units are in micrograms

#####


###############################################################################
# For NN_IP samples
# CTL has ERCC Mix 2; Treatments NIT and THI have ERCC Mix 1
# Same instrument different days for library prep and sequencing for the CTL and Treatment samples

# For NN_OG samples
# CTL has ERCC Mix 2; Treatments 3ME,MET,NAP, have ERCC Mix 1
# Same instrument same library prep and sequencing dates for the CTL and Treatment samples
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

## Run loadExpMeas function to process data file
# Transcripts <- loadExpMeas(siteName = siteName, platform = platform, analysis = analysis)

# Just process data directly, skip loadExpMeas
Transcripts <- read.delim("Data/RatToxSEQC/SEQC_TGx_GeneCounts_JMEEHAN.txt")
# force names to be ERCC- and first column name to Feature
names(Transcripts)[1] = "Feature"
Transcripts$Feature = gsub("ERCC_","ERCC-",Transcripts$Feature)
Transcripts$Feature = gsub(":Gene_AceView08","",Transcripts$Feature)
Transcripts$Feature = gsub(":Gene_RefSeq","",Transcripts$Feature)

# get the total reads per sample (from sequence files prior to mapping)
ratToxReads <- read.csv("Data/RatToxSEQC/RatToxTotalReads54subset21.csv")

###############################################################################################################

# Run loadERCCInfo function to obtain idCols and MixDef information
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
if(analysis == "RatTox"){
  designMat = generateDesignMat(TranscriptsERCCOnly, factorList = c("Site","PI","Flowcell","Barcode","Tissue","Chemical","Vehicle","Route","SeqBarcode","Lane"), patternSplit = '_')
}
###############################################################################################################
filenameRoot = paste(siteName,analysis,platform,sample1,sample2,deliveryMet,sep = ".")
# Open PDF file to write results
if (printPDF == T){
  pdf(file = paste(filenameRoot,"pdf",sep="."),onefile=T,width=7,height = 7)  
}

###############################################################################################################

# Subset to get just the A and B samples for the Count Matrix and for the designMatrix
if (deliveryMet == "NN_IP"){
  select1 = subset(designMat, (Chemical == sample1)&(Vehicle == "NN")&(Route == "IP"))
  print(select1)
  select2 = subset(designMat, (Chemical == sample2)&(Flowcell == "AB029JACXX")&(Vehicle == "NN")&(Route == "IP")&((Lane == "s_4")|(Lane == "s_6")|(Lane == "s_1")))
  print(select2)
  select = rbind(select1,select2)  
}
if(deliveryMet == "NN_OG"){
  select1 = subset(designMat, (Chemical == sample1)&(Vehicle == "NN")&(Route == "OG"))
  print(select1)
  select2 = subset(designMat, (Chemical == sample2)&(Flowcell == "AB029JACXX")&(Vehicle == "NN")&(Route == "OG"))
  print(select2)
  select = rbind(select1,select2)
}

select <- as.data.frame(lapply(select,as.character))
select <- as.data.frame(lapply(select,as.factor))

designMatAB = select

dataAB = Transcripts[-c(1)]
TranscriptsAB = cbind(Transcripts[c(1)],dataAB[c(match(select$countSet, names(dataAB)))])

totalReads = ratToxReads$total_reads[c(match(select$countSet,ratToxReads$Alt_ID))]


dataAB <- TranscriptsAB[-c(1)]
colnames(dataAB)<-paste(rep(c(sample1,sample2),each=ncol(dataAB)/2),c(1:(ncol(dataAB)/2),1:(ncol(dataAB)/2)),sep="_")
print(colnames(dataAB))

TranscriptsAB <- cbind(Feature = TranscriptsAB$Feature, dataAB)
#filter data starting with sample 1
sample1Cols = grep(pattern = sample1, colnames(TranscriptsAB))
idxsample1 <- which((rowMeans(TranscriptsAB[sample1Cols])>1)&(rowSums(TranscriptsAB[sample1Cols]>(length(sample1Cols)))))

sample2Cols = grep(pattern = sample2, colnames(TranscriptsAB))
idxsample2 <- which((rowMeans(TranscriptsAB[sample2Cols])>1)&(rowSums(TranscriptsAB[sample2Cols]>(length(sample2Cols)-1))))



#TranscriptsAB<-subset(TranscriptsAB, (rowMeans(TranscriptsAB[-c(1)])>1)&(rowSums(TranscriptsAB[-c(1)]!=0)>2));


write.csv(TranscriptsAB, paste(filenameRoot,"Transcripts.csv",sep="."),row.names = F)
###############################################################################################################

# Library size normalize the data if libeSizeNorm = T
expressDat = TranscriptsAB

if ((libeSizeNorm == T)&(totalSeqReads == F)){
  TranscriptsAll = expressDat[-c(grep("ERCC-0", expressDat$Feature)),]  
  TranscriptMappedReadSums = colSums(TranscriptsAll[-c(1)],na.rm = T)
  libeSize = TranscriptMappedReadSums
  datCols = expressDat[-c(1)]
  libeSize = libeSize/(10^6) #per million mapped reads
  #Library size normalize the data  
  libAdjust = sweep(datCols, 2, libeSize,"/")
  expressDat = cbind(expressDat[c(1)], libAdjust)
  }
if ((libeSizeNorm == T)&(totalSeqReads == T)){
  TranscriptsAll = expressDat[-c(grep("ERCC-0", expressDat$Feature)),]  
  libeSize = totalReads
  datCols = expressDat[-c(1)]
  libeSize = libeSize/(10^6) #per million mapped reads
  #Library size normalize the data  
  libAdjust = sweep(datCols, 2, libeSize,"/")
  expressDat = cbind(expressDat[c(1)], libAdjust)
}

print("Library sizes:")
print(libeSize)

# get just ERCC data in the expression data frame
expressDat = expressDat[c(grep("ERCC-0", expressDat$Feature)),]

  #Length normalize the Expected ERCC concentrations
  lengthFactor = (idCols$Length)/(1000)
  
  # If length normalization of the expected concentrations is desired (default)
  idCols$Conc1 = (idCols$Conc1*lengthFactor)
  idCols$Conc2 = (idCols$Conc2*lengthFactor)  
  
### Calculate the per ERCC amount spiked attomoles of nt/ ng total RNA
### ERCCdilution = 1/100; spikeVol = 1; totalRNAmass = 500

#idCols$Conc1 <- idCols$Conc1*(ERCCdilution*spikeVol)/(totalRNAmass)
#idCols$Conc2 <- idCols$Conc2*(ERCCdilution*spikeVol)/(totalRNAmass)

idCols$Conc1 <- idCols$Conc1*(ERCCdilution*spikeVol)#/(totalRNAmass)
idCols$Conc2 <- idCols$Conc2*(ERCCdilution*spikeVol)#/(totalRNAmass)


###############################################################################################################
# Dynamic Range Plots of Sample Biological Replicates
###############################################################################################################

print("Check for sample mRNA fraction differences...")


  ### Estimate r_m for the sample pair using a negative binomial glm
  r_m.res = glm_est_r_m(site = filenameRoot, cnt = TranscriptsAB, libeSizeNorm = libeSizeNorm,totalSeqReads = totalSeqReads,totalReads = totalReads, idCols = idCols, FCcode = FCcode, sample1 = sample1,sample2=sample2, legendLabels = legendLabels,myXLim=myXLim, ERCCxlabel=ERCCxlabel)
  
  #Return the mean estimate of r_m
  r_m.mn = r_m.res$r_m.mn

  #Return the upper bound of the 95% confidence interval for r_m
  r_m.upper = r_m.res$r_m.upper

  #Return the lower bound of the 95% confidence interval for r_m
  r_m.lower = r_m.res$r_m.lower

  # Return the r_m adjusted concentrations in idCols
  idColsAdj = r_m.res$idCols
  
  #expressDatAandBERCCs = merge(idCols[c(1,4)],expressDatSampleAB)
  r_m.mnlog = exp(r_m.mn)
  
  designMatSum <- generateDesignMat(TranscriptsAB, factorList = c("Sample","Replicate"), patternSplit = '_')

print("Signal-Abundance Plots for dynamic range estimation...")

  dynRangeDat = dynRangePlot(expressDat = expressDat, designMat = designMatSum,idCols = idColsAdj, sampleNames = c(sample1,sample2),FCcode = FCcode, legendLabels = legendLabels,myXLim = myXLim, myYLim = NULL, ERCCxlabel = ERCCxlabel)


  print(dynRangeDat$fit.coeff)
  fit.coeff = dynRangeDat$fit.coeff
  
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
  
  
  libeSize <- sampleLibeSums
  datCols = TranscriptsAB[-c(1)]
  libAdjust = sweep(datCols, 2, libeSize,"/")
  sampleLibeDataNorm = cbind(TranscriptsAB[c(1)],libAdjust)
  myDataERCC = sampleLibeDataNorm
  expressDat = myDataERCC[-c(1)] 
  sampleNameList = c(sample1,sample2)
  libenum = c("1","2","3") # set for n = 3, will need to make flexible if n > 3
  expressDatAandBERCCs = merge(idColsAdj[c(1,4)],myDataERCC)
  expressDatAandBERCCs$Feature <- as.factor(as.character(expressDatAandBERCCs$Feature))
  expressDat = expressDatAandBERCCs[-c(1:2)]
  expressDatSampleAB = expressDatAandBERCCs[c(1:2)]
  
  # Create long data frame of counts for each library in a for loop
  for (libe in 1:length(libenum)){
    if(libe == 1){
      select = subset(designMatSum, (Replicate == libenum[libe]))
      select <- as.data.frame(lapply(select,as.character))
      select <- as.data.frame(lapply(select,as.factor))
      expressDatForSum = expressDat[c(match(select$countSet, names(expressDat)))]
      expressDatSampleAdd = cbind(expressDatSampleAB, expressDatForSum)
      expressDatSampleAdd$Library = libe
      names(expressDatSampleAdd)[3:4]= c(sample1,sample2)
    }else{
      select = subset(designMatSum, (Replicate == libenum[libe]))
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
  
  # Test for differential expression if DEtest == T
  # else use existing pvalues of the data in csv files
  if (DEtest == T){
    testDE(filenameRoot = filenameRoot, cnt = TranscriptsAB, info = designMatAB, idCols = idColsAdj, FCcode = FCcode, r_m.mn = r_m.mn,totalSeqReads = totalSeqReads, totalReads = totalReads, legendLabels = legendLabels)
  } 
deRes <- read.csv(file = paste(filenameRoot,"quasiSeq.res.csv",sep="."))
p.thresh<-.1

if(any(deRes$qvals<choseFDR)) p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
print("Threshold P-value")
print(p.thresh)
if (p.thresh > .1){
  print(paste("Threshold P-value is high for the chosen FDR of ", as.character(choseFDR)))
  print("The sample comparison indicates a large amount of differential expression in the measured transcript populations")
  #print("Resetting Threshold to 0.1")
  #p.thresh<- 0.1
}

  # Find LODR estimates using the ERCC data pvalues
  lodr.ERCC = findLODR(pval.cutoff=p.thresh, prob=.9,filenameRoot=filenameRoot, kind="ERCC", FCcode = FCcode, legendLabels = legendLabels)
  
  # Find LODR estimates using data pvalues simulated from endogenous transcripts
  lodr.Sim = findLODR(pval.cutoff=p.thresh, prob=.9,filenameRoot=filenameRoot, kind="Sim", FCcode = FCcode, legendLabels = legendLabels)  
  
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
  
  ### Generate MA plots of erccs coded by concentrations from LODR
  
  maPlotABnoLODR = maConcPlot(idCols = idColsAdj,siteName = siteName, analysis = analysis, countPair = expressDatSampleAdd, r_m.mn = r_m.mnlog, sample1 = sample1,sample2 = sample2, myYLim = myYLimMA,myXLim = myXLim, replicate = T, cutoffs = NULL,FCcode = FCcode, legendLabels = legendLabels, ERCCxlabel = ERCCxlabel)
  
  maPlotAB = maConcPlot(idCols=idColsAdj,siteName = siteName, analysis = analysis, countPair = expressDatSampleAdd, r_m.mn = r_m.mnlog, sample1 = sample1,sample2 = sample2, myYLim = myYLimMA, myXLim = myXLim, replicate = T, cutoffs = ConcEst[which(!(is.na(ConcEst)))], FCcode = FCcode, legendLabels = legendLabels, ERCCxlabel = ERCCxlabel)

  print("Plotting MA plot with coding")
  print(maPlotAB[[1]])
  sdGlobal = maPlotAB$sdGlobal

  #print("Plotting MA plot without coding")
  #print(maPlotABnoLODR[[1]])
  
  erccROC.res = erccROC(filenameRoot = filenameRoot, kind = "ERCC", legendLabels = legendLabels)
  print(erccROC.res$ROCplot)
  AUCdat = erccROC.res$AUCdat
  
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
  
if (printPDF == T){
  dev.off()  
}

to.save <- ls()

save(list = to.save[grepl(pattern = filenameRoot,x=to.save)],file=paste(filenameRoot,"RData",sep = "."))

###########  END OF SCRIPT  ###############

