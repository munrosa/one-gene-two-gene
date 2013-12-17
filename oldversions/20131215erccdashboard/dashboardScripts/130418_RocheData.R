# Analysis of Roche Data
# Duplicate measurements made at 3 different SEQC sites
# Treat these as 6 sample replicates each for UHRR and HBRR 

### Import data,use same filter criteria as for other RNA-Seq samples

### Run this commented code line here and then source the rest of the script
# study = "SEQC_Main";siteName = "ALLROCHE";platform = "ROC";analysis = "WEHI";libeSizeNorm = T;printPDF = T; DEtest = T; totalSeqReads = T; sample1 = "A";sample2 = "B";  FCcode = data.frame(Ratio = c("a","b","c","d"),FC =  c(4,1,.667,.5)); myYLimMA = c(-3.5,3.5); myXLim = c(-10,15); choseFDR = 0.01; legendLabels = c("4:1","1:1","1:1.5","1:2"); ERCCdilution = 1; spikeVol = 50; totalRNAmass = 2.5*10^(3)

#####
# Input parameter adjustments
# For platform = "ILM" use analysis = "NCTR", siteName may be equal to "MAY","BGI","CNL","NVS","AGR",or "COH"
# For platform = "LIF" use analysis = "NIST", siteName may be equal to "SQW","NWU",or "PSU"
# For platform = "ROC" use analysis = "WEHI"

# ERCC dilution corresponds to dilution of the Ambion ERCC mixes prior to spiking, if no dilution was performed then ERÇCdilution = 1
# spikeVol = the amount of diluted ERCC mix spiked into the total RNA sample, units are in microliters
# totalRNAmass = mass of total RNA spiked with the diluted ERCC mix, units are in micrograms

###############################################################################
# Set all warnings to report as they appear
options(warn = 1)

# Get required libraries
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(edgeR)
library(locfit)
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
source("R/printLODRres.R") # print LODR results table
source("R/dynRangePlotLODR.R") # plot dyn Range with LODR
source("R/QL.results.R")# load the updated QL.results function for QuasiSeq 

###############################################################################################################
filenameRoot = paste(siteName,analysis,platform,sample1,sample2,sep = ".")
theme_set(theme_bw(base_size=12))
# Open PDF file to write results
if (printPDF == T){
  pdf(file = paste(filenameRoot,"pdf",sep="."),onefile=T,width=7,height = 7)  
}
###############################################################################################################

# Run loadERCCInfo function to obtain ERCC information
erccInfo <- loadERCCInfo()
idCols = erccInfo$idCols
MixDef = erccInfo$MixDef

###############################################################################################################
# Run loadExpMeas function to process count data file
expDat <- loadExpMeas(study = study, siteName = siteName, platform = platform, analysis = analysis, totalSeqReads = totalSeqReads, sample1 = sample1, sample2 = sample2, idCols = idCols)

TranscriptsAB <- expDat[[1]]
designMatAB <- expDat[[2]]
designMatAll <- expDat[[3]]
sample1 <- expDat[[4]]
sample2 <- expDat[[5]]
idCols <- expDat[[6]]
totalReads <- expDat[[7]]

# write Transcript csv file to directory
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
  print("Total Mapped Reads Library size summary:")
  print(summary(libeSize))
}
if ((libeSizeNorm == T)&(totalSeqReads == T)){
  TranscriptsAll = expressDat
  libeSize = totalReads
  datCols = expressDat[-c(1)]
  libeSize = libeSize/(10^6) #per million mapped reads
  #Library size normalize the data  
  libAdjust = sweep(datCols, 2, libeSize,"/")
  expressDat = cbind(expressDat[c(1)], libAdjust)
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

### Calculate the per ERCC amount spiked attomoles of nt / ug total RNA
### ERCCdilution = 1; spikeVol = 50; totalRNAmass = 2.5*10^(3)
spikeFraction <- (ERCCdilution*spikeVol)/(totalRNAmass)
idCols$Conc1 <- idCols$Conc1*spikeFraction
idCols$Conc2 <- idCols$Conc2*spikeFraction

ERCCxlabelIndiv = expression(paste("Log2 ERCC Spike Amount(attomol nt/",mu,"g total RNA)",sep = ""))
ERCCxlabelAve = expression(paste("Log2 Average ERCC Spike Amount(attomol nt/",mu,"g total RNA)",sep = ""))

###############################################################################################################
### Estimate r_m for the sample pair using a negative binomial glm
r_m.res = glm_est_r_m(site = filenameRoot, cnt = TranscriptsAB, libeSizeNorm = libeSizeNorm,totalSeqReads = totalSeqReads,totalReads = totalReads, idCols = idCols, FCcode = FCcode,sample1 = sample1, sample2 = sample2,legendLabels=legendLabels,myXLim=myXLim,avexlabel = ERCCxlabelAve, printPlot = T)

#Return the mean estimate of r_m
r_m.mn = r_m.res$r_m.mn

#Return the upper bound of the 95% confidence interval for r_m
r_m.upper = r_m.res$r_m.upper

#Return the lower bound of the 95% confidence interval for r_m
r_m.lower = r_m.res$r_m.lower

# Return the r_m adjusted concentrations in idCols
idColsAdj = r_m.res$idCols

###############################################################################################################

### Evaluate dynamic range
print("Signal-Abundance Plots for dynamic range estimation...") 
dynRangeDat = dynRangePlot(expressDat = TranscriptsAB, designMat = designMatAB,idCols = idColsAdj, sampleNames = c(sample1,sample2),myXLim=myXLim,myYLim=NULL, legendLabels=legendLabels,FCcode=FCcode, indivxlabel = ERCCxlabelIndiv, avexlabel= ERCCxlabelAve)

print(dynRangeDat$fit.coeff)
fit.coeff = dynRangeDat$fit.coeff
dynRangeRes <- dynRangeDat[[4]]


###############################################################################################################

### Estimate average totalSeqReads to use for mnLibeFactor estimation
if(totalSeqReads == F){
  sampleLibeSums = colSums(TranscriptsAB[-c(1)],na.rm =T)
  mnLibeFactor = (mean(as.vector(sampleLibeSums)))/10^6
  print("Using total mapped reads mean library size = ")
  print(mnLibeFactor)
}else{
  mnLibeFactor = (mean(as.vector(totalReads)))/10^6
  print("Using total sequencing reads mean library size = ")
  print(mnLibeFactor)
  sampleLibeSums = totalReads
}
###############################################################################################################

#####
# After observing that Library 5 is an outlier it is excluded from further analysis
# The four remaining library replicates for each sample are collected into a data frame with the header names:
# Feature Ratio Sample1 Sample2 Library
#####
expressDat <- TranscriptsAB[-c(1)]

colnames(expressDat)<-paste(rep(c(sample1,sample2),each=ncol(expressDat)/2),c(1:(ncol(expressDat)/2),1:(ncol(expressDat)/2)),sep="_")

expressDat <- cbind(TranscriptsAB[c(1)], expressDat)

designMatSum <- generateDesignMat(expressionData=expressDat,factorList=c("Sample","Replicate"),patternSplit="_")

myDataERCC = expressDat
expressDat = myDataERCC[-c(1)] 
sampleNameList = c(sample1,sample2)
libenum = c("1","2","3","4","5","6") 
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

###############################################################################################################

# Test for differential expression if DEtest == T
# else use existing pvalues of the data in csv files
if (DEtest == TRUE){
  suppressWarnings(testDE(filenameRoot = filenameRoot, cnt = TranscriptsAB, info = designMatAB, idCols = idColsAdj, FCcode = FCcode, r_m.mn = r_m.mn, totalSeqReads = totalSeqReads, totalReads = totalReads, legendLabels = legendLabels))
} 
if (DEtest == FALSE){
  print("Not testing for differential expression, no dispersion plots will be produced and pre-existing DE test result files will be used for additional analysis")
}
deRes <- read.csv(file = paste(filenameRoot,"quasiSeq.res.csv",sep="."))
p.thresh<-.1

if(any(deRes$qvals<choseFDR)) p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
print("Threshold P-value")
print(p.thresh)
if (p.thresh > .1){
  print(paste("Threshold P-value is high for the chosen FDR of ", as.character(choseFDR)))
  print("The sample comparison indicates a large amount of differential expression in the measured transcript populations")
}
###############################################################################################################

# Find LODR estimates using the ERCC data pvalues
lodr.ERCC = findLODR(pval.cutoff=p.thresh, prob=.9,filenameRoot=filenameRoot, kind="ERCC", FCcode = FCcode, legendLabels = legendLabels)

# Find LODR estimates using data pvalues simulated from endogenous transcripts
lodr.Sim = findLODR(pval.cutoff=p.thresh, prob=.9,filenameRoot=filenameRoot, kind="Sim", FCcode = FCcode, legendLabels = legendLabels)  

LODR.print.res.ERCC <- printLODRres(lodr.res = lodr.ERCC, fit.coeff = fit.coeff, mnLibeFactor = mnLibeFactor,FCcode = FCcode, legendLabels = legendLabels)
print(LODR.print.res.ERCC[[1]])

###############################################################################################################
### Annotate the dynamic range plots to show how LODR is used to estimate spike in concentration
#Move legend
dynRangePlotLODR(dynRangeRes = dynRangeRes,lodrDat = LODR.print.res.ERCC[[1]], folds = FCcode, legendLabels = legendLabels)


###############################################################################################################

### Generate MA plots of erccs coded by concentrations from LODR

maPlotAB = maConcPlot(idCols=idColsAdj,siteName = siteName, analysis = analysis, countPair = expressDatSampleAdd, r_m.mn = exp(r_m.mn),r_m.UL = exp(r_m.upper),r_m.LL = exp(r_m.lower), sample1 = sample1,sample2 = sample2, myYLim = myYLimMA, myXLim=myXLim, replicate = T, cutoffs = LODR.print.res.ERCC[[2]], FCcode = FCcode, legendLabels=legendLabels, avexlabel = ERCCxlabelAve, spikeFraction = spikeFraction)

print("Plotting MA plot with coding")
print(maPlotAB[[1]])
sdGlobal = maPlotAB$sdGlobal
modRatVar = maPlotAB$stdevCoef
###############################################################################################################
# Generate ROC curves for the selected differential ratios
erccROC.res = erccROC(filenameRoot = filenameRoot, kind = "ERCC", legendLabels = legendLabels, idCols = idColsAdj)
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

nam <- paste(filenameRoot, "modRatVar",sep = ".")
assign(nam,modRatVar)

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

### R_m estimate

### dynamic Range

### MA plots

### ROC curves
