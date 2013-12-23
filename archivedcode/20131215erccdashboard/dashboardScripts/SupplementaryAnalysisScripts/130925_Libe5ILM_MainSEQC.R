##############################################################################
# Analysis of per site results from SEQC from a given analysis pipeline
##############################################################################
sampleInfo = list(study = "SEQC_Main", #
                  siteName = "MAY",
                  platform = "ILM",
                  analysis = "NCTR",
                  printPDF = F,
                  DEtest = T,
                  totalSeqReads = T,
                  libeSizeNorm = T,
                  sample1Name = "A",
                  sample2Name = "B",
                  FCcode = data.frame(Ratio = c("a","b","c","d"),
                                      FC =  c(4,1,.667,.5)),
                  myYLimMA = c(-3.5,3.5),
                  myXLim = c(-10,15),
                  myYLim = NULL,
                  xlimEffects = NULL,
                  choseFDR = 0.01,
                  legendLabels = c("4:1","1:1","1:1.5","1:2"),
                  ERCCdilution = 1,
                  spikeVol = 50,
                  totalRNAmass = 2.5*10^(3)
)


###############################################################################
# Input parameter adjustments
# For platform = "ILM" use analysis = "NCTR", siteName may be equal to "MAY","BGI","CNL","NVS","AGR",or "COH"
# For platform = "LIF" use analysis = "NIST", siteName may be equal to "SQW","NWU",or "PSU"
# ERCC dilution corresponds to dilution of the Ambion ERCC mixes prior to spiking, if no dilution was performed then ERÇCdilution = 1
# spikeVol = the amount of diluted ERCC mix spiked into the total RNA sample, units are in microliters
# totalRNAmass = mass of total RNA spiked with the diluted ERCC mix, units are in micrograms

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
library(locfit)
library(QuasiSeq)
library(grid)

# Source package scripts
source("R/combineTechReps.R") # combine technical replicates
source("R/compareReps.R") # compare technical replicates for a sample 
source("R/getDesignMat.R") # build factor matrix from sample names
source("R/dynRangePlot.R") # dose response curve for dynamic range
source("R/dynRangePlotLODR.R") # plot dyn Range with LODR
source("R/est_LODR.R") # estimate LODR
source("R/est_r_m.R") # estimate r_m
source("R/loadERCCInfo.R") # load ERCC definition files
#source("R/loadExpMeas.R") # load expression measures from file
source("R/loadExpMeasEdit.R") # load expression measures from file
source("R/maConcPlot.R") # MA plot coded with LODR
source("R/multiplot.R") # plot multiple plots on the same page
source("R/printLODRres.R") # print LODR results table
source("R/QL.results.R")# load the updated QL.results function for QuasiSeq 
source("R/rocPlot.R") # ROC Plot
source("R/testDE.R") # test for DE, obtain p and q values
source("R/dashboardPDF.R") #
source("R/prepERCCDat.R") #
source("R/libeSizeNorm.R")
source("R/estMnLibeFactor.R")
source("R/meltedExpDatReshape.R")
source("R/geneExprTest.R")
source("R/printLODRres.R")
source("R/saveResults.R")

###############################################################################
# Set up PDF file to print
sampleInfo <- dashboardPDF(sampleInfo)

###############################################################################

# Run loadERCCInfo function to obtain ERCC information
erccInfo <- loadERCCInfo(erccMixes="Ambion4plex")

###############################################################################
# Run loadExpMeas function to process count data file
expDat <- loadExpMeas(sampleInfo,erccInfo)

###############################################################################
# library size normalize the data
expDat <- libeSizeNorm(sampleInfo, expDat)

###############################################################################
# length normalize the ERCC concentrations
expDat <- prepERCCDat(sampleInfo, expDat)

###############################################################################


# Revert to raw counts to sum fluidic replicates and then sum by library replicate
# exclude library 5

#sampleLibeData = combineTechReps(type = "sum", sampleInfo, expDat, designMat = expDat$designMatAB, sampleNameList = c(expDat$sample1,expDat$sample2), libeList = levels(expDat$designMatAB$Library)[-5])
#processReps<- function(){
  expDat = combineTechReps(type = "sum", sampleInfo, expressDat = expDat$TranscriptsAB, designMat = expDat$designMatAB, sampleNameList = c(expDat$sample1,expDat$sample2), libeList = levels(expDat$designMatAB$Library))
  
  sampleLibeDataNorm = expDat$expressDatCombinedSum
  
  sampleLibeDataNoNorm = expDat$allSampleLibeData
  
  sampleLibeSums = colSums(sampleLibeDataNoNorm[-c(1)],na.rm =T)
  
  if(sampleInfo$totalSeqReads==T){
    sampleLibeSums = expDat$totalReadSum
    expDat$totalReads = sampleLibeSums
  } 
  
  designMatNorm = getDesignMat(sampleLibeDataNorm, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')
  
Libe5UHRR <- sampleLibeDataNoNorm[,c(1,6)]
Libe5HBRR <- sampleLibeDataNoNorm[,c(1,11)]
sampleInfo$siteName <- "BGI"

###############################################################################
# Run loadExpMeas function to process count data file
expDat <- loadExpMeas(sampleInfo,erccInfo)

###############################################################################
# library size normalize the data
expDat <- libeSizeNorm(sampleInfo, expDat)

###############################################################################
# length normalize the ERCC concentrations
expDat <- prepERCCDat(sampleInfo, expDat)

###############################################################################


# Revert to raw counts to sum fluidic replicates and then sum by library replicate
# exclude library 5

#sampleLibeData = combineTechReps(type = "sum", sampleInfo, expDat, designMat = expDat$designMatAB, sampleNameList = c(expDat$sample1,expDat$sample2), libeList = levels(expDat$designMatAB$Library)[-5])
#processReps<- function(){
expDat = combineTechReps(type = "sum", sampleInfo, expressDat = expDat$TranscriptsAB, designMat = expDat$designMatAB, sampleNameList = c(expDat$sample1,expDat$sample2), libeList = levels(expDat$designMatAB$Library))

sampleLibeDataNorm = expDat$expressDatCombinedSum

sampleLibeDataNoNorm = expDat$allSampleLibeData

sampleLibeSums = colSums(sampleLibeDataNoNorm[-c(1)],na.rm =T)

if(sampleInfo$totalSeqReads==T){
  sampleLibeSums = expDat$totalReadSum
  expDat$totalReads = sampleLibeSums
} 

designMatNorm = getDesignMat(sampleLibeDataNorm, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')

Libe5UHRR <- merge(Libe5UHRR, sampleLibeDataNoNorm[,c(1,6)],by="Feature")
Libe5HBRR <- merge(Libe5HBRR, sampleLibeDataNoNorm[,c(1,11)],by="Feature")


sampleInfo$siteName <- "CNL"

###############################################################################
# Run loadExpMeas function to process count data file
expDat <- loadExpMeas(sampleInfo,erccInfo)

###############################################################################
# library size normalize the data
expDat <- libeSizeNorm(sampleInfo, expDat)

###############################################################################
# length normalize the ERCC concentrations
expDat <- prepERCCDat(sampleInfo, expDat)

###############################################################################


# Revert to raw counts to sum fluidic replicates and then sum by library replicate
# exclude library 5

#sampleLibeData = combineTechReps(type = "sum", sampleInfo, expDat, designMat = expDat$designMatAB, sampleNameList = c(expDat$sample1,expDat$sample2), libeList = levels(expDat$designMatAB$Library)[-5])
#processReps<- function(){
expDat = combineTechReps(type = "sum", sampleInfo, expressDat = expDat$TranscriptsAB, designMat = expDat$designMatAB, sampleNameList = c(expDat$sample1,expDat$sample2), libeList = levels(expDat$designMatAB$Library))

sampleLibeDataNorm = expDat$expressDatCombinedSum

sampleLibeDataNoNorm = expDat$allSampleLibeData

sampleLibeSums = colSums(sampleLibeDataNoNorm[-c(1)],na.rm =T)

if(sampleInfo$totalSeqReads==T){
  sampleLibeSums = expDat$totalReadSum
  expDat$totalReads = sampleLibeSums
} 

designMatNorm = getDesignMat(sampleLibeDataNorm, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')

Libe5UHRR <- merge(Libe5UHRR, sampleLibeDataNoNorm[,c(1,6)],by="Feature")
Libe5HBRR <- merge(Libe5HBRR, sampleLibeDataNoNorm[,c(1,11)],by="Feature")

#######################

libe5ILMDat <- cbind(Libe5UHRR,Libe5HBRR[-c(1)]) 

sampleInfo$totalSeqReads = F

###############################################################################
# Dynamic Range Plots of Sample Biological Replicates
###############################################################################

### Estimate r_m for the sample pair using a negative binomial glm
expDat <- est_r_m(sampleInfo, expDat, cnt = libe5ILMDat, printPlot = T)


dynRangeDat = dynRangePlot(sampleInfo, expDat, expressDat = libe5ILMDat,
                           designMat = designMatNorm)

###############################################################################

expDat <- estMnLibeFactor(sampleInfo, expDat, cnt = libe5ILMDat)

expDat <- meltedExpDat(sampleInfo, expDat, cnt = libe5ILMDat, designMat = designMatNorm)

# # Test for differential expression if DEtest == T
# # else use existing pvalues of the data in csv files
# expDat <- geneExprTest(sampleInfo, expDat, cnt = sampleLibeDataNoNorm, info = designMatNorm )
# 
# # Find LODR estimates using the ERCC data pvalues
# lodr.ERCC = findLODR(sampleInfo,expDat,kind = "ERCC", prob=0.9)
# 
# # Find LODR estimates using data pvalues simulated from endogenous transcripts
# lodr.Sim = findLODR(sampleInfo, expDat, kind = "Sim", prob = 0.9)  
# 
# #### Edit starting here
# LODR.print.res.ERCC <- printLODRres(sampleInfo, expDat, dynRangeDat,
#                                     lodr.res = lodr.ERCC)
# 
# ###############################################################################
# # Annotate the dynamic range plots to show how LODR is used to estimate spike 
# #  in concentration
# # Move legend
# dynRangePlotLODR(dynRangeRes = dynRangeDat$dynRangePlotRes,
#                  lodrDat = LODR.print.res.ERCC$LODRtable, folds = sampleInfo$FCcode, 
#                  legendLabels = sampleInfo$legendLabels)
# 
# ###############################################################################

# Generate MA plots of erccs coded by concentrations from LODR

maPlotAB = maConcPlot(sampleInfo, expDat, countPair = expDat$expressDat_l , cutoffs = NULL, alphaPoint = 0.8,r_mAdjust = T, replicate = T,ReplicateName="Site")

print("Plotting MA plot with coding")
print(maPlotAB[[1]])
ILMexpDat <- expDat
ILMmaPlotAB <- maPlotAB
save(list = c("ILMexpDat","ILMmaPlotAB"),file="ILMLibe5Res.RData")
# ###############################################################################
# # Generate ROC curves for the selected differential ratios
# erccROC.res = erccROC(filenameRoot = sampleInfo$filenameRoot, kind = "ERCC",
#                       legendLabels = sampleInfo$legendLabels,
#                       folds = sampleInfo$FCcode,idCols = expDat$idColsAdj)
# print(erccROC.res$ROCplot)
# #AUCdat = erccROC.res$AUCdat
# 
# ###############################################################################
# saveResults(sampleInfo, erccROC.res, maPlotAB, lodr.ERCC, expDat)

###############################################################################
### End of Script