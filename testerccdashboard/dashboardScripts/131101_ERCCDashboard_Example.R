


##############################################################################
# Analysis of per site results from SEQC from a given analysis pipeline
##############################################################################
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
source("R/initDat.R") # combine technical replicates
source("R/sumTechReps.R") # combine technical replicates
source("R/getDesignMat.R") # build factor matrix from sample names
source("R/dynRangePlot.R") # dose response curve for dynamic range
source("R/dynRangePlotLODR.R") # plot dyn Range with LODR
source("R/estLODR.R") # estimate LODR
source("R/est_r_m.R") # estimate r_m
source("R/loadERCCInfo.R") # load ERCC definition files
source("R/loadExpMeas.R") # load expression measures from file
source("R/maConcPlot.R") # MA plot coded with LODR
source("R/multiplot.R") # plot multiple plots on the same page
source("R/printLODRres.R") # print LODR results table
source("R/QL.results.R")# load the updated QL.results function for QuasiSeq 
source("R/erccROC.R") # ROC Plot
source("R/testDE.R") # test for DE, obtain p and q values
source("R/dashboardPDF.R") #
source("R/prepERCCDat.R") #
source("R/libeSizeNorm.R")
source("R/estMnLibeFactor.R")
source("R/meltedExpDat.R")
source("R/geneExprTest.R")
source("R/printLODRres.R")
source("R/saveResults.R")


###############################################################################
# Initialize expDat list
###############################################################################
# Input parameter adjustments
# expDat <- initDat(designFactors = c("Study","Platform","Site","Sample",
#                                     "Library"),siteName = "Lab5", 
#                   platform = "LIF", sample1Name = "UHRR", sample2Name = "HBRR",
#                   choseFDR = 0.01, ERCCdilution = 1, spikeVol = 50, 
#                   totalRNAmass = 2.5*10^(3),  
#                   printPDF = T, DEtest = T,
#                   totalSeqReads = T, libeSizeNorm = T, myYLimMA = c(-3.5,3.5),
#                   myXLim = c(-10,15),myYLim = NULL, filenameRoot = NULL)
# siteName choices are "Lab1" to "Lab6" platform = "ILM", analysis = "NCTR" 
# siteName choices are "Lab7" to "Lab9" platform = "LIF", analysis = "NIST"

# expDat <- initDat(study = "SEQC_RatTox", siteName = "COH", platform = "ILM",
#                   analysis = "RatTox", sample1Name = "NIT", sample2Name = "CTL",
#                   choseFDR = 0.1, ERCCdilution = 1/100, spikeVol = 1, 
#                   totalRNAmass = 0.500, filenameRoot = "myFilenameIsGood", 
#                   printPDF = F, DEtest = T,
#                   totalSeqReads = T, libeSizeNorm = T, myYLimMA = c(-3.5,3.5),
#                   myXLim = c(-10,15),myYLim = NULL)
# sample1Name choices are "NAP","THI","NIT","MET", & "3ME"

# ERCC dilution = the dilution of the Ambion ERCC mixes prior to spiking, if no 
# dilution was performed then ERCCdilution = 1
# spikeVol = the amount of diluted ERCC mix spiked into the total RNA sample, 
# units are in microliters
# totalRNAmass = mass of total RNA spiked with the diluted ERCC mix, 
# units are in micrograms
###############################################################################
# expDat <- initDat(designFactors = c("Study","Platform","Site","Sample",
#                                     "Library"),siteName = "Lab5", 
#                   platform = "LIF", sample1Name = "UHRR", sample2Name="HBRR",
#                   choseFDR = 0.01, ERCCdilution = 1, spikeVol = 50, 
#                   totalRNAmass = 2.5*10^(3),  
#                   printPDF = T, DEtest = T,
#                   totalSeqReads = T, libeSizeNorm = T, myYLimMA = c(-3.5,3.5),
#                   myXLim = c(-10,15),myYLim = NULL, filenameRoot = NULL)


expDat <- initDat(designFactors = c("Sample",
                                    "Library"),siteName = "COH", 
                  platform = "ILM", sample1Name = "THI", sample2Name = "CTL",
                  choseFDR = 0.1, erccMixes = "Ambion4plexPair",
                  ERCCdilution = 1/100, spikeVol = 1, totalRNAmass = 0.500,  
                  printPDF = F, DEtest = T, totalSeqReads = T, libeSizeNorm = T,
                  myYLimMA = c(-3.5,3.5), myXLim = c(-10,15),myYLim = NULL,
                  filenameRoot = NULL)


###############################################################################
# Run loadERCCInfo function to obtain ERCC information
expDat <- loadERCCInfo(expDat, erccMixes="Ambion4plexPair")

###############################################################################
# # Run loadExpMeas function to process count data file
# expDat <- loadExpMeasFromRaw(expDat)

# Assume user has created data frame countTable and totalReads vector
# process those data files to add to expDat structure
expDat <- loadExpMeas(expDat, countTable = COH.RatTox.ILM.THI.CTL.countTable,
                      designFactors = expDat$sampleInfo$designFactors, 
                      totalReads = COH.RatTox.ILM.THI.CTL.totalReads)

###############################################################################
# library size normalize the data
expDat <- libeSizeNorm(expDat)

###############################################################################
# length normalize the ERCC concentrations
expDat <- prepERCCDat(expDat)

###############################################################################
### Estimate r_m for the sample pair using a negative binomial glm
#expDat <- est_r_m(expDat, cnt = expDat$expressDatSumNoNorm, printPlot = T)
expDat <- est_r_m(expDat, cnt = expDat$Transcripts, printPlot = T)
# ###############################################################################

### Evaluate dynamic range
#print("Signal-Abundance Plots for dynamic range estimation...")

dynRangeDat = dynRangePlot(expDat, expressDat = expDat$expressDat,
                           designMat = expDat$designMat, noErrorBars = F)

###############################################################################
# Estimate the mean library size factor for estimating Spike-in Concentration
# corresponding to LODR

expDat <- estMnLibeFactor(expDat, cnt = expDat$Transcripts)

expDat <- meltedExpDat(expDat, cnt = expDat$Transcripts, 
                       designMat = expDat$designMat)

# Test for differential expression if DEtest == T
# else use existing pvalues of the data in csv files
expDat <- geneExprTest(expDat, cnt = expDat$Transcripts,
                       designMat = expDat$designMat )

# Find LODR estimates using the ERCC data pvalues
lodr.ERCC = estLODR(expDat,kind = "ERCC", prob=0.9)

# Find LODR estimates using data pvalues simulated from endogenous transcripts
lodr.Sim = estLODR(expDat, kind = "Sim", prob = 0.9)  

# Get LODR annotations for adding to plots
LODR.annot.ERCC <- printLODRres(expDat, dynRangeDat,
                                    lodr.res = lodr.ERCC)

###############################################################################
# Annotate the dynamic range plots to show how LODR is used to estimate spike 
#  in concentration
dynRangePlotLODR(dynRangeRes = dynRangeDat$dynRangePlotRes,
                 LODR.annot.ERCC = LODR.annot.ERCC)

###############################################################################

# Generate MA plots of erccs coded by concentrations from LODR
maPlotAB = maConcPlot(expDat, countPair = expDat$expressDat_l,LODR.annot.ERCC,
                      alphaPoint = 0.8,r_mAdjust = T, replicate = T,
                      ReplicateName="Library")

###############################################################################
# Generate ROC curves for the selected differential ratios
erccROC.res = erccROC(expDat)

###############################################################################
saveResults(expDat, erccROC.res, maPlotAB, lodr.ERCC)

###############################################################################
### End of Script