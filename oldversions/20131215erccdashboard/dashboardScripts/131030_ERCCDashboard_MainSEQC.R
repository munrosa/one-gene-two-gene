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
source("R/loadExpMeasFromRaw.R") # load expression measures from file
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
source("R/saveInputs.R")
source("R/saveResults.R")


###############################################################################
# Initialize expDat list
###############################################################################
# Input parameter adjustments
# For platform = "ILM" use analysis = "NCTR", siteName may be equal to 
# "MAY","BGI","CNL","NVS","AGR",or "COH"
# For platform = "LIF" use analysis = "NIST", siteName may be equal to
# "SQW","NWU",or "PSU"
# ERCC dilution = the dilution of the Ambion ERCC mixes prior to spiking, if no 
# dilution was performed then ERCCdilution = 1
# spikeVol = the amount of diluted ERCC mix spiked into the total RNA sample, 
# units are in microliters
# totalRNAmass = mass of total RNA spiked with the diluted ERCC mix, 
# units are in micrograms
###############################################################################
expDat <- initDat(study = "SEQC_Main", siteName = "Lab1", platform = "ILM",
                  analysis = "NCTR", sample1Name = "UHRR", sample2Name = "HBRR",
                  choseFDR = 0.05, ERCCdilution = 1, spikeVol = 50, 
                  totalRNAmass = 2.5*10^(3), printPDF = F, DEtest = F,
                  totalSeqReads = T, libeSizeNorm = T, myYLimMA = c(-3.5,3.5),
                  myXLim = c(-10,15),myYLim = NULL)

###############################################################################
# Run loadERCCInfo function to obtain ERCC information
expDat <- loadERCCInfo(expDat, erccMixes="Ambion4plexPair")

###############################################################################
# Run loadExpMeas function to process count data file
expDat <- loadExpMeasFromRaw(expDat)

### To Do Sarah Munro set this up so that we can provide the smallest possible
### set of data to users 

###############################################################################
# library size normalize the data
expDat <- libeSizeNorm(expDat)

###############################################################################
# length normalize the ERCC concentrations
expDat <- prepERCCDat(expDat)

###############################################################################
# Revert to raw counts to sum fluidic replicates, exclude Library 5
expDat = sumTechReps(expDat,libeList = levels(expDat$designMatAB$Library)[-5])

###############################################################################
### Save input data, only use for SEQC data prep 
saveInputs(expDat)

###############################################################################

### Estimate r_m for the sample pair using a negative binomial glm
expDat <- est_r_m(expDat, cnt = expDat$expressDatSumNoNorm, printPlot = T)

###############################################################################

### Evaluate dynamic range
#print("Signal-Abundance Plots for dynamic range estimation...")

dynRangeDat = dynRangePlot(expDat, expressDat = expDat$expressDatSumNorm,
                           designMat = expDat$designMatSum, noErrorBars = F)

###############################################################################

expDat <- estMnLibeFactor(expDat, cnt = expDat$expressDatSumNoNorm)

expDat <- meltedExpDat(expDat, cnt = expDat$expressDatSumNoNorm, 
                       designMat = expDat$designMatSum)

# Test for differential expression if DEtest == T
# else use existing pvalues of the data in csv files
expDat <- geneExprTest(expDat, cnt = expDat$expressDatSumNoNorm,
                       designMat = expDat$designMatSum )

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