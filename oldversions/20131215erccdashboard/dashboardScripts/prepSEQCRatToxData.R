##############################################################################
# Prepare Rat Tox Data for use in erccdashboard package
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
expDat <- initDat(study = "SEQC_RatTox", siteName = "COH", platform = "ILM",
                  analysis = "RatTox", sample1Name = "3ME", sample2Name = "CTL",
                  choseFDR = 0.1, ERCCdilution = 1/100, spikeVol = 1, 
                  totalRNAmass = 0.500, printPDF = F, DEtest = T,
                  totalSeqReads = T, libeSizeNorm = T, myYLimMA = c(-3.5,3.5),
                  myXLim = c(-10,15),myYLim = NULL)

###############################################################################
# Run loadERCCInfo function to obtain ERCC information
expDat <- loadERCCInfo(expDat, erccMixes="Ambion4plexPair")

###############################################################################
# Run loadExpMeas function to process count data file
expDat <- loadExpMeasFromRaw(expDat)

###############################################################################
### Save input data, only use for SEQC data prep 
saveInputs(expDat)