###############################################################################
# Analysis of Rat Toxicogenomics data from SEQC
###############################################################################

sampleInfo = list(study = "SEQC_RatTox", #
                  siteName = "COH", # 
                  platform = "ILM", #
                  analysis = "RatTox", #
                  sample1Name = "NIT", # 
                  sample2Name = "CTL", #
                  DEtest = T, #
                  totalSeqReads = T, #
                  choseFDR = 0.1, #
                  printPDF = F, #
                  FCcode = data.frame(Ratio = c("a","b","c","d"), 
                                      FC =  c(4,1,.667,.5)), # 
                  myYLimMA = c(-3.5,3.5), #
                  myXLim = c(-10,15), 
                  myYLim = NULL,# 
                  legendLabels = c("4:1","1:1","1:1.5","1:2"), # 
                  ERCCdilution = 1/100, # 
                  spikeVol = 1, #
                  totalRNAmass = 0.500 # 
                  )


###############################################################################
# For NN_IP samples
# CTL has ERCC Mix 2; Treatments NIT and THI have ERCC Mix 1
# Same instrument different days for library prep and sequencing for the CTL 
#  and Treatment samples

# For NN_OG samples
# CTL has ERCC Mix 2; Treatments 3ME,MET,NAP, have ERCC Mix 1
# Same instrument same library prep and sequencing dates for the CTL and 
#  Treatment samples
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
source("R/sumTechReps.R") # combine technical replicates
#source("R/compareReps.R") # compare technical replicates for a sample 
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
source("R/meltedExpDat.R")
source("R/geneExprTest.R")
source("R/printLODRres.R")

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
# Dynamic Range Plots of Sample Biological Replicates
###############################################################################

  ### Estimate r_m for the sample pair using a negative binomial glm
  expDat <- est_r_m(sampleInfo, expDat,cnt = expDat$TranscriptsAB, printPlot = T)

###############################################################################

### Evaluate dynamic range
  #print("Signal-Abundance Plots for dynamic range estimation...")

  dynRangeDat = dynRangePlot(sampleInfo, expDat, expressDat = expDat$expressDat,
                             designMat = expDat$designMatSum,noErrorBars = F)
  
###############################################################################


  expDat <- estMnLibeFactor(sampleInfo, expDat, cnt = expDat$TranscriptsAB)

  expDat <- meltedExpDat(sampleInfo, expDat, cnt = expDat$TranscriptsAB, designMat = expDat$designMatSum)
  
  #expDat <- estMnLibeFactor(sampleInfo, expDat)
  #expDat <- meltedExpDat(sampleInfo, expDat)



  # Test for differential expression if DEtest == T
  # else use existing pvalues of the data in csv files
  expDat <- geneExprTest(sampleInfo, expDat)

# Find LODR estimates using the ERCC data pvalues
  lodr.ERCC = findLODR(sampleInfo,expDat,kind = "ERCC", prob=0.9)
  
# Find LODR estimates using data pvalues simulated from endogenous transcripts
  lodr.Sim = findLODR(sampleInfo, expDat, kind = "Sim", prob = 0.9)  

#### Edit starting here
  LODR.annot.ERCC <- printLODRres(sampleInfo, expDat, dynRangeDat,
                                      lodr.res = lodr.ERCC)
  
###############################################################################
# Annotate the dynamic range plots to show how LODR is used to estimate spike 
#  in concentration
# Move legend
  dynRangePlotLODR(dynRangeRes = dynRangeDat$dynRangePlotRes,
                   LODR.annot.ERCC = LODR.annot.ERCC, folds = sampleInfo$FCcode, 
                   legendLabels = sampleInfo$legendLabels)
 
###############################################################################

# Generate MA plots of erccs coded by concentrations from LODR
    
  maPlotAB = maConcPlot(sampleInfo, expDat, LODR.annot.ERCC, alphaPoint = 0.8,r_mAdjust = T, replicate = T)

  print("Plotting MA plot with coding")
  print(maPlotAB[[1]])
  
###############################################################################
# Generate ROC curves for the selected differential ratios
  erccROC.res = erccROC(filenameRoot = sampleInfo$filenameRoot, kind = "ERCC",
                        legendLabels = sampleInfo$legendLabels,
                        folds = sampleInfo$FCcode,idCols = expDat$idColsAdj)
  print(erccROC.res$ROCplot)
  #AUCdat = erccROC.res$AUCdat

###############################################################################
saveResults(sampleInfo, erccROC.res, maPlotAB, lodr.ERCC, expDat)

###############################################################################




###########  END OF SCRIPT  ###############

