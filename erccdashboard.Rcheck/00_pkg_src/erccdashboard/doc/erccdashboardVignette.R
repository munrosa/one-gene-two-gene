### R code from vignette source 'erccdashboardVignette.Rnw'

###################################################
### code chunk number 1: erccdashboardVignette.Rnw:14-17
###################################################
options(width=60, continue = "  ")
#options(SweaveHooks=list(fig=function()
#                par(mar=c(5.1,4.1,1.1,2.1))))


###################################################
### code chunk number 2: initializeData
###################################################
library( "erccdashboard" )


###################################################
### code chunk number 3: loadRatData
###################################################
load(file = system.file("data/SEQC.RatTox.Example.RData", 
                        package = "erccdashboard"))


###################################################
### code chunk number 4: inspectRatData
###################################################
head(COH.RatTox.ILM.MET.CTL.countTable)
COH.RatTox.ILM.MET.CTL.totalReads


###################################################
### code chunk number 5: defineInputData
###################################################
countTable <- COH.RatTox.ILM.MET.CTL.countTable
totalReads <- COH.RatTox.ILM.MET.CTL.totalReads


###################################################
### code chunk number 6: definefileRoot
###################################################
filenameRoot = "COH.ILM"


###################################################
### code chunk number 7: defineInputSpikes
###################################################
sample1Name = "MET"
sample2Name = "CTL"
ERCCMixes = "Ambion4plexPair"
ERCCdilution = 1/100
spikeVol = 1
totalRNAmass = 0.500


###################################################
### code chunk number 8: defineFDR
###################################################
choseFDR = 0.1


###################################################
### code chunk number 9: initializeData
###################################################
expDat <- initDat(countTable, totalReads, filenameRoot, sample1Name,
                  sample2Name, ERCCMixes, ERCCdilution, spikeVol, totalRNAmass,
                  choseFDR)


###################################################
### code chunk number 10: initializeData
###################################################
summary(expDat)


###################################################
### code chunk number 11: est_r_m
###################################################
expDat <- est_r_m(expDat)


###################################################
### code chunk number 12: geneExprTest
###################################################
expDat <- geneExprTest(expDat)


###################################################
### code chunk number 13: erccROC
###################################################
expDat = erccROC(expDat)


###################################################
### code chunk number 14: estLODRERCC
###################################################
expDat = estLODR(expDat,kind = "ERCC", prob=0.9)


###################################################
### code chunk number 15: estLODRSim
###################################################
expDat = estLODR(expDat, kind = "Sim", prob = 0.9)  


###################################################
### code chunk number 16: dynRangeDat
###################################################
expDat <- dynRangePlot(expDat, errorBars = T)


###################################################
### code chunk number 17: printPanelD
###################################################
expDat$Figures$plotdynRange


###################################################
### code chunk number 18: dynRangeLODRAnnot
###################################################
expDat <- annotLODR(expDat)


###################################################
### code chunk number 19: printPanelA
###################################################
expDat$Figures$ROCplot


###################################################
### code chunk number 20: printPanelB
###################################################
expDat$Figures$plotLODR.ERCC


###################################################
### code chunk number 21: printPanelC
###################################################
expDat$Figures$dispPlot


###################################################
### code chunk number 22: printPanelD
###################################################
expDat$Figures$plotdynRangeAnnot


###################################################
### code chunk number 23: printPanelE
###################################################
expDat$Figures$plotERCCeffects


###################################################
### code chunk number 24: printPanelF
###################################################
expDat$Figures$plotRatioAnnot


###################################################
### code chunk number 25: savePlots
###################################################
savePlots(expDat)


###################################################
### code chunk number 26: saveResults
###################################################
saveResults(expDat)


###################################################
### code chunk number 27: defineMainStudyInputs
###################################################
load(file = system.file("data/SEQC.Main.Example.RData", 
                        package = "erccdashboard"))
countTable <- Lab5.ILM.UHRR.HBRR.countTable
totalReads <- Lab5.ILM.UHRR.HBRR.totalReads
filenameRoot = "Lab5"
sample1Name = "UHRR"
sample2Name = "HBRR"
ERCCMixes = "Ambion4plexPair"
ERCCdilution = 1
spikeVol = 50
totalRNAmass = 2.5*10^(3)
choseFDR = 0.1
expDat <- initDat(countTable, totalReads, filenameRoot, sample1Name,
                  sample2Name, ERCCMixes, ERCCdilution, spikeVol, totalRNAmass,
                  choseFDR)


###################################################
### code chunk number 28: est_r_m
###################################################
expDat <- est_r_m(expDat)


###################################################
### code chunk number 29: geneExprTest
###################################################
expDat <- geneExprTest(expDat)


###################################################
### code chunk number 30: erccROC
###################################################
expDat = erccROC(expDat)


###################################################
### code chunk number 31: estLODRERCC
###################################################
expDat = estLODR(expDat,kind = "ERCC", prob=0.9)


###################################################
### code chunk number 32: estLODRSim
###################################################
expDat = estLODR(expDat, kind = "Sim", prob = 0.9)  


###################################################
### code chunk number 33: dynRangeDat
###################################################
expDat <- dynRangePlot(expDat,errorBars=T)


###################################################
### code chunk number 34: printPanelD
###################################################
expDat$Figures$plotdynRange


###################################################
### code chunk number 35: dynRangeLODRAnnot
###################################################
expDat <- annotLODR(expDat)


###################################################
### code chunk number 36: printPanelASEQC
###################################################
expDat$Figures$ROCplot


###################################################
### code chunk number 37: printPanelBSEQC
###################################################
expDat$Figures$plotLODR.ERCC


###################################################
### code chunk number 38: printPanelCSEQC
###################################################
expDat$Figures$dispPlot


###################################################
### code chunk number 39: printPanelDSEQC
###################################################
expDat$Figures$plotdynRangeAnnot


###################################################
### code chunk number 40: printPanelESEQC
###################################################
expDat$Figures$plotERCCeffects


###################################################
### code chunk number 41: printPanelFSEQC
###################################################
expDat$Figures$plotRatioAnnot


###################################################
### code chunk number 42: savePlots
###################################################
savePlots(expDat)


###################################################
### code chunk number 43: saveResults
###################################################
saveResults(expDat)


