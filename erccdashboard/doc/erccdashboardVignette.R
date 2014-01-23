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
ERCCmixes = "RatioPair"
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
                  sample2Name, ERCCmixes, ERCCdilution, spikeVol, totalRNAmass,
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
### code chunk number 17: dynRangeLODRAnnot
###################################################
expDat <- annotLODR(expDat)


###################################################
### code chunk number 18: printPanelA
###################################################
expDat$Figures$dynRangePlot


###################################################
### code chunk number 19: printPanelB
###################################################
expDat$Figures$rocPlot


###################################################
### code chunk number 20: printPanelC
###################################################
expDat$Figures$lodrERCCPlot


###################################################
### code chunk number 21: printPanelD
###################################################
expDat$Figures$maPlot


###################################################
### code chunk number 22: savePlots
###################################################
savePlots(expDat)


###################################################
### code chunk number 23: saveExpDat
###################################################
save(expDat,file=paste0(expDat$sampleInfo$filenameRoot,".RData"))


###################################################
### code chunk number 24: defineMainStudyInputs
###################################################
load(file = system.file("data/SEQC.Main.Example.Simple.RData", 
                        package = "erccdashboard"))
countTable <- Lab5.ILM.UHRR.HBRR.countTable
totalReads <- Lab5.ILM.UHRR.HBRR.totalReads
filenameRoot = "Lab5"
sample1Name = "UHRR"
sample2Name = "HBRR"
ERCCMixes = "RatioPair"
ERCCdilution = 1
spikeVol = 50
totalRNAmass = 2.5*10^(3)
choseFDR = 0.01
expDat <- initDat(countTable, totalReads, filenameRoot, sample1Name,
                  sample2Name, ERCCMixes, ERCCdilution, spikeVol, totalRNAmass,
                  choseFDR)


###################################################
### code chunk number 25: estRMForLab5
###################################################
expDat <- est_r_m(expDat)


###################################################
### code chunk number 26: geneExprTestLab5
###################################################
expDat <- geneExprTest(expDat)


###################################################
### code chunk number 27: getPlots
###################################################
expDat <- erccROC(expDat)
expDat <- estLODR(expDat,kind = "ERCC", prob=0.9)
expDat <- estLODR(expDat, kind = "Sim", prob = 0.9)  
expDat <- dynRangePlot(expDat,errorBars=T)
expDat <- annotLODR(expDat)
#expDat <- maSignal(expDat, alphaPoint = 0.8,  r_mAdjust = T, replicate = T)


###################################################
### code chunk number 28: printPanelASEQC
###################################################
expDat$Figures$dynRangePlot


###################################################
### code chunk number 29: printPanelBSEQC
###################################################
expDat$Figures$rocPlot


###################################################
### code chunk number 30: printPanelCSEQC
###################################################
expDat$Figures$lodrERCCPlot


###################################################
### code chunk number 31: printPanelDSEQC
###################################################
expDat$Figures$maPlot


###################################################
### code chunk number 32: savePlots
###################################################
savePlots(expDat)


###################################################
### code chunk number 33: saveResults
###################################################
save(expDat,file=paste0(expDat$sampleInfo$filenameRoot,".RData"))


###################################################
### code chunk number 34: sessionData
###################################################
sessionInfo()


