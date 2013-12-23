#siteName = "MAY"; platform = "ILM"; analysis = "WEHI"; sample1 = "A";sample2 = "B"


#
# # Analysis of WEHI ERCC reads data
# erccSubread = read.table("/Users/sarahmunro/Downloads/16runs-3sites-ERCC/ERCC-Raw-Count-Run16.table")
# erccSubreadMAY = erccSubread[,c(grep(pattern= "MAY",x=colnames(erccSubread)))]
# erccSubreadBGI = erccSubread[,c(grep(pattern= "BGI",x=colnames(erccSubread)))]
# erccSubreadCNL = erccSubread[,c(grep(pattern= "CNL",x=colnames(erccSubread)))]


# Set all warnings to report as they appear
options(warn = 1)

# Get required libraries
require(ggplot2)
require(reshape2)
require(plyr)
require(directlabels)
require(scales)
require(edgeR)
require(DESeq)
require(QuasiSeq)
require(genefilter)

# Source package scripts
source("~/Documents/NISTMunro/Projects/ERCCDashboard/loadExpMeas.R") # load expression measures from file
source("~/Documents/NISTMunro/Projects/ERCCDashboard/loadERCCInfo.R") # load ERCC definition files
source("~/Documents/NISTMunro/Projects/ERCCDashboard/generateDesignMat.R") # build factor matrix from sample names
source("~/Documents/NISTMunro/Projects/ERCCDashboard/compareReplicates.R") # compare technical replicates for a sample 
source("~/Documents/NISTMunro/Projects/ERCCDashboard/dynRangePlot.R") # dose response curve for dynamic range
source("~/Documents/NISTMunro/Projects/ERCCDashboard/120919 maConcPlot.R") # MA plot coded with LODR
source('~/Documents/NISTMunro/Projects/ERCCDashboard/glm_est_r_m.R') # estimate r_m
source('~/Documents/NISTMunro/Projects/ERCCDashboard/120917 LODR.R') # estimate LODR
source('~/Documents/NISTMunro/Projects/ERCCDashboard/testDE.R') # obtain pvals
source("~/Documents/NISTMunro/Projects/ERCCDashboard/erccROC.R") # ROC curves
source("~/Documents/NISTMunro/Projects/ERCCDashboard/multiplot.R") # plot multiple plots on the same page
source("~/Documents/NISTMunro/Projects/ERCCDashboard/combineTechReps.R")

Transcripts <- loadExpMeas(siteName = siteName, platform = platform, analysis = analysis)

###############################################################################################################

# Run loadERCCInfo function to obtain idCols and MixDef information
erccInfo <- loadERCCInfo()
idCols = erccInfo$idCols
MixDef = erccInfo$MixDef

TranscriptsERCCOnly = Transcripts[c(grep("ERCC-0", Transcripts$Feature)),]
# Check for and remove ERCCs in the definition file that are not in the count data file
idCols = idCols[match(TranscriptsERCCOnly$Feature,idCols$Feature),]

# Remove ERCCs without a Ratio
idCols = idCols[which(is.finite(idCols$Ratio)),]

designMat = generateDesignMat(TranscriptsERCCOnly, factorList = c("Sample","Library","Site"), patternSplit = '_')

filenameRoot = paste(siteName,analysis,platform,sample1,sample2,sep = ".")

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


designMatSum <- generateDesignMat(sampleLibeDataNorm, factorList = c("Sample","Library","Site"), patternSplit = '_')

testDE(filenameRoot = filenameRoot, cnt = sampleLibeDataNoNorm, info = designMatAB, idCols = idColsAdj, FCcode = FCcode, r_m.mn = r_m.mn, totalSeqReads = totalSeqReads, totalReads = totalReads,legendLabels=legendLabels) 

deRes <- read.csv(file = paste(filenameRoot,"quasiSeq.res.csv",sep="."))
p.thresh<-.1
if(any(deRes$qvals<choseFDR)) p.thresh<-max(deRes$pvals[deRes$qvals<choseFDR])
print("Threshold P-value")
print(p.thresh)
if (p.thresh > .1){
  print(paste("Threshold P-value is high for the chosen FDR of ", as.character(choseFDR)))
  print("The sample comparison indicates a large amount of differential expression in the measured transcript populations")
  print("Resetting Threshold to 0.1")
  p.thresh<- 0.1
}

# Find LODR estimates using the ERCC data pvalues
lodr.ERCC = findLODR(pval.cutoff=p.thresh, prob=.9,filenameRoot=filenameRoot, kind="ERCC", FCcode = FCcode, legendLabels=legendLabels)