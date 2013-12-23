##############################################################################
# Analysis of per site results from SEQC from a given analysis pipeline
##############################################################################
# siteName = "MAY";platform = "ILM";analysis = "NCTR";libeSizeNorm = T;printPDF = F; DEtest = T; totalSeqReads = F; sample1 = "A";sample2 = "B";  FCcode = data.frame(Ratio = c("a","b","c","d"),FC =  c(4,1,.667,.5)); myYLimMA = c(-3,3); myXLim = c(-2,18); choseFDR = 0.05; legendLabels = c("4:1","1:1","1:1.5","1:2")

# A:B ratios are default FCcode matrix
# A:C 
# FCcode = data.frame(Ratio = c("a","b","c","d"),FC = c(1.23,1,.89,.80))
# B:D 
# FCcode = data.frame(Ratio = c("a","b","c","d"),FC = c(0.57,1,1.09,1.14))

##Create a custom color scale
#myColors <- c("#339966","#FF9900", "#6699CC", "#CC6666")
#names(myColors) <- levels(FCcode$Ratio)
#colScale <- scale_colour_manual(name = "Ratio",values = myColors, breaks = breakSet)
#fillScale <- scale_fill_manual(name = "Ratio", values = myColors, breaks = breakSet)

###############################################################################
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
###############################################################################################################

# Run loadExpMeas function to process data file
Transcripts <- loadExpMeas(siteName = siteName, platform = platform, analysis = analysis)

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

print("ERCC % of total reads (endo + ERCC):")
ERCCofTotal = ERCCMappedReadSums/TranscriptMappedReadSums
print(head(ERCCofTotal*100))

# Check for and remove ERCCs in the definition file that are not in the count data file
idCols = idCols[match(TranscriptsERCCOnly$Feature,idCols$Feature),]

# Remove ERCCs without a Ratio
idCols = idCols[which(is.finite(idCols$Ratio)),]

# Remove ERCCs from count data and idCols that are not present in the experiment
TranscriptsERCCOnly = TranscriptsERCCOnly[match(idCols$Feature,TranscriptsERCCOnly$Feature),]
Transcripts = rbind(TranscriptsERCCOnly, TranscriptsHumanOnly)

###############################################################################################################

# Generate the Design Matrix for the table example name is SEQC_ILM_BGI_A_1_L01_ATCACG_AC0AYTACXX
if (analysis == "NCTR"){
  designMat = generateDesignMat(TranscriptsERCCOnly, factorList = c("Study","Platform","Site","Sample","Library","Lane","Barcode","Flowcell"), patternSplit = '_')
}
if(analysis == "NIST"){
 designMat = generateDesignMat(TranscriptsERCCOnly, factorList = c("Study","Platform","Site","Sample","Library","Lane","Barcode","Flowcell","SiteInfo1","SiteInfo2"), patternSplit = '_')
}
if (analysis == "MAGIC"){
  designMat = generateDesignMat(TranscriptsERCCOnly, factorList = c("Study","Platform","Site","Sample","Library","Lane","Barcode","Flowcell","SiteInfo1","SiteInfo2"), patternSplit = '_')
}

###############################################################################################################
# Get total Reads (from the original fastq files)
if(totalSeqReads == T){
mainReads = read.delim(paste("all_qc_results/",platform,"_",siteName,"_qc_results.txt", sep=""))
mainReads = mainReads[c(1:9)]
print(names(mainReads))
}

###############################################################################################################
filenameRoot = paste(siteName,analysis,platform,sample1,sample2,sep = ".")
# Open PDF file to write results
if (printPDF == T){
  pdf(file = paste(filenameRoot,"%03d","pdf",sep="."),onefile=T,width=7,height = 7)  
}

###############################################################################################################

# Subset to get just the A and B samples for the Count Matrix and for the designMatrix
select = subset(designMat, (Sample == sample1)|(Sample == sample2))
select <- as.data.frame(lapply(select,as.character))
select <- as.data.frame(lapply(select,as.factor))
designMatAB = select
dataAB = Transcripts[-c(1)]
TranscriptsAB = cbind(Transcripts[c(1)],dataAB[c(match(select$countSet, names(dataAB)))])

# write csv file to directory
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
  
  totalReads = mainReads$numReads[c(match(select$countSet,mainReads$sampleName))]
  
  libeSize = totalReads
  datCols = expressDat[-c(1)]
  libeSize = libeSize/(10^6) #per million mapped reads
  #Library size normalize the data  
  libAdjust = sweep(datCols, 2, libeSize,"/")
  expressDat = cbind(expressDat[c(1)], libAdjust)
}

print("Library sizes:")
print(libeSize)

if(analysis == "MAGIC"){
  #sinhMagic = apply(expressDat[-c(1)],MARGIN = 2, FUN = sinh)
  # First, transform the data back to normalized counts from asinh normalized counts
  # Idx = log2(z + sqrt((1+(z^2)))) - 1,
  
  transformMagic <-function(Idx){
    
    #zCounts = sinh(Idx)
    ##original Magic eqn most likely incorrect
    #zCounts = sqrt( ( ( ( 2^( Idx+1 ) )^2 ) - 1 )/2 ) 
    
    ##incorrect Magic eqn undefined at 0
    #zCounts = sqrt( ( ( ( ( 2^( Idx ) ) + 1 )^2)-1)/2)
    
    ## rederived correct Magic
    zCounts = (( (2^(Idx+1))*(2^(Idx+1)) )-1 )/ (2*(2^(Idx+1)))
    
    ## rederived incorrect Magic eqn
    #zCounts = ((((2^(Idx)) + 1)^2) -1)/(2*((2^(Idx))+1))
    
    #zlog2Trans = log2(zCounts)
    #return(zlog2Trans)
    return(zCounts)
  }
  log2Magic = apply(expressDat[-c(1)],MARGIN = 2, FUN = transformMagic)
  
  expressDat = cbind(expressDat[c(1)],log2Magic)
  #expressDat = magicNoTransform
  TranscriptsAB = expressDat
}

# get just ERCC data in the expression data frame
expressDat = expressDat[c(grep("ERCC-0", expressDat$Feature)),]

if (analysis != "MAGIC"){
  #Length normalize the Expected ERCC concentrations
  lengthFactor = (idCols$Length)/(1000)
  
  # If length normalization of the expected concentrations is desired (default)
  idCols$Conc1 = (idCols$Conc1*lengthFactor)
  idCols$Conc2 = (idCols$Conc2*lengthFactor)  
  
}

###############################################################################################################
# Repeatability for Lanes and Libraries
###############################################################################################################

# Evaluate possible lane and flow cell effects
print("Evaluating Repeatability...")
print("lane comparisons of sample-library measurements")
# Evaluate the deviance of Lanes from mean values, calculate the mean for 16 lane flow cell combinations for Sample A Library 1 and Sample B Library 1
print("SEQC A")

print("LIBRARY 1")
select = subset(designMat, (Sample == sample1)&(Library == "1"))
SEQCAmeanSampleL1 = compareReplicates(repDat = expressDat, idCols, siteName, select, repType = "Lane", SampleName = sample1)


print("SEQC B")
select = subset(designMat, (Sample == sample2)&(Library == "1"))
print("LIBRARY 1")
SEQCBmeanSampleL1 = compareReplicates(repDat = expressDat, idCols, siteName, select, repType = "Lane", SampleName = sample2)


SEQCAmeanLanePlotL1 = SEQCAmeanSampleL1[[1]]
SEQCBmeanLanePlotL1 = SEQCBmeanSampleL1[[1]]

SEQCAmeanData = SEQCAmeanSampleL1[[2]]
SEQCBmeanData = SEQCBmeanSampleL1[[2]]

SEQCAmeanData$Sample = sample1
SEQCBmeanData$Sample = sample2

# Print the plots to the pdf
print(SEQCAmeanLanePlotL1)

print(SEQCBmeanLanePlotL1)

# If no flow cells or lanes looked like outliers we can sum across the entire experiment to obtain sample-library pairs and then compare library replicates

# Revert to raw counts to sum fluidic replicates and then sum by library replicate
expressDatSummed = combineTechReps(type = "sum",libeSizeNorm = T, expressDat = TranscriptsAB, designMat = designMatAB, sampleNameList = c(sample1,sample2), libeList = levels(designMatAB$Library))

expressDatSummed = expressDatSummed[[1]]

# Generate the Design Matrix for the table example name is SEQC_ILM_BGI_A_1
designMatSum = generateDesignMat(expressDatSummed, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')

# Evaluate the deviance of Libraries from mean values, calculate the mean for 5 library preparations for Sample A Library 1 and Sample B Library 1

select = subset(designMatSum, (Sample == sample1))

SEQCAmeanSample = compareReplicates(repDat = expressDatSummed, idCols, siteName, select, repType = "Library", SampleName = sample1)

select = subset(designMatSum, (Sample == sample2))

SEQCBmeanSample = compareReplicates(repDat = expressDatSummed, idCols, siteName, select, repType = "Library", SampleName = sample2)

SEQCAmeanSamplePlot = SEQCAmeanSample[[1]]
SEQCBmeanSamplePlot = SEQCBmeanSample[[1]]
SEQCAmeanData = data.frame(SEQCAmeanSample[[2]])
SEQCBmeanData = data.frame(SEQCBmeanSample[[2]])
SEQCAmeanData$Sample = sample1
SEQCBmeanData$Sample = sample2

multiplot(SEQCAmeanSamplePlot,SEQCBmeanSamplePlot, cols=2)

# Now look at the 4 libraries without Library 5
select = subset(designMatSum, (Sample == sample1)&(Library != "5"))
select <- as.data.frame(lapply(select,as.character))
select <- as.data.frame(lapply(select,as.factor))

SEQCAmeanSample = compareReplicates(repDat = expressDatSummed, idCols, siteName, select, repType = "Library", SampleName = sample1)

select = subset(designMatSum, (Sample == sample2)&(Library!= "5"))
select <- as.data.frame(lapply(select,as.character))
select <- as.data.frame(lapply(select,as.factor))

SEQCBmeanSample = compareReplicates(repDat = expressDatSummed, idCols, siteName, select, repType = "Library", SampleName = sample2)

SEQCAmeanSamplePlot = SEQCAmeanSample[[1]]
SEQCBmeanSamplePlot = SEQCBmeanSample[[1]]
SEQCAmeanData = data.frame(SEQCAmeanSample[[2]])
SEQCBmeanData = data.frame(SEQCBmeanSample[[2]])
SEQCAmeanData$Sample = sample1
SEQCBmeanData$Sample = sample2

multiplot(SEQCAmeanSamplePlot,SEQCBmeanSamplePlot, cols=2)

###############################################################################################################
# Dynamic Range Plots of Sample-Library Combinations
###############################################################################################################

print("Signal-Abundance Plots for dynamic range estimation...")

# Revert to raw counts to sum fluidic replicates and then sum by library replicate
# exclude library 5

sampleLibeData = combineTechReps(type = "sum", libeSizeNorm = T, expressDat = TranscriptsAB, designMat = designMatAB, sampleNameList = c(sample1,sample2), libeList = levels(designMatAB$Library)[-5])

sampleLibeDataNorm = sampleLibeData[[1]]

sampleLibeDataNoNorm = sampleLibeData[[2]]

sampleLibeSums = colSums(sampleLibeDataNoNorm[-c(1)],na.rm =T)

print(head(sampleLibeDataNorm))
print(head(sampleLibeDataNoNorm))

# create Data frame of the Sample library pairs to return to the workspace if desired
write.csv(sampleLibeDataNoNorm,file=paste(siteName,analysis,platform,sample1,sample2,"samplelibedataNoLibeSize.csv",sep = "."))

# create Data frame of the Sample library pairs to return to the workspace if desired
write.csv(sampleLibeDataNorm,file=paste(siteName,analysis,platform,sample1,sample2,"samplelibedataNorm.csv",sep = "."))

designMatNorm = generateDesignMat(sampleLibeDataNorm, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')

expressDatSampleMean = combineTechReps(type = "mean", libeSizeNorm = T, expressDat = sampleLibeDataNorm, designMat = designMatNorm, sampleNameList = c(sample1,sample2), libeList = levels(designMatAB$Library)[-5])

expressDatSampleMean = expressDatSampleMean[[1]]
expressDatMeanERCC = expressDatSampleMean[c(grep("ERCC-0", expressDatSampleMean$Feature)),]

write.csv(expressDatMeanERCC,file=paste(filenameRoot,"expressDatMeanERCC.csv",sep = "."))


###############################################################################################################
#if (analysis != "MAGIC"){
  ### Estimate r_m for the sample pair using a negative binomial glm
  r_m.res = glm_est_r_m(site = filenameRoot, cnt = sampleLibeDataNoNorm, libeSizeNorm = libeSizeNorm,totalSeqReads = totalSeqReads,totalReads = totalReads, idCols = idCols, FCcode = FCcode,sample1 = sample1, sample2 = sample2,legendLabels=legendLabels)
  
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
  
  # Return the r_m adjusted concentrations in idCols
  idColsAdj = r_m.res$idCols
  
  #expressDatAandBERCCs = merge(idCols[c(1,4)],expressDatSampleAB)
  r_m.mnlog = exp(r_m.mn)
  
  designMatSum <- generateDesignMat(sampleLibeDataNorm, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')
  
  dynRangeDat = dynRangePlot(expressDat = sampleLibeDataNorm, designMat = designMatSum,idCols = idColsAdj, sampleNames = c(sample1,sample2),myXLim=myXLim,myYLim=NULL, legendLabels=legendLabels,FCcode=FCcode)
  
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
  
  #####
  # After observing that Library 5 is an outlier it is excluded from further analysis
  # The four remaining library replicates for each sample are collected into a data frame with the header names:
  # Feature Ratio Sample1 Sample2 Library
  #####
  myDataERCC = sampleLibeDataNorm
  expressDat = myDataERCC[-c(1)] 
  sampleNameList = c(sample1,sample2)
  libenum = c("1","2","3","4") # Library 5 is not included
  #if(siteName == "NWU")libenum = c("1","2","3")
  expressDatAandBERCCs = merge(idColsAdj[c(1,4)],myDataERCC)
  expressDatAandBERCCs$Feature <- as.factor(as.character(expressDatAandBERCCs$Feature))
  expressDat = expressDatAandBERCCs[-c(1:2)]
  expressDatSampleAB = expressDatAandBERCCs[c(1:2)]
  
  # Create long data frame of counts for each library in a for loop
  for (libe in 1:length(libenum)){
    if(libe == 1){
      select = subset(designMatSum, (Library == libenum[libe]))
      select <- as.data.frame(lapply(select,as.character))
      select <- as.data.frame(lapply(select,as.factor))
      expressDatForSum = expressDat[c(match(select$countSet, names(expressDat)))]
      expressDatSampleAdd = cbind(expressDatSampleAB, expressDatForSum)
      expressDatSampleAdd$Library = libe
      names(expressDatSampleAdd)[3:4]= c(sample1,sample2)
    }else{
      select = subset(designMatSum, (Library == libenum[libe]))
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
  
  print(str(expressDatSampleAdd))
  
  # Test for differential expression if DEtest == T
  # else use existing pvalues of the data in csv files
  if (DEtest == T){
    testDE(filenameRoot = filenameRoot, cnt = sampleLibeDataNoNorm, info = designMatAB, idCols = idColsAdj, FCcode = FCcode, r_m.mn = r_m.mn, totalSeqReads = totalSeqReads, totalReads = totalReads,legendLabels=legendLabels)
  } 
  
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
  
  # Find LODR estimates using data pvalues simulated from endogenous transcripts
  lodr.Sim = findLODR(pval.cutoff=p.thresh, prob=.9,filenameRoot=filenameRoot,kind="Sim", FCcode = FCcode, legendLabels=legendLabels)  
  
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
  
  maPlotABnoLODR = maConcPlot(idCols = idColsAdj,siteName = siteName, analysis = analysis, countPair = expressDatSampleAdd, r_m.mn = r_m.mnlog, sample1 = sample1,sample2 = sample2, myYLim = myYLimMA, myXLim=myXLim, replicate = T, cutoffs = NULL,FCcode = FCcode, legendLabels=legendLabels)
  
  maPlotAB = maConcPlot(idCols=idColsAdj,siteName = siteName, analysis = analysis, countPair = expressDatSampleAdd, r_m.mn = r_m.mnlog, sample1 = sample1,sample2 = sample2, myYLim = myYLimMA, myXLim=myXLim, replicate = T, cutoffs = ConcEst[which(!(is.na(ConcEst)))], FCcode = FCcode, legendLabels=legendLabels)

  print("Plotting MA plot with coding")
  print(maPlotAB[[1]])
  sdGlobal = maPlotAB$sdGlobal  
#print("Plotting Ma plot without coding")
  #print(maPlotABnoLODR[[1]])
  
  erccROC.res = erccROC(filenameRoot = filenameRoot, kind = "ERCC", folds = FCcode, legendLabels = legendLabels)
  print(erccROC.res$ROCplot)
  AUCdat = erccROC.res$AUCdat
  
  
  
# }else{
#   # If Magic data is used - or any data that is already normalized...
#   r_m.mnlog = 1
#   r_m.mn = log(1)
#   
#   designMatSum = generateDesignMat(sampleLibeData, factorList = c("Study","Platform","Site","Sample","Library"), patternSplit = '_')
#   
#   write.csv(sampleLibeData,file=paste(siteName,analysis,platform,"samplelibedata.csv",sep = "."))
#   
#   dynRangeDat = dynRangePlot(expressDat = sampleLibeData, designMat = designMatSum,idCols = idColsAdj, sampleNames = c(sample1,sample2), siteName = siteName, analysis = analysis, xlimEffects=c(-20,20))
#   
#   print(dynRangeDat$fit.coeff)
#   fit.coeff = dynRangeDat$fit.coeff
#   
#   #####
#   # After observing that Library 5 is an outlier it is excluded from further analysis
#   # The four remaining library replicates for each sample are collected into a data frame with the header names:
#   # Feature Ratio Sample1 Sample2 Library
#   #####
#   sampleLibeData = sampleLibeDataNoNorm
#   myDataERCC = sampleLibeData
#   expressDat = myDataERCC[-c(1)] 
#   sampleNameList = c(sample1,sample2)
#   libenum = c("1","2","3","4") # Library 5 is not included
#   #if (siteName == "NWU") libenum = c("1","2","3")
#   expressDatAandBERCCs = merge(idColsAdj[c(1,4)],myDataERCC)
#   expressDatAandBERCCs$Feature <- as.factor(as.character(expressDatAandBERCCs$Feature))
#   expressDat = expressDatAandBERCCs[-c(1:2)]
#   expressDatSampleAB = expressDatAandBERCCs[c(1:2)]
#   
#   # Create long data frame of counts for each library in a for loop
#   for (libe in 1:length(libenum)){
#     if(libe == 1){
#       select = subset(designMatSum, (Library == libenum[libe]))
#       select <- as.data.frame(lapply(select,as.character))
#       select <- as.data.frame(lapply(select,as.factor))
#       expressDatForSum = expressDat[c(match(select$countSet, names(expressDat)))]
#       expressDatSampleAdd = cbind(expressDatSampleAB, expressDatForSum)
#       expressDatSampleAdd$Library = libe
#       names(expressDatSampleAdd)[3:4]= c("SEQC_A","SEQC_B")
#     }else{
#       select = subset(designMatSum, (Library == libenum[libe]))
#       select <- as.data.frame(lapply(select,as.character))
#       select <- as.data.frame(lapply(select,as.factor))
#       expressDatForSum = expressDat[c(match(select$countSet, names(expressDat)))]
#       #meanCounts = data.frame(rowMeans(expressDatForSum))
#       addRows = cbind(expressDatSampleAB[c(1:2)],expressDatForSum)
#       addRows$Library = libe
#       #print(addRows)
#       names(addRows)[3:4]= c("SEQC_A","SEQC_B")
#       expressDatSampleAdd = rbind(expressDatSampleAdd,addRows )  
#     } 
#   }
#   expressDatSampleAdd$Library = as.factor(expressDatSampleAdd$Library)
#   expressDatSampleAdd$Feature = as.factor(expressDatSampleAdd$Feature)
#   
#   print(str(expressDatSampleAdd))
#   r_m.mnlog = 1
#   print(r_m.mnlog)
#   print(head(expressDatSampleAdd))
#   print(summary(expressDatSampleAdd))
#   maPlotABnoLODR = maConcPlot(idColsAdj,siteName = siteName, analysis = analysis, countPair = expressDatSampleAdd, r_m.mn = r_m.mnlog, sample1 = "SEQC A",sample2 = "SEQC B", myYLim = myYLimMA, replicate = T, cutoffs = F,FCcode=FCcode, legendLabels=legendLabels)
#   print(maPlotABnoLODR)[[1]]
#   
# }


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

rm(list = ls())
###########  END OF SCRIPT  ###############