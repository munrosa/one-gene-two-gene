
#corrPlots
require(ggplot2)
require(reshape2)
source("R/multiplot.R") # plot multiple plots on the same page

bgi = read.csv("Data/MainSEQC/BGI.NCTR.ILM.E.F.expressDatMeanERCC.csv")
may = read.csv("Data/MainSEQC/MAY.NCTR.ILM.E.F.expressDatMeanERCC.csv")
cnl = read.csv("Data/MainSEQC/CNL.NCTR.ILM.E.F.expressDatMeanERCC.csv")
bgi = bgi[-c(1)]
may = may[-c(1)]
cnl = cnl[-c(1)]
bgiE = bgi[-c(3)]
mayE = may[-c(3)]
cnlE = cnl[-c(3)]
bgiF = bgi[-c(2)]
mayF = may[-c(2)]
cnlF = cnl[-c(2)]

corrPlotsE = merge(bgiE,merge(mayE,cnlE,by="Feature"), by = "Feature")

corrPlotsF = merge(bgiF,merge(mayF,cnlF,by = "Feature"),by = "Feature")

names(corrPlotsF) = gsub("SEQC_ILM_",replacement="",x=names(corrPlotsF))
names(corrPlotsE) = gsub("SEQC_ILM_",replacement="",x=names(corrPlotsE))

pairedData <- function(expressDat = expressDat) {
  
  myData = expressDat
  myNumData = myData[-c(1)]
  # get the columns of interest based on the colSelect input
  #head(colSelect$countSet)
  
  #myNumData = myNumData[c(match(colSelect$countSet, colnames(myNumData)))]
  print(head(myNumData))
  # Get pairwise combinations of data
  getCombos = combn(as.vector(log2(myNumData)), m = 2, simplify=F)
  
  # Initialize pairCorrDat dataframe
  namedCombo = paste(names(getCombos[[1]])[1],names(getCombos[[1]])[2],sep="~")
  
  #m = gregexpr((paste(pairColour,collapse = "|")), namedCombo)
  #pairName = paste(regmatches(namedCombo,m)[[1]],collapse = "~")
  
  pairCorrDat = data.frame(Feature = myData[c(1)], Pair = namedCombo, getCombos[[1]])
  names(pairCorrDat)[3:4] = c("Counts.1", "Counts.2" )  
  
  # Fill in the rest of the pairCorrDat dataframe with a for loop
  #if (idx >= 2){
    for (i in 2:length(getCombos)){
      namedCombo = paste(names(getCombos[[i]])[1],names(getCombos[[i]])[2],sep="~")
      
      # m = gregexpr((paste(pairColour,collapse = "|")), namedCombo)
      # pairName = paste(regmatches(namedCombo,m)[[1]],collapse = "~")
      
      tempFrame = data.frame(Feature = myData[c(1)],  Pair = namedCombo, getCombos[[i]])
      names(tempFrame)[3:4] = c("Counts.1","Counts.2")
      pairCorrDat = rbind(pairCorrDat, tempFrame)
    }
  #}
  plotLabs = colsplit(pairCorrDat$Pair,pattern = "~", names = c("s1","s2"))
  pairCorrDat = cbind(pairCorrDat, plotLabs)
  
  return(pairCorrDat)
}

pairCorrDatE <- pairedData(expressDat = corrPlotsE)
pairCorrDatF <- pairedData(expressDat = corrPlotsF)

ERCCratios = read.csv(file="Data/ERCCMix1and2SEQC.csv")
#corrPlots_l$Ratio = ERCCratios$Subpool[match(corrPlots_l$Feature, ERCCratios$ERCC.AMB.Expected)]
pairCorrDatE$Ratio = ERCCratios$Subpool[match(pairCorrDatE$Feature, ERCCratios$ERCC.AMB.Expected)]
pairCorrDatF$Ratio = ERCCratios$Subpool[match(pairCorrDatF$Feature, ERCCratios$ERCC.AMB.Expected)]

FCcode = data.frame(Ratio = c("a","b","c","d"),FC =  c(4,1,.667,.5))
legendLabels = c("4:1","1:1","1:1.5","1:2")

myColors <- c("#339966","#FF9900", "#6699CC", "#CC6666")
names(myColors) <- levels(FCcode$Ratio)

colScale <- scale_colour_manual(name = "Ratio",values = myColors, labels = legendLabels)

fillScale <- scale_fill_manual(name = "Ratio", values = myColors, labels = legendLabels)

theme_set(theme_bw(base_size=16))

eplot = ggplot(pairCorrDatE) + geom_point(aes(x = Counts.1, y = Counts.2, colour = Ratio), size = 6, alpha = 0.8) + facet_wrap(s1~s2) + colScale + xlab("Log 2 Counts for SEQC E") + ylab("Log 2 Counts for SEQC E")
print(eplot)

fplot = ggplot(pairCorrDatF) + geom_point(aes(x = Counts.1, y = Counts.2, colour = Ratio), size = 6, alpha = 0.8) + facet_wrap(s1~s2) + colScale + xlab("Log 2 Counts for SEQC F") + ylab("Log 2 Counts for SEQC F")
print(fplot)

multiplot(eplot,fplot,cols=1)
