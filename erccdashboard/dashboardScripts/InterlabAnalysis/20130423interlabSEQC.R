library(reshape2)
library(ggplot2)
source("R/getDesignMat.R")
source("R/multiplot.R")
theme_set(theme_bw(base_size=16))


#pdf(file = paste("InterlabFig5pieces","pdf",sep="."),onefile=T,width=7,height = 7)

   #Create a custom color scale
  #myColors <- c("#6699CC", "#CC6666")
  myColors <- c("firebrick2", "dodgerblue3")  
  names(myColors) <- c("ILM","LIF")
  colScale <- scale_colour_manual(name = "Platform",values = myColors)
  fillScale <- scale_fill_manual(name = "Platform", values = myColors)

# FCcode <- data.frame(Ratio = c("a","b","c","d"),FC =  c(4,1,.667,.5));
# legendLabels = c("4:1","1:1","1:1.5","1:2");
# 
# myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
# #names(myColors) <- levels(legendLabels)
# myColorsDiff <- myColors[-(which(FCcode$FC == 1))]
# legendLabelsDiff <- legendLabels[-(which(FCcode$FC == 1))]
# 
# colScale <- scale_colour_manual(name = "Ratio",values = myColorsDiff, labels = legendLabelsDiff)
# fillScale <- scale_fill_manual(name = "Ratio", values = myColorsDiff, labels = legendLabelsDiff)

load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/MAY.NCTR.ILM.A.B.RData")
load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/SQW.NIST.LIF.A.B.RData")
load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/PSU.NIST.LIF.A.B.RData")
load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/NWU.NIST.LIF.A.B.RData")
load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/NVS.NCTR.ILM.A.B.RData")
load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/AGR.NCTR.ILM.A.B.RData")
load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/COH.NCTR.ILM.A.B.RData")
load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/CNL.NCTR.ILM.A.B.RData")
load("/Users/smunro/Documents/NISTMunro/Projects/ERCCDashboardDevelopment/20130617_NBTmanuscriptsubmittedresultsMicroarrayEndoCloudAnalysis/Results/130423_Results/BGI.NCTR.ILM.A.B.RData")

r_m_dat = data.frame(Site = c("BGI","CNL","MAY","AGR","COH","NVS"), Platform = "ILM")

r_m_dat2 = data.frame(Site = c("PSU","NWU","SQW"), Platform = "LIF")

r_m_dat = rbind(r_m_dat,r_m_dat2)

r_m_dat$Platform = as.factor(r_m_dat$Platform)

r_m_dat$r_m = c(BGI.NCTR.ILM.A.B.r_m,CNL.NCTR.ILM.A.B.r_m,MAY.NCTR.ILM.A.B.r_m,AGR.NCTR.ILM.A.B.r_m,COH.NCTR.ILM.A.B.r_m,NVS.NCTR.ILM.A.B.r_m,PSU.NIST.LIF.A.B.r_m,NWU.NIST.LIF.A.B.r_m,SQW.NIST.LIF.A.B.r_m)

r_m_dat$r_m_lower = c(BGI.NCTR.ILM.A.B.r_m_lower,CNL.NCTR.ILM.A.B.r_m_lower,MAY.NCTR.ILM.A.B.r_m_lower,AGR.NCTR.ILM.A.B.r_m_lower,COH.NCTR.ILM.A.B.r_m_lower,NVS.NCTR.ILM.A.B.r_m_lower,PSU.NIST.LIF.A.B.r_m_lower,NWU.NIST.LIF.A.B.r_m_lower,SQW.NIST.LIF.A.B.r_m_lower)

r_m_dat$r_m_upper = c(BGI.NCTR.ILM.A.B.r_m_upper,CNL.NCTR.ILM.A.B.r_m_upper,MAY.NCTR.ILM.A.B.r_m_upper,AGR.NCTR.ILM.A.B.r_m_upper,COH.NCTR.ILM.A.B.r_m_upper,NVS.NCTR.ILM.A.B.r_m_upper,PSU.NIST.LIF.A.B.r_m_upper,NWU.NIST.LIF.A.B.r_m_upper,SQW.NIST.LIF.A.B.r_m_upper)

labelSites <- data.frame(Original = levels(r_m_dat$Site), New = c("Lab 1","Lab 2","Lab 3","Lab 4","Lab 5","Lab 6","Lab 7","Lab 8","Lab 9"))

r_m_dat$Site <- labelSites$New[match(r_m_dat$Site,labelSites$Original)] 

r_m_dat$Site <- as.factor(as.character(r_m_dat$Site))

theme_update(legend.justification=c(0,0), legend.position=c(0,0))
empRm = 1.433
pC = ggplot(r_m_dat)+geom_errorbar(aes(x = Site,y=exp(r_m),ymin = exp(r_m_lower),ymax = exp(r_m_upper))) +geom_point(aes(x = Site, y = exp(r_m),colour = Platform),size = 6) + geom_hline(yintercept = empRm,linetype = 2)+ ylab(expression(r[m])) + colScale + annotate("text",x = 6,y = empRm - 0.02,label = paste(expression(r[m]),"==",as.character(empRm), sep=" "),parse = T) + xlab("Sequencing Site")

print(pC)
ILML5res <- load("ILMLibe5Res.RData")
LIFL5res <- load("LIFLibe5Res.RData")
# TO DO Sarah Munro
#pCL5 <- ggplot()

AUCList <- ls(pattern=glob2rx(pattern=c("*.AUC")))
AUCDatNames <-data.frame(DatName = ls(pattern=glob2rx(pattern=c("*.AUC"))))
AUCdat <- colsplit(AUCDatNames$DatName,"\\.", c("Site","Analysis","Platform","Sample1","Sample2","Data"))

for (i in 1:length(AUCList)){
  numAUC <- as.numeric(gsub("<",replacement="",as.character(get(AUCList[i])$AUC)))
  print(numAUC)
  AUCdat$x4[i] = numAUC[which(get(AUCList[i])$Ratio == "4:1")]
  AUCdat$x2[i] = numAUC[which(get(AUCList[i])$Ratio == "1:2")]
  AUCdat$x1.5[i] = numAUC[which(get(AUCList[i])$Ratio == "1:1.5")]
}
AUCdat$Site <- as.factor(AUCdat$Site)
labelSites <- data.frame(Original = levels(AUCdat$Site), New = c("Lab 1","Lab 2","Lab 3","Lab 4","Lab 5","Lab 6","Lab 7","Lab 8","Lab 9"))

AUCdat$Site <- labelSites$New[match(AUCdat$Site,labelSites$Original)] 

AUCdat$Site <- as.factor(as.character(AUCdat$Site))

AUCdat_l = melt(AUCdat)
colnames(AUCdat_l)[which(colnames(AUCdat_l) == "variable")] <- "Ratio"
levels(AUCdat_l$Ratio) <- c("4:1","1:2","1:1.5")


theme_update(legend.justification=c(0,0), legend.position=c(0,0))
pA = ggplot(AUCdat_l)+geom_point(aes(x = Site,y=value, colour = Platform,shape = Ratio),alpha = 0.8,size = 6) + ylab("AUC") + colScale + xlab("Sequencing Site")
print(pA)

######################
lodrList <- ls(pattern=glob2rx(pattern=c("*.lodr.ERCC")))
lodrDatNames <-data.frame(DatName = ls(pattern=glob2rx(pattern=c("*.lodr.ERCC"))))
LODRdat <- colsplit(lodrDatNames$DatName,"\\.", c("Site","Analysis","Platform","Sample1","Sample2","Data","Source"))
LODRdat <- LODRdat[,-c(7)]

for (i in 1:length(lodrList)){
  numLODR <- as.numeric(gsub("<",replacement="",as.character(get(lodrList[i])$Estimate)))
  #minLODR <- as.numeric(gsub("<",replacement="",as.character(get(lodrList[i])[,4])))
  #maxLODR <- as.numeric(gsub("<",replacement="",as.character(get(lodrList[i])[,5])))
   LODRdat$x4[i] = numLODR[which(get(lodrList[i])$Ratio == "4:1")]
#   LODRdat$max4[i] = maxLODR[which(get(lodrList[i])$Ratio == "4:1")]
#   LODRdat$min4[i] = minLODR[which(get(lodrList[i])$Ratio == "4:1")]
#   
   LODRdat$x2[i] = numLODR[which(get(lodrList[i])$Ratio == "1:2")]
#   LODRdat$max2[i] = maxLODR[which(get(lodrList[i])$Ratio == "4:1")]
#   LODRdat$min2[i] = minLODR[which(get(lodrList[i])$Ratio == "4:1")]
#   
   LODRdat$x1.5[i] = numLODR[which(get(lodrList[i])$Ratio == "1:1.5")]
#   LODRdat$max1.5[i] = maxLODR[which(get(lodrList[i])$Ratio == "4:1")]
#   LODRdat$min1.5[i] = minLODR[which(get(lodrList[i])$Ratio == "4:1")]
}
LODRdat$Site <- as.factor(LODRdat$Site)
labelSites <- data.frame(Original = levels(LODRdat$Site), New = c("Lab 1","Lab 2","Lab 3","Lab 4","Lab 5","Lab 6","Lab 7","Lab 8","Lab 9"))

LODRdat$Site <- labelSites$New[match(LODRdat$Site,labelSites$Original)] 

LODRdat$Site <- as.factor(as.character(LODRdat$Site))


LODRdat_l = melt(LODRdat,id.vars=c("Site","Analysis","Platform","Sample1","Sample2","Data"))
colnames(LODRdat_l)[which(colnames(LODRdat_l) == "variable")] <- "Ratio"
levels(LODRdat_l$Ratio) <- c("4:1","1:2","1:1.5")



theme_update(legend.justification=c(0,1), legend.position=c(0,1))
pB <- ggplot(LODRdat_l)+geom_point(aes(x = Site,y=value, colour = Platform,shape = Ratio),alpha = 0.8,size = 6) + ylab("LODR") + colScale + xlab("Sequencing Site") + scale_y_log10()
print(pB)

###################
ratVarList <- ls(pattern=glob2rx(pattern=c("*RatVar")))
ratVarDatNames <-data.frame(DatName = ls(pattern=glob2rx(pattern=c("*RatVar"))))
ratVardat <- colsplit(ratVarDatNames$DatName,"\\.", c("Site","Analysis","Platform","Sample1","Sample2","Data"))

for (i in 1:length(ratVarList)){
  ratVardat$MinSDEst[i] = get(ratVarList[i])[1,"Estimate"]
  ratVardat$MinStErr[i] = get(ratVarList[i])[1,"Std. Error"]
  ratVardat$MaxSDEst[i] = get(ratVarList[i])[2,"Estimate"]
  ratVardat$MaxStErr[i] = get(ratVarList[i])[2, "Std. Error"]
  ratVardat$Lambda[i] = get(ratVarList[i])[3, "Estimate"]
  ratVardat$LambdaStErr[i]= get(ratVarList[i])[3, "Std. Error"]
}
ratVardat$Site <- as.factor(ratVardat$Site)
labelSites <- data.frame(Original = levels(ratVardat$Site), New = c("Lab 1","Lab 2","Lab 3","Lab 4","Lab 5","Lab 6","Lab 7","Lab 8","Lab 9"))

ratVardat$Site <- labelSites$New[match(ratVardat$Site,labelSites$Original)] 

ratVardat$Site <- as.factor(as.character(ratVardat$Site))

ratVardat_l = melt(ratVardat)
# colnames(ratVardat_l)[which(colnames(ratVardat_l) == "variable")] <- "Ratio"
# levels(ratVardat_l$Ratio) <- c("4:1","1:2","1:1.5")

theme_update(legend.justification=c(1,1), legend.position=c(1,1))
pD <- ggplot(ratVardat, aes(x = MinSDEst, y = MaxSDEst)) + geom_point(aes( colour = Platform),size = 7)+ geom_text(aes(label = Site), vjust = -1.5 )  + ylab(expression(V[max])) + xlab(expression(V[min])) + geom_errorbar(aes(ymin = MaxSDEst - MaxStErr, ymax = MaxSDEst + MaxStErr), colour = "grey30")+ geom_errorbarh(aes(xmin = MinSDEst - MinStErr, xmax = MinSDEst + MinStErr), colour = "grey30") + colScale #+ coord_cartesian(xlim=c(0,0.8),ylim= c(0,0.8))
print(pD)
LIFlib5var <- as.data.frame(LIFmaPlotAB$modRatVar)
names(LIFlib5var)[2]<- "Error"
ILMlib5var <- as.data.frame(ILMmaPlotAB$modRatVar)
names(ILMlib5var)[2]<- "Error"
pDOpen = pD + geom_point(data = LIFlib5var, aes(x = Estimate[1], y = Estimate[2]), colour = "dodgerblue3" ,size = 7, shape = 21) + geom_point(data = ILMlib5var, aes(x = Estimate[1], y = Estimate[2]), colour = "firebrick2" , size = 7, shape = 21) #+ geom_errorbar(data = ILMlib5var,aes(x = Estimate[1], ymin = Estimate[2] - Error[2], ymax = Estimate[2] + Error[2]), colour = "grey30")
# +
#   geom_errorbarh(LIFlib5var, aes(x = Estimate[1], y = Estimate[2], xmin = Estimate[1] - Error[1], xmax = Estimate[1] - Error[1]), colour = "grey30")+
#   geom_errorbar(ILMlib5var, aes(x = Estimate[1], y = Estimate[2], ymin = Estimate[2] - Error[2], ymax = Estimate[2] + Error[2]), colour = "grey30")+
#   geom_errorbarh(ILMlib5var,aes(x = Estimate[1], y = Estimate[2],xmin = Estimate[1] - Error[1], xmax = Estimate[1] - Error[1]), colour = "grey30")

pdOpenerr = pDOpen + geom_errorbar(data = ILMlib5var,aes(x = Estimate[1], ymin = Estimate[2] - Error[2], ymax = Estimate[2] + Error[2]), colour = "grey30")
#dev.off()
#multiplot(p1,p2,p,p3,cols=2)