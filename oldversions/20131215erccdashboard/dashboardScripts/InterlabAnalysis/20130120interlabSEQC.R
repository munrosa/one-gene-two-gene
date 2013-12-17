require(reshape2)
require(ggplot2)
source("R/multiplot.R")
theme_set(theme_bw(base_size=16))
#theme_update(legend.justification=c(0,1), legend.position=c(0,1))

pdf(file = paste("InterlabFig5pieces","pdf",sep="."),onefile=T,width=7,height = 7)

#Create a custom color scale
##myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
 myColors <- c("#6699CC", "#CC6666")
 names(myColors) <- c("ILM","LIF")
 colScale <- scale_colour_manual(name = "Platform",values = myColors)
 fillScale <- scale_fill_manual(name = "Platform", values = myColors)

# #FCcode = data.frame(Ratio = c("4:1","1:1","1:1.5","1:2"),FC =  c(4,1,.667,.5))
# legendLabels = c("4:1","1:1","1:1.5","1:2")
# 
# myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
# names(myColors) <- FCcode$Ratio
# myColorsDiff <- myColors[-(which(FCcode$FC == 1))]
# legendLabelsDiff <- legendLabels[-(which(FCcode$FC == 1))]
# colScale <- scale_colour_manual(name = "Ratio",values = myColorsDiff)
# fillScale <- scale_fill_manual(name = "Ratio", values = myColorsDiff)

load("MAY.NCTR.ILM.A.B.RData")
load("SQW.NIST.LIF.A.B.RData")
load("PSU.NIST.LIF.A.B.RData")
load("NWU.NIST.LIF.A.B.RData")
load("NVS.NCTR.ILM.A.B.RData")
load("AGR.NCTR.ILM.A.B.RData")
load("COH.NCTR.ILM.A.B.RData")
load("CNL.NCTR.ILM.A.B.RData")
load("BGI.NCTR.ILM.A.B.RData")
r_m_dat = data.frame(Site = c("BGI","CNL","MAY","AGR","COH","NVS"), Platform = "ILM")

r_m_dat2 = data.frame(Site = c("PSU","NWU","SQW"), Platform = "LIF")
r_m_dat = rbind(r_m_dat,r_m_dat2)

r_m_dat$Platform = as.factor(r_m_dat$Platform)
r_m_dat$r_m = c(BGI.NCTR.ILM.A.B.r_m,CNL.NCTR.ILM.A.B.r_m,MAY.NCTR.ILM.A.B.r_m,AGR.NCTR.ILM.A.B.r_m,COH.NCTR.ILM.A.B.r_m,NVS.NCTR.ILM.A.B.r_m,PSU.NIST.LIF.A.B.r_m,NWU.NIST.LIF.A.B.r_m,SQW.NIST.LIF.A.B.r_m)

r_m_dat$r_m_lower = c(BGI.NCTR.ILM.A.B.r_m_lower,CNL.NCTR.ILM.A.B.r_m_lower,MAY.NCTR.ILM.A.B.r_m_lower,AGR.NCTR.ILM.A.B.r_m_lower,COH.NCTR.ILM.A.B.r_m_lower,NVS.NCTR.ILM.A.B.r_m_lower,PSU.NIST.LIF.A.B.r_m_lower,NWU.NIST.LIF.A.B.r_m_lower,SQW.NIST.LIF.A.B.r_m_lower)

r_m_dat$r_m_upper = c(BGI.NCTR.ILM.A.B.r_m_upper,CNL.NCTR.ILM.A.B.r_m_upper,MAY.NCTR.ILM.A.B.r_m_upper,AGR.NCTR.ILM.A.B.r_m_upper,COH.NCTR.ILM.A.B.r_m_upper,NVS.NCTR.ILM.A.B.r_m_upper,PSU.NIST.LIF.A.B.r_m_upper,NWU.NIST.LIF.A.B.r_m_upper,SQW.NIST.LIF.A.B.r_m_upper)

p = ggplot(r_m_dat)+geom_errorbar(aes(x = Site,y=exp(r_m),ymin = exp(r_m_lower),ymax = exp(r_m_upper),colour = Platform), show_guide = F) + geom_hline(yintercept = 1.433,linetype = 2)+colScale + ylab(expression(r[m])) #+ geom_point(aes(x=Site, y = exp(r_m)), shape = 4) #+ geom_text(data = NULL,x = 7,y = (log(1.433)-0.05), label = expression(log(r[m])),parse=T))
print(p)
#p3 = ggplot(SDdat)+geom_point(aes(x = Site,y=ratioSD,colour = Platform),size = 6,shape = 4) + ylab("Overall Site Ratio SD")+colScale

AUCdat = r_m_dat[c(1,2)]

  AUCdat$x4 = rbind(BGI.NCTR.ILM.A.B.AUC$AUC[which(BGI.NCTR.ILM.A.B.AUC$RatioLab == "4:1")],CNL.NCTR.ILM.A.B.AUC$AUC[which(CNL.NCTR.ILM.A.B.AUC$RatioLab == "4:1")],MAY.NCTR.ILM.A.B.AUC$AUC[which(MAY.NCTR.ILM.A.B.AUC$RatioLab == "4:1")],AGR.NCTR.ILM.A.B.AUC$AUC[which(AGR.NCTR.ILM.A.B.AUC$RatioLab == "4:1")],COH.NCTR.ILM.A.B.AUC$AUC[which(COH.NCTR.ILM.A.B.AUC$RatioLab == "4:1")],NVS.NCTR.ILM.A.B.AUC$AUC[which(NVS.NCTR.ILM.A.B.AUC$RatioLab == "4:1")],PSU.NIST.LIF.A.B.AUC$AUC[which(PSU.NIST.LIF.A.B.AUC$RatioLab == "4:1")],NWU.NIST.LIF.A.B.AUC$AUC[which(NWU.NIST.LIF.A.B.AUC$RatioLab == "4:1")],SQW.NIST.LIF.A.B.AUC$AUC[which(SQW.NIST.LIF.A.B.AUC$RatioLab == "4:1")])

AUCdat$x2 = rbind(BGI.NCTR.ILM.A.B.AUC$AUC[which(BGI.NCTR.ILM.A.B.AUC$RatioLab == "1:2")],CNL.NCTR.ILM.A.B.AUC$AUC[which(CNL.NCTR.ILM.A.B.AUC$RatioLab == "1:2")],MAY.NCTR.ILM.A.B.AUC$AUC[which(MAY.NCTR.ILM.A.B.AUC$RatioLab == "1:2")],AGR.NCTR.ILM.A.B.AUC$AUC[which(AGR.NCTR.ILM.A.B.AUC$RatioLab == "1:2")],COH.NCTR.ILM.A.B.AUC$AUC[which(COH.NCTR.ILM.A.B.AUC$RatioLab == "1:2")],NVS.NCTR.ILM.A.B.AUC$AUC[which(NVS.NCTR.ILM.A.B.AUC$RatioLab == "1:2")],PSU.NIST.LIF.A.B.AUC$AUC[which(PSU.NIST.LIF.A.B.AUC$RatioLab == "1:2")],NWU.NIST.LIF.A.B.AUC$AUC[which(NWU.NIST.LIF.A.B.AUC$RatioLab == "1:2")],SQW.NIST.LIF.A.B.AUC$AUC[which(SQW.NIST.LIF.A.B.AUC$RatioLab == "1:2")])

AUCdat$x1.5 = rbind(BGI.NCTR.ILM.A.B.AUC$AUC[which(BGI.NCTR.ILM.A.B.AUC$RatioLab == "1:1.5")],CNL.NCTR.ILM.A.B.AUC$AUC[which(CNL.NCTR.ILM.A.B.AUC$RatioLab == "1:1.5")],MAY.NCTR.ILM.A.B.AUC$AUC[which(MAY.NCTR.ILM.A.B.AUC$RatioLab == "1:1.5")],AGR.NCTR.ILM.A.B.AUC$AUC[which(AGR.NCTR.ILM.A.B.AUC$RatioLab == "1:1.5")],COH.NCTR.ILM.A.B.AUC$AUC[which(COH.NCTR.ILM.A.B.AUC$RatioLab == "1:1.5")],NVS.NCTR.ILM.A.B.AUC$AUC[which(NVS.NCTR.ILM.A.B.AUC$RatioLab == "1:1.5")],PSU.NIST.LIF.A.B.AUC$AUC[which(PSU.NIST.LIF.A.B.AUC$RatioLab == "1:1.5")],NWU.NIST.LIF.A.B.AUC$AUC[which(NWU.NIST.LIF.A.B.AUC$RatioLab == "1:1.5")],SQW.NIST.LIF.A.B.AUC$AUC[which(SQW.NIST.LIF.A.B.AUC$RatioLab == "1:1.5")])



AUCdat_l = melt(AUCdat)
#names(AUCdat_l)[which(names(AUCdat_l) == "variable")] <- "Ratio"
p1 = ggplot(AUCdat_l,show_guide = F)+geom_point(aes(x = Site,y=value,colour = Platform), show_guide = F,shape = 1,size = 10) +geom_text(aes(x = Site,y=value,label = gsub(pattern="x",replacement="",variable), colour = Platform),size = 4,show_guide = F) + ylab("AUC") + colScale

#p1 = ggplot(AUCdat_l)+geom_bar(aes(x=Site,y=value,colour = variable, fill = Platform), position = "dodge",stat = "identity") + ylab("AUC") #+ colScale

print(p1)

LODRdat = r_m_dat[c(1,2)]

  LODRdat$x4 = rbind(BGI.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(BGI.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "2")],CNL.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(CNL.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "2")],MAY.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(MAY.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "2")],AGR.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(AGR.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "2")],COH.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(COH.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "2")],NVS.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(NVS.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "2")],PSU.NIST.LIF.A.B.lodr.ERCC$Estimate[which(PSU.NIST.LIF.A.B.lodr.ERCC[c(1)] == "2")],NWU.NIST.LIF.A.B.lodr.ERCC$Estimate[which(NWU.NIST.LIF.A.B.lodr.ERCC[c(1)] == "2")],SQW.NIST.LIF.A.B.lodr.ERCC$Estimate[which(SQW.NIST.LIF.A.B.lodr.ERCC[c(1)] == "2")])

LODRdat$x2 = rbind(BGI.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(BGI.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "1")],CNL.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(CNL.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "1")],MAY.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(MAY.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "1")],AGR.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(AGR.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "1")],COH.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(COH.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "1")],NVS.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(NVS.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "1")],PSU.NIST.LIF.A.B.lodr.ERCC$Estimate[which(PSU.NIST.LIF.A.B.lodr.ERCC[c(1)] == "1")],NWU.NIST.LIF.A.B.lodr.ERCC$Estimate[which(NWU.NIST.LIF.A.B.lodr.ERCC[c(1)] == "1")],SQW.NIST.LIF.A.B.lodr.ERCC$Estimate[which(SQW.NIST.LIF.A.B.lodr.ERCC[c(1)] == "1")])

LODRdat$x1.5 = rbind(BGI.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(BGI.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "0.584")],CNL.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(CNL.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "0.584")],MAY.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(MAY.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "0.584")],AGR.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(AGR.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "0.584")],COH.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(COH.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "0.584")],NVS.NCTR.ILM.A.B.lodr.ERCC$Estimate[which(NVS.NCTR.ILM.A.B.lodr.ERCC[c(1)] == "0.584")],PSU.NIST.LIF.A.B.lodr.ERCC$Estimate[which(PSU.NIST.LIF.A.B.lodr.ERCC[c(1)] == "0.584")],NWU.NIST.LIF.A.B.lodr.ERCC$Estimate[which(NWU.NIST.LIF.A.B.lodr.ERCC[c(1)] == "0.584")],SQW.NIST.LIF.A.B.lodr.ERCC$Estimate[which(SQW.NIST.LIF.A.B.lodr.ERCC[c(1)] == "0.584")])

LODRdat_l = melt(LODRdat)
#names(LODRdat_l)[which(names(LODRdat_l) == "variable")] <- "Ratio"

p2 = ggplot(LODRdat_l, show_guide = F)+geom_point(aes(x=Site,y=value,colour = Platform),shape = 1,size = 10, show_guide = F)+geom_text(aes(x = Site,y=value,label = gsub(pattern="x",replacement="",variable), colour = Platform),size = 4,show_guide = F) + ylab("LODR") +colScale

#p2 = ggplot(LODRdat_l)+geom_bar(aes(x=Site,y=value,colour = Ratio)) + ylab("LODR") + colScale

print(p2)

SDdat = r_m_dat[c(1,2)]

SDdat$ratioSD = c(BGI.NCTR.ILM.A.B.sdGlobal,CNL.NCTR.ILM.A.B.sdGlobal,MAY.NCTR.ILM.A.B.sdGlobal,AGR.NCTR.ILM.A.B.sdGlobal,COH.NCTR.ILM.A.B.sdGlobal,NVS.NCTR.ILM.A.B.sdGlobal,PSU.NIST.LIF.A.B.sdGlobal,NWU.NIST.LIF.A.B.sdGlobal,SQW.NIST.LIF.A.B.sdGlobal)

p3 = ggplot(SDdat, show_guide = F) +geom_point(aes(x = Site,y=ratioSD,colour = Platform),size = 6,shape = 4, show_guide = F) + ylab("Overall Site Ratio SD")+colScale
print(p3)

dev.off()
#multiplot(p1,p2,p,p3,cols=2)