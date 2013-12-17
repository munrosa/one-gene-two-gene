# Compare percent mapped to ERCCs for ILM and LIF in the SEQC study
library(ggplot2)
library(reshape2)
theme_set(theme_bw(base_size=16))


#pdf(file = paste("InterlabFig5pieces","pdf",sep="."),onefile=T,width=7,height = 7)

#Create a custom color scale
#myColors <- c("#6699CC", "#CC6666")
myColors <- c("firebrick2", "dodgerblue3")  
names(myColors) <- c("ILM","LIF")
colScale <- scale_colour_manual(name = "Platform",values = myColors)
fillScale <- scale_fill_manual(name = "Platform", values = myColors)


files <- list.files("data/MainSEQC/all_qc_results/")

qcResults <- read.delim(paste("data/MainSEQC/all_qc_results/",files[1], sep = ""))
qcResults$readType <- NULL

for (i in 2:length(files)){
  print(i)
  addResults <- read.delim(paste("data/MainSEQC/all_qc_results/",files[i], sep = ""))
  print(head(addResults,10))
  addResults$readType <- NULL
  qcResults <- rbind(qcResults, addResults) 
}
str(qcResults)


labelSites <- data.frame(Original = levels(qcResults$siteCode), New = c("Lab 1","Lab 2","Lab 3","Lab 4","Lab 5","Lab 6","Lab 7","Lab 8","Lab 9"))

qcResults$siteCode <- labelSites$New[match(qcResults$siteCode,labelSites$Original)] 

qcResults$siteCode <- as.factor(as.character(qcResults$siteCode))


qcResults$repNumber <- as.factor(qcResults$repNumber)
qcResults$percMappedERCC <- qcResults$numMappedERCC/qcResults$numReads
ggplot(subset(qcResults, sampleType != "E" & sampleType !="F"))+ geom_boxplot(aes(x =repNumber,y=percMappedERCC, colour = platform))

ggplot(subset(qcResults, sampleType != "E" & sampleType !="F"))+ geom_boxplot(aes(x =siteCode,y=percMappedERCC, colour = platform))+ facet_grid(~repNumber)

ggplot(subset(qcResults, sampleType == "A" | sampleType =="B"))+ geom_boxplot(aes(x =siteCode,y=percMappedERCC, colour = platform))+ facet_grid(~repNumber)+colScale

ggplot(subset(qcResults, sampleType == "A" | sampleType =="B"))+ geom_boxplot(aes(x =repNumber,y=percMappedERCC, colour = platform))+ facet_grid(~siteCode)+colScale
