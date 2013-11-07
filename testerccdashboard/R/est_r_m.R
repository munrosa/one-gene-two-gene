est_r_m <- function(expDat, cnt, printPlot = F){
  print("Check for sample mRNA fraction differences(r_m)...")
  sampleInfo = expDat$sampleInfo
  site <- sampleInfo$siteName
  avexlabel <- expDat$ERCCxlabelAve
  myXLim <- sampleInfo$myXLim
  idCols <- expDat$idCols
  # requires that the sample1 columns are first in the table -> need to add a stopifnot statement for this 
  theme_update(legend.justification=c(1,0), legend.position=c(1,0))
  
  #Create a custom color scale
  
  myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
  names(myColors) <- levels(sampleInfo$FCcode$Ratio)
  colScale <- scale_colour_manual(name = "Ratio",values = myColors, labels = sampleInfo$legendLabels)
  fillScale <- scale_fill_manual(name = "Ratio", values = myColors, labels = sampleInfo$legendLabels)
  dat = unique(cnt)
  Features = make.names(dat$Feature,unique=T)
  Features = gsub(".","-", Features, fixed = T)
  rownames(dat)<-Features; dat<-as.matrix(dat[,-1])
    
  colnames(dat)<-paste(rep(c(expDat$sample1,expDat$sample2),each=ncol(dat)/2),c(1:(ncol(dat)/2),1:(ncol(dat)/2)),sep="")
  
  ## Get ERCC names
  ERCC<-rownames(dat[substr(rownames(dat),1,5)=="ERCC-",])

  ## Specify Sample (A or B)
  trt<-rep(1:2,each=ncol(dat)/2)
  design.list<-list(trt,rep(1,ncol(dat)))
  
  ## Compute offset (e.g. total counts, 75% quantile, TMM, etc) could modify
  # in future to enable other methods
  log.offset<-log(colSums(dat))

  if(sampleInfo$totalSeqReads == T){
      log.offset = log(expDat$totalReads)
  }
  
  ######################################################
  #### Estimate r_m from each ERCC using NegBin GLM ####
  ######################################################
  
  ERCC.FC = idCols[c(1,4)];rownames(ERCC.FC)<-ERCC.FC[,1]
  ERCC.FC$NomRatio <- 1
  for (i in 1:nlevels(sampleInfo$FCcode$Ratio)){
    ERCC.FC$NomRatio[which(ERCC.FC$Ratio == as.character(
      sampleInfo$FCcode$Ratio[i]))] = sampleInfo$FCcode$FC[i]  
  }
  ERCC.Ratio = ERCC.FC[c(1,2)]
  
  ERCC.FC = ERCC.FC[-c(2)]
  
  library(MASS)
  print("log.offset")
  print(log.offset)
  
  r_m<-NULL
  print("Number of ERCC Controls Used in r_m estimate")
  print(length((1:nrow(dat))[substr(rownames(dat),1,5)=="ERCC-"]))
  #print(dat[(1:nrow(dat))[substr(rownames(dat),1,5)=="ERCC-"],])
  
  for( i in (1:nrow(dat))[substr(rownames(dat),1,5)=="ERCC-"]){ 
  
#### Obtain estimated log fold-change and corresponding standard error for each ERCC
    r_m<-rbind(r_m,summary(suppressWarnings(glm.nb(dat[i,]~as.factor(trt)+offset(log.offset))))$coefficients[2,1:2])
  }
  
colnames(r_m)<-c("r_m.hat","r_m.se"); rownames(r_m)<-rownames(dat)[substr(rownames(dat),1,5)=="ERCC-"]

  #### Add nominal log fold change to results
r_m <- data.frame(r_m)
  
r_m$nominal<- - log(ERCC.FC[rownames(r_m),2])
  
#### Make 95% CIs based on t-distribution
quant<-qt(.975,6)

#### Plot log r_m estimates and corresponding 95% CIs
  r_m$Ratio = ERCC.Ratio$Ratio[match(row.names(r_m), ERCC.Ratio$Feature)]
  # Subset data frame to remove missing ERCC controls
  r_m <- subset(r_m, is.finite(Ratio))
  # weighted mean of the individual means
  r_m.mn<-sum((r_m[,1]-r_m[,3])/r_m[,2]^2)/sum(r_m[,2]^-2)
  
  
  # standard deviation, sigma of the weighted mean
  r_m.mn.sd <- sqrt(1/(sum(1/(r_m[,2]^2))))
  
  r_m = cbind(as.character(row.names(r_m)), r_m)
  names(r_m)[1]= "Feature"
  
  
  r_m.mnlog = exp(r_m.mn)
  
  idCols$Conc2 = r_m.mnlog*idCols$Conc2
  idColsAve = idCols
  
  idColsAve = idCols[match(r_m$Feature,idCols$Feature),]
  
  r_m$AveConc = log2((idColsAve$Conc1+ idColsAve$Conc2)/2)
  
  r_m = r_m[order(r_m$AveConc),]
  
  
  r_m$ymax = r_m$r_m.hat - r_m$nominal + (quant)*r_m$r_m.se
  r_m$ymin = r_m$r_m.hat - r_m$nominal - (quant)*r_m$r_m.se
  textDat = subset(r_m, (r_m.mn<(r_m.hat - nominal - (quant*r_m.se)))|(r_m.mn>(r_m.hat - nominal + (quant*r_m.se))))
  
  print("Outlier ERCCs for GLM r_m Estimate:")
  print(as.character(textDat$Feature))
  
  print(paste(site,"GLM log(r_m) estimate:"))
  print(r_m.mn)
  print("GLM log(r_m) estimate standard deviation: ")
  print(r_m.mn.sd)
  print(nrow(r_m))
  r_m.upper.limit = r_m.mn + (quant)*((r_m.mn.sd)/(sqrt(nrow(r_m))))
  r_m.lower.limit = r_m.mn - (quant)*((r_m.mn.sd)/(sqrt(nrow(r_m))))
  
  print(paste(site,"GLM r_m estimate:"))
  print(exp(r_m.mn))
  
  print("upper limit")
  print(exp(r_m.upper.limit))
  
  print("lower limit")
  print(exp(r_m.lower.limit))
  
  if (nrow(textDat)>1){
    plotSiter_m = ggplot(r_m, aes(x = AveConc, y = r_m.hat - nominal, colour =  Ratio)) + geom_point(size = 6, alpha = 0.7) + geom_errorbar(aes(ymin = ymin,ymax = ymax),alpha = 0.7) +  coord_cartesian(xlim = myXLim, ylim = c(-1.5,1.5)) + xlab(avexlabel) + ylab(expression(log(r[m]))) + geom_text(data = textDat, aes(x = AveConc, y = (r_m.hat - nominal), label = gsub("ERCC-00","",Feature)),colour = "black", size = 6,show_guide = F,angle = 90) + geom_hline(yintercept = r_m.mn) + colScale  
  }else{
    plotSiter_m = ggplot(r_m, aes(x = AveConc, y = r_m.hat - nominal, colour =  Ratio)) + geom_point(size = 6, alpha = 0.7) + geom_errorbar(aes(ymin = ymin,ymax = ymax),alpha = 0.7) +  coord_cartesian(xlim = myXLim, ylim = c(-1.5,1.5)) + xlab(avexlabel) + ylab(expression(log(r[m]))) + geom_hline(yintercept = r_m.mn) + colScale
  }
  
  if(printPlot != F){
    print(plotSiter_m)  
  }
  
  expDat$idColsAdj <- idCols
  expDat$r_m.res <- list(r_m.mn = r_m.mn, r_m.upper = r_m.upper.limit,r_m.lower = r_m.lower.limit) 
  return(expDat)
}