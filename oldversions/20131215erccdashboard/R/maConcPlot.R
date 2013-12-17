# Plots with target ratios and R adjusted ratios, MA plots and Ratio Summaries
maConcPlot <-function(expDat, LODR.annot.ERCC, alphaPoint = 0.8,
                      r_mAdjust = T, replicate = T){
  
  ReplicateName = "Rep"
  # Melt the data
  expDat <- meltExpDat(expDat, cnt = expDat$Transcripts, 
                       designMat = expDat$designMat)
  
  countPair <- expDat$expressDat_l
  sampleInfo <- expDat$sampleInfo
  myYLim = sampleInfo$myYLimMA
  myXLim = sampleInfo$myXLim
  siteName = sampleInfo$siteName
  analysis = sampleInfo$analysis
  filenameRoot = sampleInfo$filenameRoot
  sample1 = expDat$sampleNames[1]
  sample2 = expDat$sampleNames[2]
 
  
  idCols = expDat$idColsAdj
  r_m.res = expDat$r_m.res
  
  cutoffs = LODR.annot.ERCC$cutoffs
  FCcode = sampleInfo$FCcode
  legendLabels=sampleInfo$legendLabels
  avexlabel = expDat$ERCCxlabelAve
  spikeFraction = expDat$spikeFraction
  
  r_m.mn = exp(r_m.res$r_m.mn)
  r_m.UL = exp(r_m.res$r_m.upper)
  r_m.LL = exp(r_m.res$r_m.lower)
  
  theme_update(legend.justification=c(1,0), legend.position=c(1,0))
  
  #Create a custom color scale
  myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
  names(myColors) <- levels(FCcode$Ratio)
  colScale <- scale_colour_manual(name = "Ratio",values = myColors, 
                                  labels = legendLabels)
  fillScale <- scale_fill_manual(name = "Ratio", values = myColors, 
                                 labels = legendLabels)
  library(gridExtra)
  
  if (length(cutoffs) > 0){
    cat("\nLODR estimates are available to code ratio-abundance plot\n")  
  }
  
  counts1 = countPair$NormCounts[which(countPair$Sample == sample1)]
  counts2 = countPair$NormCounts[which(countPair$Sample == sample2)]
  
  if (length(counts1) != length(counts2)){
    stop("Uneven replication for the sample pair")
  }
 
  
  if (replicate == T){
    Mdat = log2(counts1) - log2(counts2)
    Adat = log2((counts1 + counts2)/2)
    colIdx <- which(colnames(countPair) == ReplicateName)
    maData = data.frame(Feature = countPair$Feature[which(countPair$Sample == 
                                                            sample1)], 
                        Ratio = countPair$Ratio[which(countPair$Sample ==
                                                        sample1)],
                        M = Mdat, A = Adat, Replicate = 
                          countPair[which(countPair$Sample == sample1),colIdx])
  }
  else{
    stop("Replicate samples are needed")
    maData = data.frame(Feature = countPair$Feature, Ratio = countPair$Ratio,
                        M = countPair$Mdat, A = countPair$Adat)
  }
    
  names(maData)[3:4] = c("M","A")
  maData$Nominal = FCcode$FC[1]
  maData$Nominal[which(maData$Ratio == FCcode$Ratio[2])] = FCcode$FC[2]
  maData$Nominal[which(maData$Ratio == FCcode$Ratio[3])] = FCcode$FC[3]
  maData$Nominal[which(maData$Ratio == FCcode$Ratio[4])] = FCcode$FC[4]
  
  if(r_mAdjust == T){
    maData$Empirical = maData$Nominal/r_m.mn  
  }
  xlabel = xlab(avexlabel)
  if(replicate == T){
   # Use Average Concentration values (instead of Average Signals)
    FeatureA = data.frame(Feature= idCols$Feature, A = log2((idCols$Conc1 +
                                                               idCols$Conc2)/2))
   for (i in 1:nlevels(FeatureA$Feature)){
     maData$A[which(maData$Feature == as.character(FeatureA$Feature[i]))] =
       FeatureA$A[i]  
   }
    maDataAve <- data.frame(tapply(maData$M,maData[,"Feature"],mean))
    
    maDataAveSD = cbind(row.names(maDataAve),
                      maData$Ratio[match(row.names(maDataAve),
                                         table=maData$Feature)],
                      maData$Empirical[match(row.names(maDataAve),
                                             table=maData$Feature)],
                      maData$A[match(row.names(maDataAve),
                                     table=maData$Feature)], maDataAve[c(1)])
  
  names(maDataAveSD)[1:5] = c("Feature","Ratio","Empirical","A","M.Ave")
  
  maDataSD <- as.vector(tapply(maData$M,maData[,"Feature"],sd))
  maDataAveSD$M.SD = maDataSD
  
  cutERCCs = unique(maDataAveSD$Feature[which(is.na(maDataAveSD$M.SD))])
    
  if(length(cutERCCs) != 0){
    cat(paste("\nThese ERCCs were not included in the ratio-abundance plot, \n",
              "because not enough non-zero replicate measurements of these \n",
              "controls were obtained for both samples:\n"))
    
    cat(paste(as.character(cutERCCs),collapse="\n"))
    
    maDataAveSD = maDataAveSD[-(which(maDataAveSD$Feature %in% cutERCCs)),]
    searchCut = paste(cutERCCs,collapse="|")
    maData = maData[-(grep(searchCut,maData$Feature)),]
  }
    
    maDataAveSD$Feature = as.factor(as.character(maDataAveSD$Feature))  
    maData$Feature = as.factor(as.character(maData$Feature))
    #write.csv(maData,file = paste(filenameRoot,"maDataFinite.csv",sep = "."))
    
    # Estimate SD for all ratio measurements
    # For the ERCCs that are plotted take the log ratio data and subtract the log 
    # nominal ratio, this will normalize the data
    normRats = maData$M - log2(maData$Nominal)
    #print(normRats)
    
    # Take the sd of all of the normalized log ratio data to find a global SD for the ERCCs at this site
    sdGlobal = sd(normRats)
    cat("\nGlobal Ratio SD for this sample pair is: ")
    cat(sdGlobal, "\n")
    
    ratVarDat <- maDataAveSD
    
    ratVarDat$A <- log2((2^(ratVarDat$A))/(spikeFraction))
    
    theme_update(legend.justification=c(1,1), legend.position=c(1,1))
    
    sdRatioplot <- ggplot() + geom_point(data = ratVarDat, aes(x = A, y = M.SD, 
                                                               colour = Ratio),
                                         size = 5, alpha = alphaPoint)+ colScale
    
    
    lmfit <- lm(formula = M.SD ~ 1, data = ratVarDat, weights = 1/((M.SD)^2))
    
    stdevCoef <- summary(lmfit)$coefficients
    
    params<-list(Mo = 0.1, lambda = 2) # setup default params.
    f = coef(lmfit)[[1]]  
    #print("Using default parameters")
    #print(params)
    #print(f)
    
    while(TRUE){
      
      expmod<-NULL
      try(expmod<-nls(formula=M.SD ~ Mo*exp(-lambda * A) + f, data = ratVarDat,
                      start = params)); # does not stop in the case of error
      
      if(!is.null(expmod))break; # if nls works, then quit from the loop
      #print("Try changing initial guess for parameters")
      params<-list(Mo = sdGlobal, lambda = 0.2)# change the params for nls
      #print(params)
      
    }
    
    if(!is.null(expmod)){
      stdevCoef <- rbind(stdevCoef, summary(expmod)$coefficients)
    } 
    
    #print(stdevCoef)
    row.names(stdevCoef)<- c("Minimum SD Estimate", "Maximum SD Estimate",
                             "Lambda")
    print(stdevCoef[,c(1:2)])
    
    
    # add fitted curve
    lineDat <- data.frame("xdat" = ratVarDat$A,"ydat" = 
                            predict(lmfit, x = ratVarDat$A)) 
    if(!is.null(expmod)){
      fitDat <- data.frame("xdat" = ratVarDat$A,"ydat" = predict(expmod,
                                                                 list(x = ratVarDat$A))) 
    } 
    sdRatioplotFit <- sdRatioplot + geom_line(data = lineDat,aes(x=xdat,y=ydat))+
      ylab("sd(Log2 Ratios of Counts)") + 
      xlab(expression(paste("Log2 Average ERCC Abundance in Pure Spike Mixes ", 
                            "(attomole nt/",mu,"L Mix)",sep = "")))
    
    if(!is.null(expmod)){
      sdRatioplotFit <- sdRatioplotFit + 
        geom_line(data = fitDat, aes(x = xdat, y = ydat))
    } 
    
    theme_update(legend.justification=c(1,0), legend.position=c(1,0))
    if (length(cutoffs)>0){
      maDataAveSD$LODR = "below"
      FCcodeC = FCcode[-c(which(FCcode$FC == 1)),]
      for (i in 1:length(cutoffs)){
        maDataAveSD$LODR[which((maDataAveSD$A > cutoffs[i])&
                                 (maDataAveSD$Ratio ==
                                    FCcodeC$Ratio[i]))] = "above"
      }
      
      maDataAveSD$LODR <- as.factor(maDataAveSD$LODR)
      rm_dat = data.frame("Mean r_m" = signif(r_m.mn,digits = 4),
                          "Lower Limit" = signif(r_m.LL,digits = 4), 
                          "Upper Limit" = signif(r_m.UL,digits = 4))
      colnames(rm_dat) <- c("Mean", expression(paste("95% CI Lower Bound")), 
                            expression(paste("95% CI Upper Bound")))
      rownames(rm_dat) <- expression(r[m])
      
      # Plot ratio-signal data coding points below the LODR with open circles
      maPlot <- ggplot(maDataAveSD, aes(x = A, y = M.Ave) ) + 
        geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD, 
                          colour = Ratio),alpha = 0.3) + 
        geom_point(aes(colour = Ratio),size = 5, alpha = alphaPoint) +
        geom_point(data = subset(maDataAveSD, (LODR == "below")),
                   colour = "white",size = 2.5) + 
        geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
                   linetype = "longdash",alpha = 0.7) + 
        ylab("Log2 Ratio of Counts") + xlabel + 
        coord_cartesian(xlim = myXLim, ylim = myYLim) + colScale + 
        annotation_custom(tableGrob(rm_dat,parse=T), xmin = 
                            0.25*max(maDataAveSD$A), xmax = max(maDataAveSD$A),
                          ymin = (myYLim[2]) - 0.25*myYLim[2],  ymax = myYLim[2] ) 
    }else{
      # Plot ratio-signal data without LODR coding
      maPlot <- ggplot(maDataAveSD, aes(x = A, y = M.Ave, colour = Ratio) ) + 
        geom_errorbar(aes(ymax = M.Ave + M.SD, ymin = M.Ave - M.SD), 
                      alpha = 0.3) + geom_point(size = 5, alpha = alphaPoint) + 
        geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), 
                   linetype = "longdash",alpha = 0.7) + 
        ylab("Log2 Ratio of Counts") + xlabel + colScale+ 
        coord_cartesian(xlim = myXLim, ylim = myYLim)
    }
  }else{
    
    maDataSub = subset(maData, is.finite(maData$A))
    
    # Use average concentration values instead of average signals
    idCols = idCols[match(countPair$Feature,idCols$Feature),]
    A = log2((idCols$Conc1 + idCols$Conc2)/2)
    maData$A = A
    maDataSub = subset(maData, is.finite(maData$A))
    
    maPlot <- ggplot(maDataSub, aes(x = A, y = M)) + geom_point(aes(colour = Ratio), alpha = alphaPoint, size = 5) + geom_point(data = subset(maDataSub, (Ratio == FCcode$Ratio[2])),colour = "white", size = 3) +geom_point(data = subset(maDataSub, (A < cutoffs[1])&(Ratio == FCcode$Ratio[1])),colour = "white", size = 3) + geom_point(data = subset(maDataSub, (A < cutoffs[2])&(Ratio == FCcode$Ratio[3])),colour = "white", size = 3) + geom_point(data = subset(maDataSub, (A < cutoffs[3])&(Ratio == FCcode$Ratio[4])),colour = "white", size = 3) + geom_hline(aes(yintercept = log2(Empirical), colour = Ratio), linetype = "longdash",alpha = 0.7) + ylab("Log2 Ratio of Counts") + xlabel + scale_color_hue(l = 40)+ coord_cartesian(xlim = myXLim, ylim = myYLim)
}
  
  #print(maPlot)
  #print(sdRatioplotFit)
  
  expDat$Figures$plotRatioAnnot <- maPlot
  expDat$sdGlobal <- sdGlobal
  expDat$modRatVar <- stdevCoef
  expDat$Figures$plotSDratio <- sdRatioplotFit
  return(expDat)
  #return(list(maPlot = maPlot,sdGlobal = sdGlobal, modRatVar = stdevCoef, sdRatioplotFit = sdRatioplotFit))
}