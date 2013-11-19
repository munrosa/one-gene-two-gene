erccROC <- function(expDat){
  # get the req'd libraries
  require( ggplot2 )
  require( reshape2 )
  
  filenameRoot = expDat$sampleInfo$filenameRoot
  folds = expDat$sampleInfo$FCcode
  legendLabels = expDat$sampleInfo$legendLabels
  idCols = expDat$sampleInfo$idColsSRM
  idCols <- idCols[-which(is.na(idCols$Ratio)),]
  #Create a custom color scale
  #myColors <- c("#CC3333", "#339900","#FF9933","#66CC99")
  myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
  names(myColors) <- levels(folds$Ratio)
  myColorsDiff <- myColors[-(which(folds$FC == 1))]
  legendLabelsDiff <- legendLabels[-(which(folds$FC == 1))]
  
  
  
  colScale <- scale_colour_manual(name = "Ratio",values = myColorsDiff, labels = legendLabelsDiff)
  fillScale <- scale_fill_manual(name = "Ratio", values = myColorsDiff, labels = legendLabelsDiff)
  
  # get these handy libraries, too
  library("ROCR")
  library("pROC")
  
  # Read in the p.values from the file
  pValDat = read.csv(file=paste(filenameRoot,"ERCC","Pvals.csv"),header=T)
  names(pValDat)[1]= "Feature"
  pValDat <- pValDat[-2]
  
  pValDat$Feature <- as.character(pValDat$Feature)
  
  # now build the ROCR prediction objects
  # Format of FCcode = data.frame(Ratio = c("a","b","c","d"),FC =  c(4,1,.667,.5));
  FCcodeC = folds[-c(which(folds$FC == 1)),]
  pool.pred <- NULL
  FPR <- NULL
  TPR <- NULL
  FoldChange <- NULL
  ROCdat <- NULL
  AUCdat <-NULL
  for (i in 1:nrow(FCcodeC)){
    
    pool.predMeas = prediction(1- pValDat$Pval[ (pValDat$Fold == FCcodeC$FC[i]) | (pValDat$Fold == 1 )], pValDat$Fold[(pValDat$Fold == FCcodeC$FC[i]) | (pValDat$Fold == 1)],label.ordering = c(1, FCcodeC$FC[i]))
    pool.perf = performance(pool.predMeas, "tpr","fpr")
    pool.auc = performance(pool.predMeas, "auc")
    
    rocProc = roc(response=pValDat$Fold[(pValDat$Fold == FCcodeC$FC[i]) | (pValDat$Fold == 1)],predictor=1- pValDat$Pval[ (pValDat$Fold == FCcodeC$FC[i]) | (pValDat$Fold == 1 )],)
    
    # now build the three vectors for plotting - TPR, FPR, and FoldChange
    AUC = unlist(pool.auc@y.values)
    #print(paste("Ratio",as.character(FCcodeC$FC[i])))
    #print(paste("AUC",as.character(AUC)))
    AUCdatnew = data.frame(Ratio = legendLabelsDiff[i],AUC = round(AUC,digits = 3), Measured = length(pValDat$Fold[(pValDat$Fold == FCcodeC$FC[i])]), Spiked = length(idCols$Ratio[idCols$Ratio == FCcodeC$Ratio[i]]))
    AUCdat = rbind(AUCdat,AUCdatnew)
    FPR = c( unlist(pool.perf@x.values)) 
    TPR = c( unlist(pool.perf@y.values))
    #print(TPR)
    #print(FPR)
    Ratio = c(rep(as.character(FCcodeC$Ratio[i]), length(unlist(pool.perf@y.values))))
    ROCdatnew = data.frame(FPR = FPR, TPR = TPR, Ratio = Ratio)
    ROCdat = rbind(ROCdat,ROCdatnew)
  }
  
  AUCAnnot <- AUCdat
  cat("\nArea Under the Curve (AUC) Results:\n")
  print(AUCAnnot)
  
  AUCdat$xval = 0.7
  AUCdat$yval = seq(to = 0.25,from = 0.1,length.out=nrow(FCcodeC))
  
  ROCplot = ggplot(data = ROCdat, aes(x = FPR, y = TPR)) + geom_path(size = 2, aes(colour = Ratio), alpha = 0.7)+geom_point(size = 5, aes(colour = Ratio), alpha = 0.7) +colScale + geom_abline(intercept = 0, slope = 1, linetype = 2)+ annotation_custom(grob=tableGrob(AUCAnnot,show.rownames=F,equal.width=T,equal.height=T),xmin=0.375,xmax=1.0,ymin = 0,ymax = 0.25) + theme(legend.position=c(1,0.5)) #+geom_text(data = AUCdat,aes(x = xval, y = yval,label = paste(Ratio, "AUC =",as.character(round(AUC,digits=3))),hjust = 0)) 
  
  expDat$Figures$ROCplot <- ROCplot
  expDat$AUCdat <- AUCAnnot
  return(expDat)
 
}