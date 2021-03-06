estLODR <- function(expDat,kind = "ERCC", prob=0.9){
  cat(paste("\nEstimating Limit of Detection of Ratios (LODR) for ",
            kind, " data...\n"))
  ##############
  #### LODR ####
  ##############
  ### Select fitting parameters
  # pval.cutoff default is .001
  # probability cutoff is .9
  # kind can be ERCC or Sim
  sampleInfo <- expDat$sampleInfo
  pval.cutoff<- expDat$p.thresh 
  filenameRoot<-sampleInfo$filenameRoot
  FCcode <- sampleInfo$FCcode
  #FCcode = data.frame(Ratio = c("4:1","1:1","1.5:1","2:1"),FC =  c(4,1,.667,.5))
  legendLabels = sampleInfo$legendLabels

  #Create a custom color scale
  myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
  names(myColors) <- levels(FCcode$Ratio)
  colScale <- scale_colour_manual(name = "Ratio",values = myColors, labels = legendLabels)
  fillScale <- scale_fill_manual(name = "Ratio", values = myColors, labels = legendLabels)
  
  ### Read in results for LODR ####
  ### Modify this file path to work with other ERCC data
  ### Note that you must format your provided data in the same manner as the provided example
  ### i.e. the column headings of your selected datamust match those of the provided examples
  pval.res<-read.csv(paste(filenameRoot,kind,"Pvals.csv"),header=TRUE) 
  names(pval.res)[1]= "Feature"
  
  #library(locfit)
  #library(gridExtra)
  #### Function used to estimate LODR and its uncertainty
  LODR<-function(pval,mn,cutoff,prob){
    cutoff<-log10(cutoff)
    fit<-locfit(log10(pval)~lp(log10(mn)))
    X<-preplot(fit,band="pred",newdata=log10(mn))
    
    ### See plot of fit
    #plot(fit,band="pred",get.data=TRUE,xlim=range(log10(mn)))
    
    find.mn<-function(mn,fit,cutoff,prob){
      X<-preplot(fit,newdata=mn,band="pred")
      (X$fit+qnorm(prob)*X$se.fit-cutoff)^2
    }
   
    rng.mn<-range(log10(mn))
    
  #  ### Search in sections to get first crossing 
   # segmented.search<-function(fit){
# 	ppp<-.2
#     t.lodr<-optimize(find.mn,c(rng.mn[1],sum(rng.mn*c(1-ppp,ppp))),fit=fit,cutoff=cutoff,prob=prob)
#     while(t.lodr$objective>.0001&ppp<=1){
#         t.lodr<-optimize(find.mn,c(sum(rng.mn*c(1-ppp+.2,ppp-.2)),sum(rng.mn*c(1-ppp,ppp))),fit=fit,cutoff=cutoff,prob=prob)
#         ppp<-ppp+.2
#     }
#     t.lodr
# }   
  ### Search in sections to get first crossing 
  segmented.search<-function(fit){
    X<-preplot(fit,newdata=min(log10(mn)),band="pred")
    if((X$fit+qnorm(prob)*X$se.fit)<cutoff){
      t.lodr<-list(minimum=min(log10(mn)),objective=0)
    } else{
      ppp<-.2
      t.lodr<-optimize(find.mn,c(rng.mn[1],sum(rng.mn*c(1-ppp,ppp))),fit=fit,cutoff=cutoff,prob=prob)
      while(t.lodr$objective>.0001&ppp<=1){
        t.lodr<-optimize(find.mn,c(sum(rng.mn*c(1-ppp+.2,ppp-.2)),sum(rng.mn*c(1-ppp,ppp))),fit=fit,cutoff=cutoff,prob=prob)
        ppp<-ppp+.2
      }
    }
    t.lodr
  }    
  
  lodr<-segmented.search(fit)
     
    ### Bootstrap to estimate uncertainty
    lodr.boot<-NULL
    for(ii in 1:500){
      y.boot<-X$fit+sample(residuals(fit),length(mn))
      #if(ii %in%c(20*(1:5))){points(log10(mn),y.boot,col=j); j<-j+1}
      fit.boot<-locfit(y.boot~lp(log10(mn)))
      lodr.boot<-c(lodr.boot,segmented.search(fit.boot)$minimum)
      if (ii %% 100 == 0){
        cat ("...")
      }
    }
    return(c(lodr$objective,lodr$minimum,quantile(lodr.boot,c(.05,.95))))
  }
    

  #### Apply to loaded data ####
  lodr.resPlot<-NULL; set.seed(1)
  lodr.resLess<-NULL; set.seed(1)
  
  #lodr.res<-NULL; set.seed(1)
  lineDat <-NULL
  
  pval.res$Ratio = "Ratio"
  for (i in 1:nlevels(FCcode$Ratio)){
    pval.res$Ratio[which(pval.res$Fold == FCcode$FC[i])] = as.character(levels(FCcode$Ratio)[i])
  }
  
  pval.res$Ratio <- as.factor(pval.res$Ratio)
  for(i in 1:nlevels(pval.res$Ratio)){  ### Analyze each fold change subgroup separately
    grp<-which(pval.res$Ratio==levels(pval.res$Ratio)[i])  ### identify subgroup members

    x<-pval.res$MnCnt[grp][pval.res$Pval[grp]!=0]; y<-pval.res$Pval[grp][pval.res$Pval[grp]!=0]
    
    fit<-locfit(log10(y)~lp(log10(x)))
    #x.new<-log10(pval.res$MnCnt[grp])
    x.new<-seq(min(log10(x)),max(log10(x)),length.out=100)
    X<-preplot(fit,band="pred",newdata=x.new)

    x.new<-10^x.new
    fitLine = 10^(X$fit)
    fitUpper = 10^(X$fit+qnorm(prob)*X$se.fit)
    fitLower = 10^(X$fit-qnorm(prob)*X$se.fit)
    fitData = data.frame(x.new,fitLine,fitUpper,fitLower,Ratio = levels(pval.res$Ratio)[i])
    lineDat = rbind(lineDat,fitData)

# ### Estimate LODR
#lodr.res<-rbind(lodr.res,c(abs(log2(FCcode$FC[i])),LODR(y,x,cutoff=pval.cutoff,prob=prob)))

    ### Estimate LODR
    if(FCcode$FC[i]!=1){
      t.res<-LODR(y,x,cutoff=pval.cutoff,prob=prob)
      t.res[-1]<-signif(10^t.res[-1],2)
      if(t.res[1]>.01){
        t.res[2]<-Inf
        t.res[3:4]<-NA
      }
      t.resLess <- t.res
      t.resLess[-1][t.resLess[-1]==signif(min(x),2)]<-paste("<",signif(min(x),2),sep="")
      t.res[-1][t.res[-1]==signif(min(x),2)] <- Inf
    }
    if(FCcode$FC[i]==1){
      t.res<-rep(NA,4)
      t.resLess<-rep(NA,4)
    } 
    lodr.resPlot<-rbind(lodr.resPlot,c(round(abs(log2(FCcode$FC[i])),3),t.res))
    lodr.resLess <- rbind(lodr.resLess,c(round(abs(log2(FCcode$FC[i])),3),t.resLess))
    
}
  pval.res = subset(pval.res,pval.res$Pval != 0)
  #print(lodr.resLess)
  #print(lodr.resPlot)
  
  colnames(lodr.resLess)[1:3]<-c("|log2(Fold)|","MinError","Estimate")
  #print(paste(filenameRoot,kind,"LODR Estimates for P-value Threshold:",pval.cutoff,"; Probability:",prob))
  
  colnames(lodr.resPlot)[1:3]<-c("Ratio","MinError","Estimate")
  colnames(lodr.resLess)[1:3]<-c("Ratio","MinError","Estimate")
  lodr.resPlot <- as.data.frame(lodr.resPlot)
  lodr.resLess <- as.data.frame(lodr.resLess)
  lodr.resPlot$Ratio <- as.character(legendLabels)
  lodr.resLess$Ratio <- as.character(legendLabels) 
  annoTable = lodr.resLess[-c(2)]
  colnames(annoTable) <- c("Ratio",expression("LODR Estimate"), expression("90% CI Lower Bound"), expression("90% CI Upper Bound"))
  annoTable <- annoTable[-which(annoTable$Ratio == "1:1"),]
  cat("\n")
  print(annoTable)
  
  arrowDat = data.frame(Ratio = FCcode$Ratio, FC = FCcode$FC, x = lodr.resPlot[,3], y = pval.cutoff, xend = lodr.resPlot[,3], yend = 0)
  arrowDat$x[grep('<',lodr.resLess[,3])] <- Inf
  arrowDat = arrowDat[-which(arrowDat$FC == "1"),]
  arrowDat = arrowDat[which(is.finite(arrowDat$x)),]
  
  if(dim(arrowDat)[1] == 0){
    cat(paste("\nError! Estimated distribution of p-values does not cross \n",
              "threshold p-value, may be due to insufficient data quantity\n"))
    break
  }
  
  LODRplot = ggplot(pval.res) + geom_point(aes(x = MnCnt, y = Pval,colour = Ratio),size = 6) + scale_x_log10(limits = c(1,max(pval.res$MnCnt))) + scale_y_log10(breaks = c(1e-12,1e-10,1e-8,1e-5,1e-4,1e-3,1e-2,1e-1,1e0)) + geom_ribbon(data = lineDat, aes(x = x.new, y = fitLine, ymin=fitLower, ymax=fitUpper,fill = Ratio), alpha = 0.3,show_guide =F) + geom_line(data = lineDat,aes(x = x.new, y=fitLine, colour = Ratio),show_guide = F) + colScale + fillScale + xlab("Mean Counts") + ylab("DE Test P-values") + geom_hline(yintercept = pval.cutoff, linetype = 2, size = 2 ) + geom_segment(data = arrowDat, aes(x = x,y = y,xend = xend , yend = yend, colour = Ratio), lineend = "round",arrow = arrow(length =unit(0.5,"cm")), size = 2, alpha = 0.6) + theme(legend.justification=c(0,0), legend.position=c(0,0)) #+ annotation_custom(tableGrob(annoTable, show.rownames=F, parse = F), xmin = min(log(pval.res$MnCnt)), xmax = 0.25*log(mean(pval.res$MnCnt)), ymin =log(1e-12),ymax = 0  ) 
    
  ## create inset table
  my_table <- tableGrob(annoTable,show.rownames=F,gpar.coretext =gpar(fontsize=14),gpar.coltext=gpar(fontsize=14), gpar.rowtext=gpar(fontsize=14))
  
  
  ### final result 
  
#   Layout <- grid.layout(nrow = 2, ncol = 1, heights = unit(c(2, 0.5), c("null", "null")))
#   #grid.show.layout(Layout)
#   vplayout <- function(...) {
#     grid.newpage()
#     pushViewport(viewport(layout = Layout))
#   }
#   
#   subplot <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
#   
#   mmplot <- function(a, b) {
#     vplayout()
#     print(a, vp = subplot(1, 1))
#     print(b, vp = subplot(2, 1))
#   }
#   
#   annotLODRplot <- grob(mmplot(LODRplot, arrangeGrob(my_table)))
  #vplayout()
  annotLODRplot <- arrangeGrob(LODRplot, arrangeGrob(my_table), ncol = 1, heights = c(2,0.5))
  print(annotLODRplot)
  
  nam <- paste("plotLODR",kind, sep = ".")
  expDat$Figures$plotLODR <- annotLODRplot
  names(expDat$Figures)[which(names(expDat$Figures) == "plotLODR")] <- nam
  
  nam <- paste("lodr.res",kind,sep = ".")
  expDat$lodr.res <- lodr.resLess
  names(expDat)[which(names(expDat) == "lodr.res")] <- nam
  
  return(expDat)
}