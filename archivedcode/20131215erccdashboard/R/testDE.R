testDE<- function(sampleInfo, expDat, cnt = cnt, info = info){
  library(QuasiSeq)
  #library(DESeq)
  library(edgeR)
  
  filenameRoot = sampleInfo$filenameRoot
  legendLabels = sampleInfo$legendLabels
  FCcode = sampleInfo$FCcode
  totalSeqReads = sampleInfo$totalSeqReads
  
  
  idCols = expDat$idColsAdj
  r_m.mn = expDat$r_m.res$r_m.mn 
  totalReads = expDat$totalReads 
  sample1 = expDat$sampleNames[1]
  sample2 = expDat$sampleNames[2]
  
  #Create a custom color scale
  myColors <- c("#339966","#FF9900", "#6699CC", "#CC6666")
  names(myColors) <- levels(FCcode$Ratio)
  colScale <- scale_colour_manual(name = "Ratio",values = myColors, labels = legendLabels)
  fillScale <- scale_fill_manual(name = "Ratio", values = myColors, labels = legendLabels)
  #Move legend
  theme_update(legend.justification=c(1,1), legend.position=c(1,1))
  
  cnt = unique(cnt)
  Features = make.names(cnt$Feature,unique=T)
  Features = gsub(".","-", Features, fixed = T)
  #print(str(cnt))
  #print(dim(cnt))
  rownames(cnt)<-Features; cnt<-as.matrix(cnt[,-1])
  
  colnames(cnt)<-paste(rep(c(sample1,sample2),each=ncol(cnt)/2),c(1:(ncol(cnt)/2),1:(ncol(cnt)/2)),sep="")
  ## Get ERCC names
  ERCC<-rownames(cnt[substr(rownames(cnt),1,5)=="ERCC-",])
  
  ## Specify Sample (A or B)
  trt<-rep(1:2,each=ncol(cnt)/2)
  
  design.list<-list(trt,rep(1,ncol(cnt)))
  
  ## Compute offset (e.g. total counts, 75% quantile, TMM, etc)
  if(totalSeqReads == F){
    log.offset<-log(colSums(cnt))
    cat("\nUsing Mapped Reads\n")
    cat(colSums(cnt),"\n")  
  }else{
    log.offset<- log(totalReads)
    cat("\nUsing Total Reads\n")
    cat(totalReads,"\n")
  }
  
  ERCC.FC = idCols[c(1,4)];rownames(ERCC.FC)<-ERCC.FC[,1]
  
  ERCC.FC$NumRatio <- 1
  for (i in 1:nlevels(FCcode$Ratio)){
    ERCC.FC$NumRatio[which(ERCC.FC$Ratio == as.character(FCcode$Ratio[i]))] = FCcode$FC[i]  
  }
  ERCC.Ratio = ERCC.FC[c(1,2)]
  ERCC.FC = ERCC.FC[-c(2)]
    
  group <- as.factor(trt)
  d <- DGEList(counts=cnt,group=group)
  # use log.offset for the library size
  d$samples$lib.size <- exp(log.offset) 
  #print(d$samples$lib.size)
  #Dispersion trend
  design <- model.matrix(~group)
  d1 <- estimateGLMCommonDisp(d,design,verbose=TRUE)
  d1 <- estimateGLMTrendedDisp(d1,design)
  d1 <- estimateGLMTagwiseDisp(d1,design)

  #######################################################
  ###  Simulate ERCC data from negative binomial fit ####
  ###  (to be used for sim-based LODR)               ####
  #######################################################
  # Define function simcnt.lodr
  simcnt.lodr<-function(cnt,disp,trt,fold,log.offset){
    #### 'cnt' is the matrix of counts for endogenous genes
    #### 'disp' is the central trend fitted to the estimated dispersions (vs. expression)
    #### 'trt' is a vector specifying to which treatment each column in 'cnt' belongs
    #### 'fold' is the desired fold change (trt 2/trt 1)
    #### 'log.offset' is used to account for differences in library size
    
    sim.ind<-round(nrow(cnt)*c(1:49)/50)  ### mimick every 797th gene (roughly), when sorted by total count
    norm.cnt<-t(t(cnt)/exp(log.offset))
    sim.mn<-matrix(sort(rowMeans(norm.cnt))[sim.ind],nrow=length(sim.ind),ncol=ncol(cnt))
    sim.mn<-sim.mn*2/(1+fold)
    sim.mn[,trt==2]<-sim.mn[,trt==2]*fold
    sim.mn<-t(t(sim.mn)*exp(log.offset))
    
    sim.disp<-disp[order(rowMeans(norm.cnt))][sim.ind]
    size<-matrix(1/sim.disp,length(sim.disp),ncol(cnt))
    simcnt<-matrix(rnbinom(length(sim.disp)*ncol(cnt),mu=sim.mn,size=size),length(sim.disp),ncol(cnt))
    rownames(simcnt)<-paste("Sim",fold,"Fold",1:length(sim.ind),sep="")
    return(simcnt)
  }
  
  #### Simulate data for each of the fold changes used in the ERCCs
  simcnt<-NULL
  for(fold in unique(ERCC.FC[!is.na(ERCC.FC[,2]),2])){
    #simcnt<-rbind(simcnt,simcnt.lodr(cnt,disp=fitInfo(cds)$fittedDispEsts,trt=trt,fold=fold,log.offset=log.offset))
    simcnt<-rbind(simcnt,simcnt.lodr(cnt,disp=d1$trended.dispersion,trt=trt,fold=fold,log.offset=log.offset))
  }
  
  ##### Analyze combination of observed and simulated data with edgeR
  group <- as.factor(trt)
  d2 <- DGEList(counts=rbind(cnt,simcnt),group=group)
  d2$samples$lib.size <- exp(log.offset)
  
  design <- model.matrix(~group)
  d2 <- estimateGLMCommonDisp(d2,design,verbose=TRUE)
  d2 <- estimateGLMTrendedDisp(d2,design)
  d2 <- estimateGLMTagwiseDisp(d2,design)
 
  NBdisp<-d2$tagwise.dispersion;names(NBdisp)<-rownames(rbind(cnt,simcnt))
  NBdisptrend<-d2$trended.dispersion;names(NBdisptrend)<-rownames(rbind(cnt,simcnt))
  
  ### Plot estimated dispersions 
  dispcnt = data.frame(x = rowMeans(rbind(cnt,simcnt)), y = NBdisp)
  dispSimcnt = data.frame(xsim = rowMeans(simcnt), ysim = NBdisp[-c(1:nrow(cnt))])
  dispERCCcnt = data.frame(xERCC = rowMeans(cnt)[ERCC],yERCC = NBdisp[ERCC])
  
  dispPlot = ggplot(dispcnt) + geom_point(aes(x = x, y = y)) + scale_x_log10()+ xlab("Average Count") + ylab("Estimated Dispersion") + geom_point(data = dispSimcnt,aes(x = xsim, y = ysim),colour = 3) + geom_point(data = dispERCCcnt, aes(x = xERCC, y = yERCC),colour = 2)
   #print(dispPlot) 
  
  ### Run QuasiSeq
  use.fit<-QL.fit(rbind(cnt,simcnt),design.list,log.offset=log.offset,NBdisp=NBdisptrend)# using trended dispersion from edgeR
  use.res<-QL.results(use.fit, Plot = F)
  
  useRescnt <- data.frame(use.res$P.values) ### added to allow next command to access the Spline P-values
  
  ### Collect results for simulated data; to be passed along to LODR function
  sim.pval.res<-cbind(row.names(simcnt), rowMeans(simcnt),use.res$P.values$QLSpline[-(1:nrow(cnt))],use.res$log.P.values$QLSpline[-(1:nrow(cnt))],use.res$F.stat$QLSpline[-(1:nrow(cnt))],
  rep(unique(ERCC.FC[!is.na(ERCC.FC[,2]),2]),each=49),rep(use.res$d0[2],nrow(simcnt)))
  colnames(sim.pval.res)<-c("Feature","MnCnt","Pval","LogPval","F.stat","Fold",
                            "Den.df")
  rownames(sim.pval.res)<-rownames(simcnt)

  write.csv(sim.pval.res,file=paste(filenameRoot,"Sim Pvals.csv"))
  
  ## remove results for simulated data
  use.fit2<-use.fit
  use.fit2$LRT<-matrix(use.fit2$LRT[1:nrow(cnt),],nrow(cnt),1)
  use.fit2$phi.hat.dev<-use.fit2$phi.hat.dev[1:nrow(cnt)]
  use.fit2$phi.hat.pearson<-use.fit2$phi.hat.pearson[1:nrow(cnt)]
  use.fit2$mn.cnt<-use.fit2$mn.cnt[1:nrow(cnt)]
  use.fit2$NB.disp<-use.fit2$NB.disp[1:nrow(cnt)]
  use.fit2$fitted.values<-use.fit2$fitted.values[1:nrow(cnt),]
  use.fit2$coefficients<-use.fit2$coefficients[1:nrow(cnt),]
  
  use.res2<-QL.results(use.fit2, Plot = F)
  
  ###################################
  #### Examine results for ERCCs ####
  ###################################
  pvals<-use.res2$P.values$QLSpline
  names(pvals)<-rownames(cnt)
  ERCC.pvals<-pvals[ERCC]
  
  rownames(use.fit2$coefficients)<-rownames(cnt)
  est.FC<-use.fit2$coefficients[ERCC,2]
  
  ## Reanalyze ERCC transcripts using adjusted offsets to center fold change estimates
  adj <- r_m.mn #### Use r_m estimated from NegBin GLM
  use.fit.adj<-QL.fit(cnt[ERCC,],design.list,log.offset=log.offset-rep(c(adj,0),each=ncol(cnt)/2),NBdisp=NBdisptrend[ERCC])

  rownames(use.fit.adj$coefficients)<-ERCC;
  est.FC.adj<-use.fit.adj$coefficients[ERCC,2]
  use.fit3<-use.fit2
  rownames(use.fit3$LRT)<-rownames(cnt);use.fit3$LRT[ERCC,]<-use.fit.adj$LRT
  use.res.adj<-QL.results(use.fit3, Plot = F)

  pvals<-use.res.adj$P.values$QLSpline

  log.pvals<-use.res.adj$log.P.values$QLSpline
  F.stat<-use.res.adj$F.stat$QLSpline
  names(F.stat)<-names(log.pvals)<-names(pvals)<-rownames(cnt)

  qvals<-use.res.adj$Q.values$QLSpline
  quasiSeq.res = data.frame(Feature = names(pvals), pvals = pvals, qvals = qvals, log.pvals=log.pvals, F.stat=F.stat, den.df=rep(use.res.adj$d0[2], length(pvals)))  

  write.csv(quasiSeq.res, file = paste(filenameRoot,"quasiSeq.res.csv",sep="."),row.names = F)
  
  ERCC.pvals.adj<-pvals[ERCC]
  ERCC.F.stat.adj<-F.stat[ERCC]
  ERCC.log.pvals.adj<-log.pvals[ERCC]

  ### Collect results for ERCCs; to be passed along to LODR function
  pval.res<-data.frame(row.names(cnt[ERCC,]),rowMeans(cnt[ERCC,]),ERCC.pvals.adj,ERCC.log.pvals.adj,ERCC.F.stat.adj,ERCC.FC[ERCC,2],rep(use.res.adj$d0[2],length(ERCC.pvals.adj)))
  colnames(pval.res)<-c("Feature","MnCnt","Pval","LogPval","F.stat","Fold",
                        "Den.df")
  print(str(pval.res))
  row.names(pval.res) <- NULL
  write.csv(pval.res,file=paste(filenameRoot,"ERCC Pvals.csv"))
  print("Finished DE testing")
  
  expDat$quasiSeq.res <- quasiSeq.res
  expDat$ERCCpvals <- pval.res
  #################################################
  ###Examine ERCC blending with endogenous genes###
  #################################################
  
  ### Plot estimated dispersions and central trend estimates ###
  
  NBdisp1 <- d1$tagwise.dispersion;names(NBdisp1)<-rownames(cnt) # use tagwise
  NBdisp2<-d1$trended.dispersion;names(NBdisp2)<-rownames(cnt) # use trended?

  dispcnt = data.frame(Feature = rownames(cnt), x = rowMeans(cnt), y = NBdisp1)
  dispcntSort = data.frame(xsort = sort(rowMeans(cnt)), ysort = NBdisp2[order(rowMeans(cnt))])
    
  dispERCCcnt = data.frame(meanCnt = rowMeans(cnt)[ERCC],disp = NBdisp1[ERCC], Ratio = ERCC.Ratio$Ratio[match(ERCC,ERCC.Ratio$Feature)] )
  dispERCCcnt$Ratio <- as.factor(dispERCCcnt$Ratio) 
  
  dispPlot = ggplot(dispcnt) +  geom_point(data = dispcnt, aes(x = x, y = y),colour = "grey85",size = 5, alpha = 0.6) + geom_point(data = dispERCCcnt, aes(x = meanCnt, y = disp, colour = Ratio), size = 5, alpha = 0.6) + xlab("Mean Counts") + ylab("Negative Binomial Dispersion") + stat_smooth(data = dispcntSort,aes(x = xsort, y = ysort),colour = "black", alpha = 0.8)  + colScale + scale_x_log10()
  #print(dispPlot)
  
  
  phi.hat<-use.fit$phi.hat.dev[1:nrow(cnt)]
  mn.cnt<-rowMeans(cnt);den.df=length(trt)-length(unique(trt))
  y<-log(phi.hat);names(y)<-rownames(cnt); y[y==-Inf]<-min(y[y!=-Inf]); y[y==Inf]<-max(y[y!=Inf])
  spline.fit<-sreg(log(mn.cnt), y)
  fit.method<-"spline"
  
  if(spline.fit$eff.df==0){
    fit.method<-"locfit curve"
    print("Spline fitting failed.  Using locfit instead.")
    fit<-locfit(y~lp(log(mn.cnt)))
    spline.fit<-list(fitted.values=predict(fit,newdata=log(mn.cnt)),eff.df=fit$dp["df1"])
  }
  
  ### We use Smyth's (2004) approach from LIMMA to estimate parameters for prior distributions 
  ### of gene specific dispersion estimates. The function below also provides resulting 
  ### shrunken point estimates of dispersion 
  
  #### Code for implementing Smyth's approach begins here ####
  shrink.phi<-function(phi.hat,den.df){
    phi.hat[phi.hat<=0]<-min(phi.hat[phi.hat>0])
    z<-log(phi.hat); z[z==Inf]<-max(z[z!=Inf]); z[z==-Inf]<-min(z[z!=-Inf]);mnz<-mean(z)
    
    ## solve for d0 and phi0
    d0arg<-var(z)-trigamma((den.df)/2)
    if(d0arg>0){
      dif<-function(x,y) abs(trigamma(x)-y)
      inverse.trigamma<-function(y) optimize(dif,interval=c(0,10000),y=y)$minimum
      d0<-2*inverse.trigamma(d0arg)
      phi0<-exp(mnz-digamma((den.df)/2)+digamma(d0/2)- log(d0/(den.df)))
      
      ## compute shrunken phi's
      phi.shrink<-((den.df)*phi.hat+d0*phi0)/(den.df+d0)  }
    else{phi.shrink<-rep(exp(mnz),length(z)); d0<-Inf; phi0<-exp(mnz) }
    return(list(phi.shrink=phi.shrink,d0=d0,phi0=phi0))  }
  #### Code for implementing Smyth's approach ends here ####
  
  ### Obtain estimate for prior degrees of freedom and scaling factor after adjusting for cubic spline fit
  y2<-phi.hat/exp(spline.fit$fitted.values)
  shrink<-shrink.phi(y2,den.df)
  D0<-shrink[[2]]; 
  phi0<-shrink[[3]]; print(paste("Spline scaling factor:",phi0))
  
  ### plot the resulting cubic spline fit
  #mean = log(mn.cnt)
  mean = mn.cnt
  dispcnt = data.frame(mean = mean,y = y)
  dispcntSort = data.frame(xsort = sort(mean), ysort = spline.fit$fitted.values[order(mn.cnt)])
  dispERCC = data.frame(mean = mean[ERCC],y = y[ERCC], Ratio= ERCC.Ratio$Ratio[match(ERCC,ERCC.Ratio$Feature)] )
  dispERCC$Ratio <- as.factor(dispERCC$Ratio) 
  
  ###Compare quantiles from estimated theoretical and empirical dispersion distributions
  RR<-NULL; sort.mn<-log(sort(mn.cnt)); qq<-c(.05,.95);ord.F<-y[order(mn.cnt)]
  bins<-c(1+round(length(y)*(0:19)/20),length(y))
  for(ii in 1:(length(bins)-1)){
    ind<-bins[ii]:bins[ii+1]
    RR<-rbind(RR,c(quantile(ord.F[ind],qq),median(sort.mn[ind])))
  }
  dispcntInterval = data.frame(x=sort(log(mn.cnt)),upperY = spline.fit$fitted.values[order(mn.cnt)]+log(phi0*qf(.95,den.df,D0)), lowerY =spline.fit$fitted.values[order(mn.cnt)]+log(phi0*qf(.05,den.df,D0)) )
  
  dispcntRR = data.frame(x = RR[,3], upperY = RR[,2], lowerY = RR[,1])
  
  quasiDispPlot = ggplot() + geom_point(data = dispcnt, aes(x = mean, y = y),colour = "grey85", size = 5, alpha = 0.6) + geom_point(data = dispERCC, aes(x = mean, y = y,colour = Ratio), size = 5, alpha = 0.6) + xlab("Mean Counts") + ylab("log(Quasi Dispersion)") + stat_smooth(data = dispcntSort,aes(x = xsort, y = ysort),colour = "black") + scale_x_log10() + colScale
  expDat$Figures$dispPlot <- quasiDispPlot
  #save(quasiDispPlot, file=paste(filenameRoot,"DispPlot.RData", sep = "."))
  cat("\nFinished examining dispersions\n")
  return(expDat)

}