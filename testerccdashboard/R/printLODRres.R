printLODRres <- function(expDat, dynRangeDat, lodr.res){
  sampleInfo <- expDat$sampleInfo
  fit.coeff <- dynRangeDat$fit.coeff
  mnLibeFactor <- expDat$mnLibeFactor
  
  FCcode = sampleInfo$FCcode
  legendLabels = sampleInfo$legendLabels
  

  lodr.res = data.frame(lodr.res)
  
  Fold = lodr.res[c(1)]
  Count = as.numeric(gsub("<", "",lodr.res$Estimate))
  Ratio = legendLabels 
  
  
  # Convert LODR count estimate to library size normalized
  logCount = log2((Count/(mnLibeFactor)))#+1)
  
  ConcEst = (log2(Count/(mnLibeFactor))-fit.coeff[1])/fit.coeff[2]
  
  LODR.print.res = data.frame(Fold, Ratio, Count, logCount, ConcEst)
  
  names(LODR.print.res)<- c("Fold","Ratio","Count","Log2Count_normalized","Log2Conc")
  print(LODR.print.res)
  
  cutoffs = ConcEst[which(!(is.na(ConcEst)))]
  
  return(list(LODRtable = LODR.print.res, cutoffs=cutoffs))
  
}