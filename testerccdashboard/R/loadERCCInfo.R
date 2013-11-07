loadERCCInfo <- function(expDat, erccMixes = "Ambion4plexPair"){
  
  if(erccMixes == "Ambion4plexPair"){
    # Get the ERCC Mix definition file provided by user and combine it with the package ERCCDef file
    ERCCMix1and2SEQC <- read.csv(system.file("data/ERCCMix1and2SEQC.csv",
                                             package = "testerccdashboard"))
    ERCCDef <- read.csv(system.file("data/ERCCDef.csv",
                                    package = "testerccdashboard"))
    
    MixDef <- ERCCMix1and2SEQC
    names(MixDef)[1:2] <- c("Feature","Ratio") 
    
    # Sort by the feature column
    MixDefSort <- MixDef[do.call(order, MixDef[c(1)]), ]
    MixDef <- MixDefSort
    
    # Check that the MixDef has the same ERCCs as the ERCCDef data frame
    ERCCmatch <- ERCCDef[match(MixDef$Feature,ERCCDef$Feature),]
    
    # Combine the Mix definition and ERCCdef data frames
    idCols <- merge(ERCCmatch,MixDef, by = "Feature")
    
    # Fix the names
    names(idCols)[5:6] <- c("Conc1", "Conc2")
    
    FCcode = data.frame(Ratio = c("a","b","c","d"), # Default for Ambion pool
                        FC =  c(4,1,.667,.5))
    legendLabels = c("4:1","1:1","1:1.5","1:2") # Default for Ambion pool
    expDat$sampleInfo <- append(expDat$sampleInfo,list(idColsSRM = idCols, MixDef = MixDef, FCcode = FCcode, legendLabels = legendLabels) )
    return(expDat)
    
  }
  else{
  print("Invalid ERCC mixture selected")
  break
  }
  
}