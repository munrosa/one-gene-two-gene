loadERCCInfo <- function(expDat, erccMixes){
  # Get the ERCC Mix definition file provided by user and combine it with the package ERCCDef file
  load(file = system.file("data/ERCC.RData", 
                          package = "erccdashboard"))
  ERCCMix1and2SEQC <- ERCCMix1and2
  ERCCDef <- ERCCDef
  MixDef <- ERCCMix1and2SEQC
  names(MixDef)[1:2] <- c("Feature","Ratio") 
  if(erccMixes == "Ambion4plexPair"){
    
    
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
    expDat$sampleInfo <- append(expDat$sampleInfo,list(idColsSRM = idCols, 
                                                       MixDef = MixDef, 
                                                       FCcode = FCcode, 
                                                       legendLabels = 
                                                         legendLabels) )
    return(expDat)
    
  }
  if(erccMixes == "AmbionSingle"){
    
    # Sort by the feature column
    MixDefSort <- MixDef[do.call(order, MixDef[c(1)]), ]
    MixDef <- MixDefSort
    
    # Check that the MixDef has the same ERCCs as the ERCCDef data frame
    ERCCmatch <- ERCCDef[match(MixDef$Feature,ERCCDef$Feature),]
    
    # Combine the Mix definition and ERCCdef data frames
    idCols <- merge(ERCCmatch,MixDef, by = "Feature")
    
    # Fix the names and column identity, assigning the mix 1 concentration 
    # values to both columns
    idCols[c(6)] <- idCols[c(5)] 
    names(idCols)[5:6] <- c("Conc1", "Conc2")

    FCcode = data.frame(Ratio = c("a","b","c","d"), # Default for Ambion pool
                        FC =  c(1,1,1,1))
    legendLabels = c("1:1a","1:1b","1:1c","1:2d") # Default for Ambion pool
    expDat$sampleInfo <- append(expDat$sampleInfo,list(idColsSRM = idCols, 
                                                       MixDef = MixDef, 
                                                       FCcode = FCcode, 
                                                       legendLabels = 
                                                         legendLabels) )
    return(expDat)
    
  }
  else{
  print("Invalid ERCC mixture selected")
  break
  }
  
}