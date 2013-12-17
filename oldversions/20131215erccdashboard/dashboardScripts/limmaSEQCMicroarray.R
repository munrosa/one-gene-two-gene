COHarray = read.delim("Data/QCmicroarrays_COH_UTSW/COHdata110909.txt")
COHarrayAB <- COHarray[-c(11:16)]
y <- COHarrayAB[c(2,5:10)]
row.names(y)<- y$TargetID
y <- y[-c(1)]
ylog <- log(y)
design <- cbind(Grp1=1,Grp1vs2=c(1,1,1,0,0,0))

fit <- lmFit(ylog,design)
fit <- eBayes(fit)
topTable(fit,coef=2)

