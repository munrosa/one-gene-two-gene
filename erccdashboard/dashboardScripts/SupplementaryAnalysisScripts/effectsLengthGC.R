
sampleInfo = list(study = "SEQC_Main", #
                  siteName = "MAY",
                  platform = "ILM",
                  analysis = "NCTR",
                  printPDF = F,
                  DEtest = T,
                  totalSeqReads = T,
                  libeSizeNorm = T,
                  sample1Name = "E",
                  sample2Name = "F",
                  FCcode = data.frame(Ratio = c("a","b","c","d"),
                                      FC =  c(4,1,.667,.5)),
                  myYLimMA = c(-3.5,3.5),
                  myXLim = c(-15,25),
                  myYLim = NULL,
                  xlimEffects = c(-15,15),
                  choseFDR = 0.01,
                  legendLabels = c("4:1","1:1","1:1.5","1:2"),
                  ERCCdilution = 1,
                  spikeVol = 1,
                  totalRNAmass = 1
)

source('~/Documents/NISTMunro/Projects/testerccdashboard/dashboardScripts/130924_PurePoolsEandF_MainSEQC.R')

MAYeffects <- dynRangeDat$effects

sampleInfo$siteName <- "BGI"
source('~/Documents/NISTMunro/Projects/testerccdashboard/dashboardScripts/130924_PurePoolsEandF_MainSEQC.R')

BGIeffects <- dynRangeDat$effects

sampleInfo$siteName <- "CNL"
source('~/Documents/NISTMunro/Projects/testerccdashboard/dashboardScripts/130924_PurePoolsEandF_MainSEQC.R')

CNLeffects <- dynRangeDat$effects

### TO DO Sarah Munro make Lines 1 - 38 into a loop so you don't have to repeat code
BGIeffects$site <- "Lab 2"
CNLeffects$site <- "Lab 3"
MAYeffects$site <- "Lab 5"
effects <- rbind(BGIeffects,CNLeffects)
effects <- rbind(effects,MAYeffects)
str(effects)
effects$site <- as.factor(effects$site)
effects_l <- melt(effects)




library(ggplot2)
theme_update(legend.justification=c(0,0), legend.position=c(0,0))
#effects <- dynRangeDat$effects
erccDef <- read.csv("data/ERCCDef.csv")

#minLabel = effects$ERCC.effect <= fivenum(effects$ERCC.effect)[2]- 1.5*IQR(effects$ERCC.effect)
#maxLabel = effects$ERCC.effect >= fivenum(effects$ERCC.effect)[4]+ 1.5*IQR(effects$ERCC.effect)

#Create a custom color scale
myColors <- c("#FF9900","#339966", "#6699CC", "#CC6666")
names(myColors) <- levels(effects$Ratio)
colScale <- scale_colour_manual(name = "Ratio",values = myColors, labels = sampleInfo$legendLabels)
fillScale <- scale_fill_manual(name = "Ratio", values = myColors, labels = sampleInfo$legendLabels)

EFeffects <- merge(erccDef,effects)

# Plot ERCC effects as a function of length
ggplot(EFeffects)+ geom_point(aes(x = Length, y = ERCC.effect, colour = Ratio), size = 6)+facet_grid(~site) + colScale + ylim(-30,30) + xlab("Log2 Average ERCC Amount in Pure Pool (attomol nt/uL)")+ ylab("Per ERCC Differences from Linear Fit Standardized by s.e., unitless)")+ geom_text(data=subset(EFeffects,(ERCC.effect <= fivenum(ERCC.effect)[2] - (1.5)*IQR(ERCC.effect))),aes(x = Length, y = ERCC.effect, label = gsub("ERCC-00","",Feature)),colour = "black",show_guide = F,angle = 45,hjust = -0.25, position = position_jitter(width=0.5)) + geom_abline(aes(intercept = 0, slope = 0), colour = "grey")