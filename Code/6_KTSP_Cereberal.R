############################################################################


rm(list = ls())

setwd("~/Documents/Research/Projects/Malaria")

renv::activate("~/Documents/Research/Projects/Malaria")

### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(enrichR)
library(mltools)
library(xtable)
library(dplyr)
library(superheat)

## Load data
load("./Objs/MalariaDataGood_NCvsC.rda")

# Quantile normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

boxplot(mixTrainMat)

boxplot(usedTrainMat)

# Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

# Load mechanistic pairs
Up <- c("PUM2", "SETX", "RABEP1", "ELF2", "ZNF197", "KRIT1", "EPHA4", "XRCC5", "LARP4", "SCN2B", "CDH8", "MREG", "TTC17", "THRB")
Down <- c("TRIP12", "MYH11", "ANK2", "CD53", "SPATS2L", "KPNA6", "CHRNA10", "ASB7", "C18orf8")

myTSPs <- expand.grid(Up, Down)
myTSPs <- as.matrix(myTSPs)
colnames(myTSPs) <- c("Up", "Down")

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))

#UsedTrainMat <- UsedTrainMat[keepGns, ]
#UsedTestMat <- UsedTestMat[keepGns, ]

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #8
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=9, featureNo= featNo, 
  FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = myTSPs)

ktspPredictorRes

#save(ktspPredictorRes, file = "./Objs/KTSP_Model_Cerebral.rda")
load("./Objs/KTSP_Model_Cerebral.rda")

###########################################################################
############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("nonCerebral", "cerebral"), direction = "<"), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("nonCerebral", "cerebral"), direction = "<"), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = T, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("nonCerebral", "cerebral"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")
ROCTrain

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) >= thr)

### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "cerebral", mode = "everything")
ConfusionTrain

MCC_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, ConfusionTrain$overall["Accuracy"], ConfusionTrain$byClass["Balanced Accuracy"], ConfusionTrain$byClass["Sensitivity"], ConfusionTrain$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
ROCTest <- roc(usedTestGroup, ktspStatsTestRes$statistics, plot = T, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("nonCerebral", "cerebral"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROCTest

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) >= thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "cerebral", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, ConfusionTest$overall["Accuracy"], ConfusionTest$byClass["Balanced Accuracy"], ConfusionTest$byClass["Sensitivity"], ConfusionTest$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
CereberalKTSP_Perf <- cbind(TrainPerf, TestPerf)

# Save
save(CereberalKTSP_Perf, file = "./Objs/CereberalKTSP_Perf.rda")

########################################################################
########################################################################
## Plot AUC of both the agnostic and mechanistic classifiers using ggplot2
### Prepare the legend
forLegend_KTSP <- apply(rbind(
  ci(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("nonCerebral", "cerebral"), direction = "<")),
  ci(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("nonCerebral", "cerebral"), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


#################################################################
### ROC curves Using ggplot2

### Training
datTrn_KTSP <- melt(data.frame(
  ## Training Group
  Training=factor(usedTrainGroup, levels = c("nonCerebral", "cerebral")),
  ## Mechanistic KTSP SUM training
  KTSP.Training=ktspStatsTrainRes$statistics))
### Change Colnames
colnames(datTrn_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")


### Testing
datTst_KTSP <- melt(data.frame(
  ## Testing group
  Testing=factor(usedTestGroup, levels = c("nonCerebral", "cerebral")),
  ## Mechanistic KTSP SUM training
  KTSP.Testing=ktspStatsTestRes$statistics))
### Change Colnames
colnames(datTst_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")

### Combine
dat_KTSP <- rbind(datTrn_KTSP, datTst_KTSP)
dat_KTSP$Status <- as.numeric(dat_KTSP$Status)-1

### Replace levels
levels(dat_KTSP$KTSP_type) <- gsub("\\.", "-", levels(dat_KTSP$KTSP_type))
levels(dat_KTSP$KTSP_type) <- paste(levels(dat_KTSP$KTSP_type), forLegend_KTSP[c(1,2)])

#################################################################
### Plot Curve
tiff("./Figs/AUCggplot_Cerebral.tiff",
    width=2800, height=2500, res=300)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "Performance of the 9-TSPs cerebral malaria signature in the training and testing data"
#legendTitle <- paste("Mechanistic (", nrow(ktspPredictorRes$TSPs), " pairs)",
#                     "Agnostic (", nrow(ktspPredictorUnRes$TSPs), " pairs)",  sep="")
### Plot
basicplot_KTSP <- ggplot(dat_KTSP, aes(d=Status, m=KTSP_sum, color=KTSP_type,
                                       linetype = KTSP_type)) +
  geom_roc(cutoffs.at = seq(1,20,1)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=14, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 14),
        legend.title = element_text(face="bold", size=0)) +
  #scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  #scale_linetype_manual(legendTitle, values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_KTSP
### Close device
dev.off()

save(basicplot_KTSP, file = "./Objs/BasicPlot_KTSP_Cerebral.rda")

##################################################################################################
### Heatmaps
##################################################################################################
#### Make a heatmap of the TSPs expression in the 3 datasets


Cerebral_Stats <- ktspStatsTestRes$comparisons
Cerebral_Stats <- Cerebral_Stats*1

# Order the stats by the sum of votes
NewOrder <- order(rowSums(Cerebral_Stats))
Cerebral_Stats <- Cerebral_Stats[NewOrder, ]

#GroupCol <- ArrayGroup[NewOrder]
#levels(GroupCol) <- c("red", "blue")

# Order the true class labels
CerebralGroup <- usedTestGroup
CerebralGroup <- CerebralGroup[NewOrder]

# Get the predicted class labels
PredClass <- usedTestPredictionRes
#PredClass <- factor(PredClass, levels = c("POS", "NEG")) 
levels(PredClass) <- c("blue", "red")
PredClass <- PredClass[NewOrder]
# Check order 
all(rownames(Cerebral_Stats) == names(PredClass)) # TRUE:: No Need to order them (Already ordered by sum of votes) 

# Plot
png(filename = "./Figs/CerebralKTSP_Heatmap.png", width = 3000, height = 2000, res = 200)
superheat(Cerebral_Stats, col.dendrogram = F, yr = rowSums(Cerebral_Stats), yr.plot.type  = "bar", left.label.text.col = c("blue", "red"), membership.rows = CerebralGroup, bottom.label.text.size = 3.5, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the testing data")
dev.off()

##########################
sessionInfo()
