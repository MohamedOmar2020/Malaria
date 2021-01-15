############################################################################


rm(list = ls())

#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/BreastChemo")

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

# Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

# Load mechanistic pairs
Up <- c("PRRX1",
        "ITGA9",
        "ANXA6",
        "CALU",
        "MAPK4",
        "CA12",
        "MMP11",
        "COL4A5",
        "EGLN3",
        "MSX1")

Down <- c("FGL2",
          "TLR8",
          "IFIT3",
          "GBP2",
          "ALAS2",
          "PCDH9",
          "CD53",
          "SLAMF7",
          "CD36",
          "CD48")

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
  usedTrainMat, usedTrainGroup, krange=10, featureNo= featNo, 
  FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = myTSPs)

ktspPredictorRes

#save(ktspPredictorRes, file = "./Objs/KTSP_Model_Cerebral.rda")
#load("./Objs/KTSP_Model_Cerebral.rda")

#Mechanistic_KTSP <- cbind(ktspPredictorRes$TSPs, ktspPredictorRes$score)
#colnames(Mechanistic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Mechanistic_KTSP, type = "latex"), file = "./Objs/KTSP/Mechanistic.tex")

## Save the mechanistic Pairs
#MechanisticKTSP_Pairs <- c("ZKSCAN1>ZW10", "HES1>NR3C1", "ZNF160>ZMIZ1", "KLHL9>PML", "NUDC>RELA", "SNCA>ZBTB7A", "ALDH7A1>NFIC", "SMC3>RNF44", "UBE2N>SUZ12", "ARL2>MAX", "EGR1>TSTA3")
#save(MechanisticKTSP_Pairs, file = "./Objs/KTSP/MechanisticKTSP_Pairs.rda")


###########################################################################
### Check consistency with biology
# keep <- ktspPredictorRes$TSPs[,1] %in% myTSPs[,"Up"] & ktspPredictorRes$TSPs[,2] %in% myTSPs[,"Down"]
# table(keep)
# 
# ###Subset
# ktspPredictorRes$name <- paste(sum(keep), "TSPs", sep="")
# ktspPredictorRes$TSPs <- ktspPredictorRes$TSPs[ keep, ]
# ktspPredictorRes$score <- ktspPredictorRes$score[ keep ]
# ktspPredictorRes$tieVote <- droplevels(ktspPredictorRes$tieVote[keep])
# 
# #### Visualize the classifier
# ktspPredictorRes
# 
# save(ktspPredictorRes, file = "./Objs/MechanisticKTSP_Filt.rda")

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
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, print.thres=thr, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("nonCerebral", "cerebral"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")
ROCTrain

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) >= 1.5)

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

#KTSP_STATs_Test_Mechanistic <- t(ktspStatsTestRes$comparisons)
#KTSP_STATs_Test_Mechanistic[KTSP_STATs_Test_Mechanistic == FALSE] <- 0

#save(KTSP_STATs_Train_Mechanistic, KTSP_STATs_Test_Mechanistic, file = "./Objs/KTSP/TNBC_KTSP_STATs_Mechanistic.rda")

## Threshold
#thr_test <- coords(roc(UsedTestGroup, ktspStatsTestRes$statistics, levels = c("Resistant", "Sensitive")), transpose = T,"best")["threshold"]
#thr_test

#####
## Print ROC curve local maximas
#coords(roc(UsedTestGroup, ktspStatsTestRes$statistics, levels = c("Resistant", "Sensitive")), transpose = T, "local maximas")

## Plot curve
ROCTest <- roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.thres=thr_test, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("nonCerebral", "cerebral"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROCTest

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) >= thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "cerebral", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Test
