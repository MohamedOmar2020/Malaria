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
load("./Objs/MalariaDataGood_Comp.rda")

# Quantile normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

# Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #11
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=10, featureNo= featNo, 
  FilterFunc = SWAP.Filter.Wilcoxon)

ktspPredictorRes

save(ktspPredictorRes, file = "./Objs/KTSP_Model.rda")
load("./Objs/KTSP_Model.rda")
#Mechanistic_KTSP <- cbind(ktspPredictorRes$TSPs, ktspPredictorRes$score)
#colnames(Mechanistic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Mechanistic_KTSP, type = "latex"), file = "./Objs/KTSP/Mechanistic.tex")

## Save the mechanistic Pairs
#MechanisticKTSP_Pairs <- c("ZKSCAN1>ZW10", "HES1>NR3C1", "ZNF160>ZMIZ1", "KLHL9>PML", "NUDC>RELA", "SNCA>ZBTB7A", "ALDH7A1>NFIC", "SMC3>RNF44", "UBE2N>SUZ12", "ARL2>MAX", "EGR1>TSTA3")
#save(MechanisticKTSP_Pairs, file = "./Objs/KTSP/MechanisticKTSP_Pairs.rda")


###########################################################################
### Check consistency with biology
# keep <- ktspPredictorRes$TSPs[,1] %in% myTSPs[,"GoodGene"] & ktspPredictorRes$TSPs[,2] %in% myTSPs[,"BadGene"]
# table(keep)
# #
# # ###Subset
# ktspPredictorRes$name <- paste(sum(keep), "TSPs", sep="")
# ktspPredictorRes$TSPs <- ktspPredictorRes$TSPs[ keep, ]
# ktspPredictorRes$score <- ktspPredictorRes$score[ keep ]
# ktspPredictorRes$tieVote <- droplevels(ktspPredictorRes$tieVote[keep])
# 
# # ### Visualize the classifier
# ktspPredictorRes
# # 
# save(ktspPredictorRes, file = "./Objs/MechanisticKTSP_Filt.rda")

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("unComplicated", "Complicated"), direction = "<"), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("unComplicated", "Complicated"), direction = "<"), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, print.thres=thr, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("unComplicated", "Complicated"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")
ROCTrain

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Complicated", mode = "everything")
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
ROCTest <- roc(usedTestGroup, ktspStatsTestRes$statistics, plot = T, print.thres=thr$threshold, print.auc=TRUE, print.auc.col="black", ci = T, levels = c("unComplicated", "Complicated"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROCTest

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) >= thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Complicated", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, ConfusionTest$overall["Accuracy"], ConfusionTest$byClass["Balanced Accuracy"], ConfusionTest$byClass["Sensitivity"], ConfusionTest$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
CompKTSP_Perf <- cbind(TrainPerf, TestPerf)

# Save
save(CompKTSP_Perf, file = "./Objs/CompKTSP_Perf.rda")

########################################################################
## Plot AUC of both the agnostic and mechanistic classifiers using ggplot2
### Prepare the legend
forLegend_KTSP <- apply(rbind(
  ci(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("unComplicated", "Complicated"), direction = "<")),
  ci(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("unComplicated", "Complicated"), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


#################################################################
### ROC curves Using ggplot2

### Training
datTrn_KTSP <- melt(data.frame(
  ## Training Group
  Training=factor(usedTrainGroup, levels = c("unComplicated", "Complicated")),
  ## Mechanistic KTSP SUM training
  KTSP.Training=ktspStatsTrainRes$statistics))
### Change Colnames
colnames(datTrn_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")


### Testing
datTst_KTSP <- melt(data.frame(
  ## Testing group
  Testing=factor(usedTestGroup, levels = c("unComplicated", "Complicated")),
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
png("./Figs/AUCggplot_Complicated.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "Performance of the complicated malaria signature in the training and testing data"
#legendTitle <- paste("Mechanistic (", nrow(ktspPredictorRes$TSPs), " pairs)",
#                     "Agnostic (", nrow(ktspPredictorUnRes$TSPs), " pairs)",  sep="")
### Plot
basicplot_KTSP <- ggplot(dat_KTSP, aes(d=Status, m=KTSP_sum, color=KTSP_type,
                                       linetype = KTSP_type)) +
  geom_roc(cutoffs.at = seq(1,20,1)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  #scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  #scale_linetype_manual(legendTitle, values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_KTSP
### Close device
dev.off()

save(basicplot_KTSP, file = "./Objs/BasicPlot_KTSP_Complicated.rda")

##################################################################################################
### Heatmaps
##################################################################################################
#### Make a heatmap of the TSPs expression in the 3 datasets


Comp_Stats <- ktspStatsTestRes$comparisons
Comp_Stats <- Comp_Stats*1

# Order the stats by the sum of votes
NewOrder <- order(rowSums(Comp_Stats))
Comp_Stats <- Comp_Stats[NewOrder, ]

#GroupCol <- ArrayGroup[NewOrder]
#levels(GroupCol) <- c("red", "blue")

# Order the true class labels
CompGroup <- usedTestGroup
CompGroup <- CompGroup[NewOrder]

# Get the predicted class labels
PredClass <- usedTestPredictionRes
#PredClass <- factor(PredClass, levels = c("POS", "NEG")) 
levels(PredClass) <- c("blue", "red")
PredClass <- PredClass[NewOrder]
# Check order 
all(rownames(Comp_Stats) == names(PredClass)) # TRUE:: No Need to order them (Already ordered by sum of votes) 

# Plot
png(filename = "./Figs/CompKTSP_Heatmap.png", width = 3000, height = 2000, res = 200)
superheat(Comp_Stats, col.dendrogram = F, yr = rowSums(Comp_Stats), yr.plot.type  = "bar", left.label.text.col = c("blue", "red"), membership.rows = CompGroup, bottom.label.text.size = 3.5, yr.num.ticks = 6, yr.axis.name	= "Sum of votes", yr.obs.col = PredClass, title = "Heatmap for the TSPs votes in the testing data")
dev.off()

