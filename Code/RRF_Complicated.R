############################################################################


rm(list = ls())

#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/BreastChemo")

### Load library
library(RRF)
require(limma)
library(randomForest)
library(boot)
library(precrec)
library(pheatmap)
library(randomForestExplainer)
library(inTrees)
library(pROC)
library(caret)
library(patchwork)

## Load data
load("./Objs/MalariaDataGood_Comp.rda")
load("./Objs/PlacentalMalaria.rda")
load("./Objs/ExtraMalaria.rda")

# Quantile normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

# Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

# Transpose the matrix
PredictorData <- t(usedTrainMat)

# Bind expression with the groups
DataTrain <- cbind(PredictorData, usedTrainGroup)
DataTrain <- as.data.frame(DataTrain)
DataTrain$usedTrainGroup <- as.factor(DataTrain$usedTrainGroup)
levels(DataTrain[, "usedTrainGroup"]) <- c("unComplicated", "Complicated")

names(DataTrain) <- make.names(names(DataTrain))

# Feature selection via RRF
lambda <- 0.8 # Both the number of features and the quality of the features are quite sensitive to lambda for RRF. A smaller lambda leads to fewer features.

# The function for bootstraping
# RF_Strap <- function(data, indices) {
#   d <- data[indices, ] # allows boot to select sample
#   rrf <- RRF(usedTrainGroup~., data = d, flagReg = 1, coefReg=lambda) # coefReg is a constant for all variables.   #either "X,as.factor(class)" or data frame like "Y~., data=data" is fine, but the later one is significantly slower. 
#   TrainData <- d
#   TrainData$usedTrainGroup <- NULL
#   subsetRRF <- rrf$feaSet
#   SelFeats <- colnames(TrainData[, subsetRRF])
#   return(as.vector(SelFeats[1:30]))
# }
# 
# set.seed(333)
# bootobject_Comp <- boot(data = DataTrain, statistic = RF_Strap, R = 100, parallel = "multicore", ncpus = 15) 
# 
# 
# save(bootobject_Comp, file = "./Objs/bootobject_Comp.rda")

load("./Objs/bootobject_Comp.rda")

OutFeat <- bootobject_Comp$t

#####################################
## Select the most frequent features
u_genes <- na.omit(unique(as.vector(as.matrix(OutFeat))))
find_gen_rep <- function(dat, gene){
  # NAs create problems in the function so we substitute that with "unknown"
  dat[is.na(dat)] <- "unknown"
  rep_rows <- sum(apply(dat, 1,  function(x) any(x == gene)))
  names(rep_rows) <- gene
  as.data.frame(rep_rows)
}
list_results <- lapply(u_genes, find_gen_rep, dat = OutFeat)
sum_result <- do.call(rbind, list_results)

sum_result$Gene <- rownames(sum_result)
sum_result <- sum_result[sum_result$rep_rows >= 5, ]
Sel <- sum_result$Gene
Sel

####################################
# Frequency figure
sum_result <- sum_result[order(sum_result$rep_rows, decreasing = T), ]

png(filename = "./Figs/CompFrequency.png", width = 2000, height = 2000, res = 300)
CompFreq <- ggplot(data=sum_result, aes(x=rep_rows, y=reorder(Gene, rep_rows))) +
  geom_col(width=0.5) +
  scale_x_continuous(limits = c(0,15), breaks = 0:15) +
  labs(y = "Gene", x = "Frequency", title = "Severe malaria signature")
CompFreq
dev.off()


##################
# Together with cerebral signature
load("./Objs/CerebralFreqPlot.png")

CerebralFreq$labels$title <- "Cerebral malaria signature"
tiff(filename = "./Figs/CombinedFrequency.tiff", width = 3000, height = 1700, res = 350)
CompFreq + CerebralFreq
dev.off()


########################################
## Heatmap
# NewOrder <- order(usedTrainGroup)
# usedTrainGroup_ord <- usedTrainGroup[NewOrder]
# 
# X <- usedTrainMat[Sel, ]
# 
# X <- X[, NewOrder]
# 
# Annot <- as.data.frame(usedTrainGroup_ord)
# rownames(Annot) <- names(usedTrainGroup_ord)
# 
# png(filename = "./Figs/CompHeatmap.png", width = 2000, height = 1500, res = 300)
# pheatmap::pheatmap(X, annotation_col = Annot, cluster_cols = F, cluster_rows = F, show_colnames = F)
# dev.off()
#####################################
## Use the selected features to build a new random forest model

usedTrainMat_Filt <- usedTrainMat[Sel, ]

usedTestMat_Filt <- usedTestMat[Sel, ]

# Transpose the matrices
PredictorData_Filt <- t(usedTrainMat_Filt)
TestingData_Filt <- t(usedTestMat_Filt)

### Sampsize
tmp <- as.vector(table(usedTrainGroup))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)

###########################
## Explain the RF

# DataTrain <- cbind(PredictorData_Filt, usedTrainGroup)
# DataTrain <- as.data.frame(DataTrain)
# DataTrain$usedTrainGroup <- as.factor(DataTrain$usedTrainGroup)
# levels(DataTrain$usedTrainGroup) <- c("unComplicated", "Complicated")
# 
# # MAke a violin plot
# library(tidyr)
# 
# X <- pivot_longer(
#   DataTrain,
#   cols = 1:28,
#   names_to = "Gene",
#   names_repair = "check_unique",
#   values_to = "Expression",
# )
# 
# png(filename = "./Figs/ComViolinPlot.png", width = 2000, height = 1200, res = 150)
# ggplot(X, 
#        aes(x = usedTrainGroup, 
#            y = Expression)) + 
#   geom_violin(aes(fill = usedTrainGroup),
#               scale = "count")+
#   #geom_jitter(width = 0.1, size = 0.2)+
#   facet_wrap(~Gene)
# dev.off()


# set.seed(333)
# tuneRF(x = PredictorData_Filt, y = usedTrainGroup, mtryStart = 1, ntreeTry = 500, stepFactor = 1, improve = 0.01, trace = F, plot = F)
# 
# set.seed(333)
# RF <- randomForest(usedTrainGroup~., data = DataTrain, mtry = 1, ntree = 500, trace = F, plot = F, doBest = T, sampsize = sampsizes, importance = T)
# RF
# 
# explain_forest(RF, interactions = TRUE, data = DataTrain)


##################################
## Build the random forest model
set.seed(333)
RF_Comp <- tuneRF(x = PredictorData_Filt, y = usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
RF_Comp

# treeList <- RF2List(RF_Comp)  # transform rf object to an inTrees' format
# exec <- extractRules(treeList, PredictorData_Filt)  # R-executable conditions
# exec[1:2,]
# 
# ruleMetric <- getRuleMetric(exec,PredictorData_Filt,usedTrainGroup)  # get rule metrics
# ruleMetric[1:2,]
# 
# readableRules <- presentRules(ruleMetric, colnames(PredictorData_Filt))
# readableRules[1:2, ]
# 
# freqPattern <- getFreqPattern(ruleMetric)

# Save the model
#save(RF_Comp, file = "./Objs/RF_Comp.rda")

load("./Objs/RF_Comp.rda")
################
# Predict in the training data
PredVotes_Train <- predict(RF_Comp, newdata = PredictorData_Filt, type = "vote")
PredResponse_Train <- predict(RF_Comp, PredictorData_Filt, type="response")

ROCTrain <- roc(usedTrainGroup, PredVotes_Train[,2], plot = F, print.auc = TRUE, levels = c("unComplicated", "Complicated"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTrain

confusion_test <- confusionMatrix(PredResponse_Train, usedTrainGroup, positive = "Complicated")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = PredResponse_Train, actuals = usedTrainGroup)
MCC_Train


#################
## Predict in the testing data
PredVotes_Test <- predict(RF_Comp, newdata = TestingData_Filt, type = "vote")
PredResponse_Test <- predict(RF_Comp, TestingData_Filt, type="response")

ROCTest <- roc(usedTestGroup, PredVotes_Test[,2], plot = F, print.auc = TRUE, levels = c("unComplicated", "Complicated"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

### Resubstitution peRF_Compormance in the Test set
ConfusionTest <- confusionMatrix(PredResponse_Test, usedTestGroup, positive = "Complicated", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = PredResponse_Test, actuals = usedTestGroup)
MCC_Test

# For ROC and PRC curves
sscurves_Test_Comp <- evalmod(scores = PredVotes_Test[,2], labels = usedTestGroup)
sscurves_Test_Comp
ROC_Test_Comp <- autoplot(sscurves_Test_Comp, curvetype = c("ROC")) + labs(title = "ROC curve of the severe malaria signature") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.85"), size = 4)
PRC_Test_Comp <- autoplot(sscurves_Test_Comp, curvetype = c("PRC")) + labs(title = "PRC curve of the severe malaria signature") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.95"), size = 4)

save(ROC_Test_Comp, PRC_Test_Comp, file = "./Objs/SevereSigROC_PRC.rda")
#######################################################################################
## Predict in the 2nd testing data 
# Placental malaria vs non-placental malaria
usedTestMat_Filt2 <- Expr_Test2[Sel, ]
TestingData_Filt2 <- t(usedTestMat_Filt2)

PredVotes_Test2 <- predict(RF_Comp, newdata = TestingData_Filt2, type = "vote")
PredResponse_Test2 <- predict(RF_Comp, TestingData_Filt2, type="response")

ROCTest2 <- roc(ClassComplicatedVSunComplicated, PredVotes_Test2[,2], plot = F, print.auc = TRUE, levels = c("unComplicated", "Complicated"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest2

### Resubstitution peRF_Compormance in the Test set
ConfusionTest2 <- confusionMatrix(PredResponse_Test2, ClassComplicatedVSunComplicated, positive = "Complicated", mode = "everything")
ConfusionTest2

MCC_Test2 <- mltools::mcc(pred = PredResponse_Test2, actuals = ClassComplicatedVSunComplicated)
MCC_Test2

# For ROC and PRC curves
sscurves_PM <- evalmod(scores = PredVotes_Test2[,2], labels = ClassComplicatedVSunComplicated)
sscurves_PM
ROC_PM<- autoplot(sscurves_PM, curvetype = c("ROC")) + labs(title = "ROC curve PM +ve vs PM -ve") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.70"), size = 4)
PRC_PM <- autoplot(sscurves_PM, curvetype = c("PRC")) + labs(title = "PRC curve PM +ve vs PM -ve") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.62"), size = 4)

#######################################################################################
## Predict in the 2nd testing data 
# Inflammation vs no-inflammation
usedTestMat_Filt2 <- Expr_Test2[Sel, ]
TestingData_Filt2 <- t(usedTestMat_Filt2)

PredVotes_Test2 <- predict(RF_Comp, newdata = TestingData_Filt2, type = "vote")
PredResponse_Test2 <- predict(RF_Comp, TestingData_Filt2, type="response")

ROCTest2 <- roc(ClassInflammation, PredVotes_Test2[,2], plot = F, print.auc = TRUE, levels = c("No", "Yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest2

### Resubstitution peRF_Compormance in the Test set
levels(ClassInflammation) <- c("unComplicated", "Complicated")
ConfusionTest2 <- confusionMatrix(PredResponse_Test2, ClassInflammation, positive = "Complicated", mode = "everything")
ConfusionTest2

MCC_Test2 <- mltools::mcc(pred = PredResponse_Test2, actuals = ClassInflammation)
MCC_Test2

# For ROC and PRC curves
sscurves_PM_Inflamm <- evalmod(scores = PredVotes_Test2[,2], labels = ClassInflammation)
sscurves_PM_Inflamm
ROC_PM_Inflamm <- autoplot(sscurves_PM_Inflamm, curvetype = c("ROC")) + labs(title = "ROC curve inflammation vs no inflammation") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.76"), size = 4)
PRC_PM_Inflamm <- autoplot(sscurves_PM_Inflamm, curvetype = c("PRC")) + labs(title = "PRC curve inflammation vs no inflammation") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.60"), size = 4)

ROC_PM$theme$plot.title$size <- 8
PRC_PM$theme$plot.title$size <- 8
ROC_PM_Inflamm$theme$plot.title$size <- 8
PRC_PM_Inflamm$theme$plot.title$size <- 8

#################################################################
#######################################################################################
## Predict in the 2nd testing data 
# mild vs severe malaria
# usedTestMat_Filt3 <- Expr_Test3[Sel, ]
# TestingData_Filt3 <- t(usedTestMat_Filt3)
# 
# PredVotes_Test3 <- predict(RF_Comp, newdata = TestingData_Filt3, type = "vote")
# PredResponse_Test3 <- predict(RF_Comp, TestingData_Filt3, type="response")
# 
# ROCTest3 <- roc(ClassComplicatedVSunComplicated2, PredVotes_Test3[,2], plot = F, print.auc = TRUE, levels = c("unComplicated", "Complicated"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
# ROCTest3
# 
# ### Resubstitution peRF_Compormance in the Test set
# ConfusionTest3 <- confusionMatrix(PredResponse_Test3, ClassComplicatedVSunComplicated2, positive = "Complicated", mode = "everything")
# ConfusionTest3
# 
# MCC_Test3 <- mltools::mcc(pred = PredResponse_Test3, actuals = ClassInflammation)
# MCC_Test3
# 
# # For ROC and PRC curves
# sscurves_PM_Inflamm <- evalmod(scores = PredVotes_Test2[,2], labels = ClassInflammation)
# sscurves_PM_Inflamm
# ROC_PM_Inflamm <- autoplot(sscurves_PM_Inflamm, curvetype = c("ROC")) + labs(title = "ROC curve inflammation vs no inflammation") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.76"), size = 5)
# PRC_PM_Inflamm <- autoplot(sscurves_PM_Inflamm, curvetype = c("PRC")) + labs(title = "PRC curve inflammation vs no inflammation") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.60"), size = 5)

#######################################################################################
## Make a combined figure for the paper
tiff(filename = "./Figs/severeSignaturesPerformance_PM.tiff", width = 2500, height = 2000, res = 350)
((ROC_PM / PRC_PM + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
    (ROC_PM_Inflamm / PRC_PM_Inflamm + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12)))
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the severe malaria signature in the placental malaria dataset',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )
dev.off()

###################
###########################

