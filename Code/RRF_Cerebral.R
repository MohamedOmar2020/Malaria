############################################################################


rm(list = ls())

#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/BreastChemo")

### Load library
library(RRF)
require(limma)
library(randomForest)
library(boot)
library(patchwork)
library(ggThemeAssist)

## Load data
load("./Objs/MalariaDataGood_NCvsC.rda")
#load("./Objs/CerebralExtraValidation.rda")
#load("./Objs/CerebralExtraValidation2.rda")

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
levels(DataTrain[, "usedTrainGroup"]) <- c("nonCerebral", "cerebral")

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
# bootobject_Cerebral <- boot(data = DataTrain, statistic = RF_Strap, R = 100, parallel = "multicore", ncpus = 15) 
# 
# save(bootobject_Cerebral, file = "./Objs/bootobject_Cerebral.rda")

load("./Objs/bootobject_Cerebral.rda")

OutFeat <- bootobject_Cerebral$t

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
# Frequency plot
# sum_result <- sum_result[order(sum_result$rep_rows, decreasing = T), ]
# 
# png(filename = "./Figs/CerebralFrequency.png", width = 2000, height = 2000, res = 300)
# CerebralFreq <- ggplot(data=sum_result, aes(x=rep_rows, y=reorder(Gene, rep_rows))) +
#   geom_col(width=0.5) + 
#   scale_x_continuous(limits = c(0,10), breaks = 0:10) +
#   labs(y = "Gene", x = "Frequency", title = " Frequency of genes in the cerebral malaria signature")
# CerebralFreq
# dev.off()

#save(CerebralFreq, file = "./Objs/CerebralFreqPlot.png")

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
# png(filename = "./Figs/CerebralHeatmap.png", width = 2000, height = 1500, res = 300)
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
# levels(DataTrain$usedTrainGroup) <- c("nonCerebral", "cerebral")
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
# png(filename = "./Figs/CerebralViolinPlot.png", width = 2000, height = 1200, res = 150)
# ggplot(X, 
#        aes(x = usedTrainGroup, 
#            y = Expression)) + 
#   geom_violin(aes(fill = usedTrainGroup),
#               scale = "count")+
#   #geom_jitter(width = 0.1, size = 0.2)+
#   facet_wrap(~Gene)
# dev.off()
# 

# set.seed(333)
# tuneRF(x = PredictorData_Filt, y = usedTrainGroup, mtryStart = 1, ntreeTry = 500, stepFactor = 1, improve = 0.01, trace = F, plot = F)
# 
# set.seed(333)
# RF <- randomForest(usedTrainGroup~., data = DataTrain, mtry = 1, ntree = 500, trace = F, plot = F, doBest = T, sampsize = sampsizes, importance = T)
# RF
# 
# explain_forest(RF, interactions = TRUE, data = DataTrain)

################
## Build the random forest model
set.seed(333)
RF_Cerebral <- tuneRF(x = PredictorData_Filt, y = usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
RF_Cerebral

# treeList <- RF2List(RF_Cerebral)  # transform rf object to an inTrees' format
# exec <- extractRules(treeList, PredictorData_Filt)  # R-executable conditions
# exec[1:2,]
# # 
# ruleMetric <- getRuleMetric(exec,PredictorData_Filt,usedTrainGroup)  # get rule metrics
# ruleMetric[1:2,]
# # 
# readableRules <- presentRules(ruleMetric, colnames(PredictorData_Filt))
# readableRules[1:2, ]
# 
# freqPattern <- getFreqPattern(ruleMetric)

# save the model
#save(RF_Cerebral, file = "./Objs/RF_Cerebral.rda")

load("./Objs/RF_Cerebral.rda")
################
# Predict in the training data
PredVotes_Train <- predict(RF_Cerebral, newdata = PredictorData_Filt, type = "vote")
PredResponse_Train <- predict(RF_Cerebral, PredictorData_Filt, type="response")

ROCTrain <- roc(usedTrainGroup, PredVotes_Train[,2], plot = F, print.auc = TRUE, levels = c("nonCerebral", "cerebral"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTrain

confusion_test <- confusionMatrix(PredResponse_Train, usedTrainGroup, positive = "cerebral")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = PredResponse_Train, actuals = usedTrainGroup)
MCC_Train


#################
## Predict in the testing data
PredVotes_Test <- predict(RF_Cerebral, newdata = TestingData_Filt, type = "vote")
PredResponse_Test <- predict(RF_Cerebral, TestingData_Filt, type="response")

ROCTest <- roc(usedTestGroup, PredVotes_Test[,2], plot = F, print.auc = TRUE, levels = c("nonCerebral", "cerebral"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(PredResponse_Test, usedTestGroup, positive = "cerebral", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = PredResponse_Test, actuals = usedTestGroup)
MCC_Test

# For ROC and PRC curves
sscurves_Test_Cerebral <- evalmod(scores = PredVotes_Test[,2], labels = usedTestGroup)
sscurves_Test_Cerebral
ROC_Test_Cerebral <- autoplot(sscurves_Test_Cerebral, curvetype = c("ROC")) + labs(title = "ROC curve 1st testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.98"), size = 5)
PRC_Test_Cerebral <- autoplot(sscurves_Test_Cerebral, curvetype = c("PRC")) + labs(title = "PRC curve 1st testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.98"), size = 5)


######################################################################
## Predict in the testing data2
TestingData_Filt2 <- t(Expr_Test2[Sel, ])

PredVotes_Test2 <- predict(RF_Cerebral, newdata = TestingData_Filt2, type = "vote")
PredResponse_Test2 <- predict(RF_Cerebral, TestingData_Filt2, type="response")

ROCTest2 <- roc(ClassCerebralVsNonCerebral, PredVotes_Test2[,2], plot = F, print.auc = TRUE, levels = c("nonCerebral", "cerebral"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest2

### Resubstitution performance in the Test set
ConfusionTest2 <- confusionMatrix(PredResponse_Test2, ClassCerebralVsNonCerebral, positive = "cerebral", mode = "everything")
ConfusionTest2

MCC_Test2 <- mltools::mcc(pred = PredResponse_Test2, actuals = ClassCerebralVsNonCerebral)
MCC_Test2

# For ROC and PRC curves
sscurves_Test_Cerebral2 <- evalmod(scores = PredVotes_Test2[,2], labels = ClassCerebralVsNonCerebral)
sscurves_Test_Cerebral2
ROC_Test_Cerebral2 <- autoplot(sscurves_Test_Cerebral2, curvetype = c("ROC")) + labs(title = "ROC curve 2nd testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.87"), size = 5)
PRC_Test_Cerebral2 <- autoplot(sscurves_Test_Cerebral2, curvetype = c("PRC")) + labs(title = "PRC curve 2nd testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.82"), size = 5)

######################################################################
## Predict in the testing data3
TestingData_Filt3 <- t(Expr_Test3[Sel, ])

PredVotes_Test3 <- predict(RF_Cerebral, newdata = TestingData_Filt3, type = "vote")
PredResponse_Test3 <- predict(RF_Cerebral, TestingData_Filt3, type="response")

ROCTest3 <- roc(ClassCerebralVsNonCerebral3, PredVotes_Test3[,2], plot = F, print.auc = TRUE, levels = c("nonCerebral", "cerebral" ), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest3

### Resubstitution performance in the Test set
ConfusionTest3 <- confusionMatrix(PredResponse_Test3, ClassCerebralVsNonCerebral3, positive = "cerebral", mode = "everything")
ConfusionTest3

MCC_Test3 <- mltools::mcc(pred = PredResponse_Test3, actuals = ClassCerebralVsNonCerebral3)
MCC_Test3

# For ROC and PRC curves
sscurves_Test_Cerebral3 <- evalmod(scores = PredVotes_Test3[,2], labels = ClassCerebralVsNonCerebral3)
sscurves_Test_Cerebral3
ROC_Test_Cerebral3 <- autoplot(sscurves_Test_Cerebral3, curvetype = c("ROC")) + labs(title = "ROC curve 3rd testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUC = 1"), size = 5)
PRC_Test_Cerebral3 <- autoplot(sscurves_Test_Cerebral3, curvetype = c("PRC")) + labs(title = "PRC curve 3rd testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 1"), size = 5)

#########################################################################
## Combine testing data
CommonGns <- intersect(rownames(usedTestMat), rownames(Expr_Test2))
CommonGns <- intersect( CommonGns, rownames(Expr_Test3))
usedTestMat <- usedTestMat[CommonGns, ]
Expr_Test2 <- Expr_Test2[CommonGns, ]
Expr_Test3 <- Expr_Test3[CommonGns, ]

Expr_MetaTest <- cbind(usedTestMat, Expr_Test2, Expr_Test3)
Expr_MetaTest <- normalizeBetweenArrays(Expr_MetaTest, method = "quantile")

Pheno_MetaTest <- c(usedTestGroup, ClassCerebralVsNonCerebral, ClassCerebralVsNonCerebral3)
table(Pheno_MetaTest)
Pheno_MetaTest <- as.factor(Pheno_MetaTest)
levels(Pheno_MetaTest) <- c("nonCerebral", "cerebral")

## Predict in the testing Meta data
MetaTestingData_Filt <- t(Expr_MetaTest[Sel, ])

PredVotes_MetaTest <- predict(RF_Cerebral, newdata = MetaTestingData_Filt, type = "vote")
PredResponse_MetaTest <- predict(RF_Cerebral, MetaTestingData_Filt, type="response")

ROC_MetaTest <- roc(Pheno_MetaTest, PredVotes_MetaTest[,2], plot = F, print.auc = TRUE, levels = c("nonCerebral", "cerebral" ), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROC_MetaTest

### Resubstitution performance in the Test set
Confusion_MetaTest <- confusionMatrix(PredResponse_MetaTest, Pheno_MetaTest, positive = "cerebral", mode = "everything")
Confusion_MetaTest

MCC_MetaTest <- mltools::mcc(pred = PredResponse_MetaTest, actuals = Pheno_MetaTest)
MCC_MetaTest

# For ROC and PRC curves
sscurves_MetaTest_Cerebral <- evalmod(scores = PredVotes_MetaTest[,2], labels = Pheno_MetaTest)
sscurves_MetaTest_Cerebral
ROC_MetaTest_Cerebral <- autoplot(sscurves_Test_Cerebral3, curvetype = c("ROC")) + labs(title = "ROC curve combined testing data") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.98"), size = 5)
PRC_MetaTest_Cerebral <- autoplot(sscurves_Test_Cerebral3, curvetype = c("PRC")) + labs(title = "PRC curve combined testing data") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.976"), size = 5)

########################################################################
##############################################
## Make a combined figure for the paper
png(filename = "./Figs/CerebralSignaturesPerformance.png", width = 1500, height = 1000, res = 100)
((ROC_Test_Cerebral / PRC_Test_Cerebral + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
    (ROC_Test_Cerebral2 / PRC_Test_Cerebral2 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
    (ROC_Test_Cerebral3 / PRC_Test_Cerebral3 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
    (ROC_MetaTest_Cerebral / PRC_MetaTest_Cerebral + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12)))
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the cerebral malaria signature in the testing data',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 17, face = "bold"))
  )
dev.off()

###################
