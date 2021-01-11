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

## Load data
load("./Objs/MalariaDataGood_Comp.rda")
#load("./Objs/CompExtraValidation.rda")
#load("./Objs/CompExtraValidation2.rda")

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
  scale_x_continuous(limits = c(0,10), breaks = 0:10) +
  labs(y = "Gene", x = "Frequency", title = " Frequency of genes in the severe malaria signature")
CompFreq
dev.off()


##################
# Together with cerebral signature
load("./Objs/CerebralFreqPlot.png")

png(filename = "./Figs/CombinedFrequency.png", width = 2000, height = 2000, res = 300)
CompFreq + CerebralFreq
dev.off()


########################################
## Heatmap
NewOrder <- order(usedTrainGroup)
usedTrainGroup_ord <- usedTrainGroup[NewOrder]

X <- usedTrainMat[Sel, ]

X <- X[, NewOrder]

Annot <- as.data.frame(usedTrainGroup_ord)
rownames(Annot) <- names(usedTrainGroup_ord)

png(filename = "./Figs/CompHeatmap.png", width = 2000, height = 1500, res = 300)
pheatmap::pheatmap(X, annotation_col = Annot, cluster_cols = F, cluster_rows = F, show_colnames = F)
dev.off()
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

################
## Build the random forest model
set.seed(333)
RF_Comp <- tuneRF(x = PredictorData_Filt, y = usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
RF_Comp

# Save the model
save(RF_Comp, file = "./Objs/RF_Comp.rda")

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
ROC_Test_Comp <- autoplot(sscurves_Test_Comp, curvetype = c("ROC")) + labs(title = "ROC curve 1st testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.85"), size = 5)
PRC_Test_Comp <- autoplot(sscurves_Test_Comp, curvetype = c("PRC")) + labs(title = "PRC curve 1st testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.95"), size = 5)

#######################################################################################
######################################################################
## Predict in the testing data2
TestingData_Filt2 <- t(Expr_Test2[Sel, ])

PredVotes_Test2 <- predict(RF_Comp, newdata = TestingData_Filt2, type = "vote")
PredResponse_Test2 <- predict(RF_Comp, TestingData_Filt2, type="response")

ROCTest2 <- roc(ClassComplicatedVSunComplicated, PredVotes_Test2[,2], plot = F, print.auc = TRUE, levels = c("unComplicated", "Complicated"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest2

### Resubstitution performance in the Test set
ConfusionTest2<- confusionMatrix(PredResponse_Test2, ClassComplicatedVSunComplicated, positive = "Complicated", mode = "everything")
ConfusionTest2

MCC_Test2 <- mltools::mcc(pred = PredResponse_Test2, actuals = ClassComplicatedVSunComplicated)
MCC_Test2

# For ROC and PRC curves
sscurves_Test_Comp2 <- evalmod(scores = PredVotes_Test2[,2], labels = ClassComplicatedVSunComplicated)
sscurves_Test_Comp2
ROC_Test_Comp2 <- autoplot(sscurves_Test_Comp2, curvetype = c("ROC")) + labs(title = "ROC curve 2nd testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.85"), size = 5)
PRC_Test_Comp2 <- autoplot(sscurves_Test_Comp2, curvetype = c("PRC")) + labs(title = "PRC curve 2nd testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.90"), size = 5)

######################################################################
## Predict in the testing data3
TestingData_Filt3 <- t(Expr_Test3[Sel, ])

PredVotes_Test3 <- predict(RF_Comp, newdata = TestingData_Filt3, type = "vote")
PredResponse_Test3 <- predict(RF_Comp, TestingData_Filt3, type="response")

ROCTest3 <- roc(ClassComplicatedVSunComplicated3, PredVotes_Test3[,2], plot = F, print.auc = TRUE, levels = c("unComplicated", "Complicated"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest3

### Resubstitution performance in the Test set
ConfusionTest3 <- confusionMatrix(PredResponse_Test3, ClassComplicatedVSunComplicated3, positive = "Complicated", mode = "everything")
ConfusionTest3

MCC_Test3 <- mltools::mcc(pred = PredResponse_Test3, actuals = ClassComplicatedVSunComplicated3)
MCC_Test3

# For ROC and PRC curves
sscurves_Test_Comp3 <- evalmod(scores = PredVotes_Test3[,2], labels = ClassComplicatedVSunComplicated3)
sscurves_Test_Comp3
ROC_Test_Comp3 <- autoplot(sscurves_Test_Comp3, curvetype = c("ROC")) + labs(title = "ROC curve 3rd testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUC = 1"), size = 5)
PRC_Test_Comp3 <- autoplot(sscurves_Test_Comp3, curvetype = c("PRC")) + labs(title = "PRC curve 3rd testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 1"), size = 5)

#########################################################################
## Combine testing data
CommonGns <- intersect(rownames(usedTestMat), rownames(Expr_Test2))
CommonGns <- intersect( CommonGns, rownames(Expr_Test3))
usedTestMat <- usedTestMat[CommonGns, ]
Expr_Test2 <- Expr_Test2[CommonGns, ]
Expr_Test3 <- Expr_Test3[CommonGns, ]

Expr_MetaTest <- cbind(usedTestMat, Expr_Test2, Expr_Test3)
Expr_MetaTest <- normalizeBetweenArrays(Expr_MetaTest, method = "quantile")

Pheno_MetaTest <- c(usedTestGroup, ClassComplicatedVSunComplicated, ClassComplicatedVSunComplicated3)
table(Pheno_MetaTest)
Pheno_MetaTest <- as.factor(Pheno_MetaTest)
levels(Pheno_MetaTest) <- c("unComplicated", "Complicated")

## Predict in the testing Meta data
MetaTestingData_Filt <- t(Expr_MetaTest[Sel, ])

PredVotes_MetaTest <- predict(RF_Comp, newdata = MetaTestingData_Filt, type = "vote")
PredResponse_MetaTest <- predict(RF_Comp, MetaTestingData_Filt, type="response")

ROC_MetaTest <- roc(Pheno_MetaTest, PredVotes_MetaTest[,2], plot = F, print.auc = TRUE, levels = c("unComplicated", "Complicated"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROC_MetaTest

### Resubstitution performance in the Test set
Confusion_MetaTest <- confusionMatrix(PredResponse_MetaTest, Pheno_MetaTest, positive = "Complicated", mode = "everything")
Confusion_MetaTest

MCC_MetaTest <- mltools::mcc(pred = PredResponse_MetaTest, actuals = Pheno_MetaTest)
MCC_MetaTest

# For ROC and PRC curves
sscurves_MetaTest_Comp <- evalmod(scores = PredVotes_MetaTest[,2], labels = Pheno_MetaTest)
sscurves_MetaTest_Comp
ROC_MetaTest_Comp <- autoplot(sscurves_MetaTest_Comp, curvetype = c("ROC")) + labs(title = "ROC curve combined testing data") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.93"), size = 5)
PRC_MetaTest_Comp <- autoplot(sscurves_MetaTest_Comp, curvetype = c("PRC")) + labs(title = "PRC curve combined testing data") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.97"), size = 5)


#######################################################################################
## Make a combined figure for the paper
png(filename = "./Figs/CompSignaturesPerformance.png", width = 1500, height = 1000, res = 100)
((ROC_Test_Comp / PRC_Test_Comp + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
    (ROC_Test_Comp2 / PRC_Test_Comp2 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
    (ROC_Test_Comp3 / PRC_Test_Comp3 + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
    (ROC_MetaTest_Comp / PRC_MetaTest_Comp + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12)))
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the complicated malaria signature in the testing data',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 17, face = "bold"))
  )
dev.off()

###################3
###########################

