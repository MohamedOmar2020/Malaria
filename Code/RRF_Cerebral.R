############################################################################


rm(list = ls())

#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/BreastChemo")

### Load library
library(RRF)
require(limma)
library(randomForest)
library(boot)

## Load data
load("./Objs/MalariaDataGood_NCvsC.rda")

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
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  rrf <- RRF(usedTrainGroup~., data = d, flagReg = 1, coefReg=lambda) # coefReg is a constant for all variables.   #either "X,as.factor(class)" or data frame like "Y~., data=data" is fine, but the later one is significantly slower. 
  TrainData <- d
  TrainData$usedTrainGroup <- NULL
  subsetRRF <- rrf$feaSet
  SelFeats <- colnames(TrainData[, subsetRRF])
  return(as.vector(SelFeats[1:30]))
}

set.seed(333)
bootobject_Cerebral <- boot(data = DataTrain, statistic = RF_Strap, R = 100, parallel = "multicore", ncpus = 15) 

save(bootobject_Cerebral, file = "./Objs/bootobject_Cerebral.rda")

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
RF_Cerebral <- tuneRF(x = PredictorData_Filt, y = usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
RF_Cerebral

# save the model
save(RF_Cerebral, file = "./Objs/RF_Cerebral.rda")
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

ROCTest <- roc(usedTestGroup, PredVotes_Test[,2], plot = F, print.auc = TRUE, levels = c("nonCerebral", "cerebral" ), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(PredResponse_Test, usedTestGroup, positive = "cerebral", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = PredResponse_Test, actuals = usedTestGroup)
MCC_Test

###################3
# RandomForest calculates an importance measures for each variable.
RF_Cerebral_importances <- randomForest::importance(RF_Cerebral, scale=F)
# RF_Cerebral_importances <- RF_Cerebral_importances[order(RF_Cerebral_importances[,4], decreasing = TRUE), ]
# RF_Cerebral_importances <- rf_importances[!(rf_importances[,4] == 0), ]
# #MechanisticRF_Pairs <- rownames(rf_importances)[1:50]
# #save(MechanisticRF_Pairs, file = "./Objs/RF/MechanisticRF_Pairs.rda")

