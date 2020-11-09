#################################################################################
### Mohamed Omar
### November 5, 2020
### GOAL: Creating a rondom forest classifier for Malaria (control vs nonCerebral vs Cerebral)
### Using: ALL genes (Agnostic)
#################################################################################

# Clean the work space
rm(list = ls())

## settng the working directory
#setwd("/Users/mohamedomar/Documents/Research/Projects/Bladder")

## Load necessary libraries
library(randomForest)
library(pROC)
library(caret)
library(limma)
library(mltools)


## Load the data
load("./Objs/MalariaDataGood_NCvsC.rda")
#load("./Objs/Correlation/RGenes.rda")



### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

#names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

#names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
#TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 74)

## Subset the expression matrix to the top DE genes only
#usedTrainMat <- usedTrainMat[TopDEgenes, ]
#usedTestMat <- usedTestMat[TopDEgenes, ]
################################################################################
################################################################################
################################################################################

## transpose the matrix
predictor_data <- t(usedTrainMat)
predictor_names <- c(as.vector(rownames(usedTrainMat))) #gene symbol
colnames(predictor_data) <- predictor_names

## Setting the variable we are trying to predict as our target variable. In this case, it is Progression status.
## train group here is just the column containing the phenptype of interest (Progression vs NoProgression) from the phenotype table

#usedTrainGroup <- ordered(usedTrainGroup, levels=c("NoProgression", "Progression"))
target <- usedTrainGroup

## Finally we run the RF algorithm. 
## NOTE: use an ODD number for ntree. This is because when the forest is used on test data, ties are broken randomly. Having an odd number of trees avoids this issue.
## Use down-sampling to attempt to compensate for unequal class-sizes (less progression than noProgression).
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)

##############

# Tunning RF model (to find the best mtry)
set.seed(333)
tuneRF(x=predictor_data, y = target, plot = TRUE, improve = 0.01, ntreeTry = 1000, proximity = TRUE, sampsize = sampsizes, na.action = na.omit)

## Ntrees = 500, mtry = 400
set.seed(333)  
RF_Agnostic <- randomForest(x =predictor_data, y=target, importance = TRUE, ntree = 1000, mtry = 13 ,proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
print(RF_Agnostic)
plot(RF_Agnostic)

#save(RF_Agnostic, file = "./Objs/RF_Model.rda")

## RandomForest calculates an importance measures for each variable.
rf_importances <- randomForest::importance(RF_Agnostic, scale=T)
rf_importances <- rf_importances[order(rf_importances[,"MeanDecreaseGini"], decreasing = T), ]
rf_importances <- rf_importances[rf_importances[,"MeanDecreaseGini"] != 0, ]

## Create a representation of the top 30 variables categorized by importance.
png("./Figs/RF_varImp_NCvsC.png", width = 2000, height = 2000, res = 300)
varImpPlot(RF_Agnostic, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

## An MDS plot provides a sense of the separation of classes.
png("./Figs/MDS_Training_NCvsC.png", width = 2000, height = 2000, res = 300)
target_labels=as.vector(target)
MDSplot(RF_Agnostic, target, k=2, pch=target_labels, palette=c("red", "blue", "green"), main="MDS plot")
dev.off()


# ROC curve in training data
train_pred_votes_Agnostic <- predict(RF_Agnostic, newdata = predictor_data, type = "vote")
roc(usedTrainGroup, train_pred_votes_Agnostic[,2], plot = F, print.auc=TRUE, levels = c("nonCerebral", "cerebral"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Predict in the training data

train_pred_Response_Agnostic <- predict(RF_Agnostic, newdata = predictor_data, type = "response")

confusion_train <- confusionMatrix(train_pred_Response_Agnostic, usedTrainGroup, positive = "cerebral")
confusion_train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = train_pred_Response_Agnostic, actuals = usedTrainGroup)
MCC_Train
################################################################################ 
################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(usedTestMat)
predictor_names2 <- c(as.vector(rownames(usedTestMat))) #gene symbol
colnames(predictor_data2) <- predictor_names2

## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Agnostic$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Agnostic <- predict(RF_Agnostic, predictor_data2, type="response")
RF_predictions_votes_Agnostic <- predict(RF_Agnostic, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Agnostic, usedTestGroup, positive = "cerebral")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Agnostic, actuals = usedTestGroup)
MCC_Test

## ROC curve and AUC
ROCTest <- roc(usedTestGroup, RF_predictions_votes_Agnostic[,2], plot = T, print.auc = TRUE, ci = T, levels = c("nonCerebral", "cerebral"), direction = "<", col = "blue", grid = TRUE)
ROCTest

# rs <- roc.multiTest[['rocs']]
# plot.roc(rs[[3]][[1]])
# sapply(1:length(rs),function(i) lines.roc(rs[[i]][[i]],col=i))

### MDS plot for testing data
## An MDS plot provides a sense of the separation of classes.
png("./Figs/MDS_Testing_NCvsC.png", width = 2000, height = 2000, res = 300)
target_labels=as.vector(usedTestGroup)
MDSplot(RF_Agnostic, usedTestGroup, k=2, pch=target_labels, palette=c("red", "blue", "green"), main="MDS plot")
dev.off()

#######################################################################
######################################################################
### Partial dependence plots


## Training data
DataTrain <- cbind(predictor_data, usedTrainGroup)

# Apparently, CRIPT is overexpressed in nonCerebral vs cerebral 
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "CRIPT", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "CRIPT", which.class = "cerebral")

# same for CDC42SE1
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "CDC42SE1", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "CDC42SE1", which.class = "cerebral")

# Same for EGR1
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "EGR1", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "EGR1", which.class = "cerebral")

# Same
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "NDUFC1", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "NDUFC1", which.class = "cerebral")

# same
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "PPP6C", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "PPP6C", which.class = "cerebral")

# same
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "C18orf8", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "C18orf8", which.class = "cerebral")

# same
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "CTNNB1", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "CTNNB1", which.class = "cerebral")

# same
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "ZFR2", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "ZFR2", which.class = "cerebral")

# ZFR2 is Overexpressed in cerebral vs non cerebral
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "ZFR2", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "ZFR2", which.class = "cerebral")

# same
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "ASB7", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "ASB7", which.class = "cerebral")

# TRIM8 is Overexpressed in cerebral vs non cerebral
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "TRIM8", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTrain, x.var = "TRIM8", which.class = "cerebral")


## Testing data
DataTest <- cbind(predictor_data2, usedTestGroup)

# Apparently, SLC39A8 is overexpressed in nonCerebral and cerebral compared to control
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "SLC39A8", which.class = "control")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "SLC39A8", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "SLC39A8", which.class = "cerebral")

# same for MTHFD2
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "MTHFD2", which.class = "control")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "MTHFD2", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "MTHFD2", which.class = "cerebral")

# BCL11B is overexpressed in control vs cerebral and nonCerebral
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "BCL11B", which.class = "control")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "BCL11B", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "BCL11B", which.class = "cerebral")

# BCL6 is overexpressed in nonCerebral and cerebral compared to control
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "BCL6", which.class = "control")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "BCL6", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "BCL6", which.class = "cerebral")

# MTF1 is overexpressed in nonCerebral and cerebral compared to control
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "MTF1", which.class = "control")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "MTF1", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "MTF1", which.class = "cerebral")

# TLR8 is overexpressed in nonCerebral and cerebral compared to control
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "TLR8", which.class = "control")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "TLR8", which.class = "nonCerebral")
partialPlot(RF_Agnostic, pred.data = DataTest, x.var = "TLR8", which.class = "cerebral")




