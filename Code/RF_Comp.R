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
load("./Objs/MalariaDataGood_Comp.rda")
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
tuneRF(x=predictor_data, y=target, plot = TRUE, improve = 0.01, ntreeTry = 1000, proximity = TRUE, sampsize = sampsizes, na.action = na.omit)

## Ntrees = 500, mtry = 400
set.seed(333)  
RF_Agnostic <- randomForest(x =predictor_data, y=target, importance = TRUE, ntree = 1000, mtry = 50, proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
print(RF_Agnostic)
plot(RF_Agnostic)

#save(RF_Agnostic, file = "./Objs/RF_Model.rda")

## RandomForest calculates an importance measures for each variable.
rf_importances <- randomForest::importance(RF_Agnostic, scale=T)
rf_importances <- rf_importances[order(rf_importances[,"MeanDecreaseGini"], decreasing = T), ]
rf_importances <- rf_importances[rf_importances[,"MeanDecreaseGini"] != 0, ]

## Create a representation of the top 30 variables categorized by importance.
png("./Figs/RF_varImp_Comp.png", width = 2000, height = 2000, res = 300)
varImpPlot(RF_Agnostic, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

## An MDS plot provides a sense of the separation of classes.
png("./Figs/MDS_Training_Comp.png", width = 2000, height = 2000, res = 300)
target_labels=as.vector(target)
MDSplot(RF_Agnostic, target, k=2, pch=target_labels, palette=c("red", "blue", "green"), main="MDS plot")
dev.off()


# ROC curve in training data
train_pred_votes_Agnostic <- predict(RF_Agnostic, newdata = predictor_data, type = "vote")
#roc(usedTrainGroup, train_pred_votes_Agnostic[,2], plot = F, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)

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
roc.multiTest <- multiclass.roc(usedTestGroup, RF_predictions_votes_Agnostic, plot = T, print.auc = TRUE, levels = c("control", "unComplicated", "Complicated"), col = "blue", grid = TRUE)
AUC <- roc.multiTest$auc
AUC

# rs <- roc.multiTest[['rocs']]
# plot.roc(rs[[3]][[1]])
# sapply(1:length(rs),function(i) lines.roc(rs[[i]][[i]],col=i))

### MDS plot for testing data
## An MDS plot provides a sense of the separation of classes.
png("./Figs/MDS_Testing_Comp.png", width = 2000, height = 2000, res = 300)
target_labels=as.vector(usedTestGroup)
MDSplot(RF_Agnostic, usedTestGroup, k=2, pch=target_labels, palette=c("red", "blue", "green"), main="MDS plot")
dev.off()

#########
## Using multiROC package
library(multiROC)
library(tidyverse)
library(dummies)
library(ggplot2)


true_label <- dummy(usedTestGroup, sep = ".")
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
final_df <- cbind(true_label, RF_predictions_votes_Agnostic)
colnames(final_df)[c(4,5,6)] <- paste(colnames(final_df)[c(4,5,6)], "_pred_RF")


roc_res <- multi_roc(final_df, force_diag=T)
plot_roc_df <- plot_roc_data(roc_res)

png(filename = "./Figs/MultiROC_TestingData_Comp.png", width = 2000, height = 2000, res = 300)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group, linetype=Method), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))

dev.off()


pr_res <- multi_pr(final_df, force_diag=T)
plot_pr_df <- plot_pr_data(pr_res)

png(filename = "./Figs/MultiPRC_TestingData_Comp.png", width = 2000, height = 2000, res = 300)
ggplot(plot_pr_df, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = Group, linetype=Method), size=1.5) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))
dev.off()


################
## Look at the results closly
unlist(roc_res$AUC)  # Each class vs the others + macro + micro
unlist(pr_res$AUC)

## With CIs
roc_auc_with_ci_res <- roc_auc_with_ci(final_df, conf= 0.95, type='basic', R = 100)
roc_auc_with_ci_res





