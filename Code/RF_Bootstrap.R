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
library(boot)

## Load the data
load("./Objs/MalariaDataGood.rda")
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

DataTrain <- cbind(predictor_data, usedTrainGroup)
DataTrain <- as.data.frame(DataTrain)
DataTrain$usedTrainGroup <- as.factor(DataTrain$usedTrainGroup)
levels(DataTrain[, "usedTrainGroup"]) <- c("control", "nonCerebral", "cerebral")

predictor_data_Test <- t(usedTestMat)

colnames(DataTrain) <- make.names(colnames(DataTrain))
colnames(predictor_data_Test) <- make.names(colnames(predictor_data_Test))

##############
# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- randomForest(usedTrainGroup~., data=d, importance = TRUE, ntree = 500, mtry = 25 ,proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
  Importance <- randomForest::importance(RF, scale = T)
  Importance <- Importance[order(Importance[,"MeanDecreaseGini"], decreasing = TRUE), ]
  Importance <- Importance[Importance[,"MeanDecreaseGini"] > 0, ]
  ImportanVariables <- rownames(Importance)[1:50]
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")  
  test_preds <- predict(RF, newdata = predictor_data_Test, type = "vote")
  roc.multiTrain <- multiclass.roc(PhenoTrain, train_preds, plot = F, print.auc = TRUE, levels = c("control", "nonCerebral", "cerebral"), col = "blue", grid = TRUE)
  roc.multiTest <- multiclass.roc(usedTestGroup, test_preds, plot = F, print.auc = TRUE, levels = c("control", "nonCerebral", "cerebral"), col = "blue", grid = TRUE)
  return(c(roc.multiTrain$auc, roc.multiTest$auc, ImportanVariables))
}


set.seed(333)
bootobjectAgnostic <- boot(data= DataTrain, statistic= RF_Strap, R= 100) 

#save(bootobjectAgnostic, file = "./Objs/bootobjectAgnostic.rda")

load("./Objs/bootobjectAgnostic.rda")

AUCs_RF <- bootobjectAgnostic$t
#colnames(AUCs_RF)[1:2] <- c("AUC_Train", "AUC_Test")

AUCs_RF_Train <- data.frame(AUC = AUCs_RF[, 1])
AUCs_RF_Train$DataType <- "Training"
AUCs_RF_Test <- data.frame(AUC = AUCs_RF[, 2])
AUCs_RF_Test$DataType <- "Testing"

AUC_RF_All <- rbind(AUCs_RF_Train, AUCs_RF_Test)
AUC_RF_All$AUC <- as.numeric(AUC_RF_All$AUC)

png("./Figs/RF_BS_AUC_Testing.png", width = 3000, height = 1500, res = 300)
AUC_Test_DistrHist <- ggplot(AUC_RF_All, aes(AUC, fill = DataType)) + 
  geom_density(alpha = 0.5, adjust = 0.01) +
  scale_x_continuous(limits = c(0.8, 1)) +
  labs(title="AUC distribution of the RF models in the testing data") 
AUC_Test_DistrHist
dev.off()











