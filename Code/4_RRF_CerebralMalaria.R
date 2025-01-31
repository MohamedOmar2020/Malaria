############################################################################


rm(list = ls())

setwd("~/Documents/Research/Projects/Malaria")

renv::activate("~/Documents/Research/Projects/Malaria")

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
library(patchwork)
library(mltools)
library(caret)

## Load data
load("./Objs/MalariaDataGood_NCvsC.rda")
load("./Objs/PlacentalMalaria.rda")
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
#lambda <- 0.8 # Both the number of features and the quality of the features are quite sensitive to lambda for RRF. A smaller lambda leads to fewer features.

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

Bootstrap_CerebralMalaria <- OutFeat
save(Bootstrap_CerebralMalaria, file = "Bootstrap_CerebralMalaria.rda")

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
sum_result <- sum_result[order(sum_result$rep_rows, decreasing = T), ]

# png(filename = "./Figs/CerebralFrequency.png", width = 2000, height = 2000, res = 300)
# CerebralFreq <- ggplot(data=sum_result, aes(x=rep_rows, y=reorder(Gene, rep_rows))) +
#   geom_col(width=0.5) +
#   scale_x_continuous(limits = c(0,12), breaks = 0:10) +
#   labs(y = "Gene", x = "Frequency", title = " Frequency of genes in the cerebral malaria signature")
# CerebralFreq
# dev.off()

#save(CerebralFreq, file = "./Objs/CerebralFreqPlot.png")

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
################
## Build the random forest model
set.seed(333)
RF_Cerebral <- tuneRF(x = PredictorData_Filt, y = usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
RF_Cerebral

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
ROC_Test_Cerebral <- autoplot(sscurves_Test_Cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.98"), size = 3)
PRC_Test_Cerebral <- autoplot(sscurves_Test_Cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.98"), size = 3)

########################################################################
##############################################
## Make a combined figure for the paper
load("./Objs/SevereSigROC_PRC.rda")

ROC_Test_Comp$theme$plot.title$size <- 8
PRC_Test_Comp$theme$plot.title$size <- 8
ROC_Test_Cerebral$theme$plot.title$size <- 8
PRC_Test_Cerebral$theme$plot.title$size <- 8


tiff(filename = "./Figs/TwoSignaturesPerformance.tiff", width = 2500, height = 2000, res = 350)
((ROC_Test_Comp / PRC_Test_Comp + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12))) | 
    (ROC_Test_Cerebral / PRC_Test_Cerebral + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 12)))  
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the two malaria signatures in the testing data',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )
dev.off()

###################

sessionInfo()
