##################################
## Test both the severe and cerebral malaria signatures on tuberculosis

rm(list = ls())

setwd("~/Documents/Research/Projects/Malaria")

renv::activate("~/Documents/Research/Projects/Malaria")

library(GEOquery)
library(randomForest)
library(pROC)
library(precrec)
library(ggplot2)

# TB3 <- getGEO("GSE62525", GSEMatrix = T, AnnotGPL = T)
# TB3 <- TB3$GSE62525_series_matrix.txt.gz
# 
# save(TB3, file = "./Data/TB3.rda")

load("./Data/TB3.rda")

Expr_TB3 <- exprs(TB3)
Pheno_TB3 <- pData(TB3)
FeatData_TB3 <- fData(TB3)

############################
## Annotation

## Expr_TB3
head(rownames(Expr_TB3))
rownames(Expr_TB3) <- FeatData_TB3$Gene_symbol
summary(is.na(rownames(Expr_TB3)))
#rownames(Expr_TB3) <- gsub("-","", rownames(Expr_TB3))
#rownames(Expr_TB3) <- gsub("_","",rownames(Expr_TB3))
sel <- which(apply(Expr_TB3, 1, function(x) all(is.finite(x)) ))
Expr_TB3 <- Expr_TB3[sel, ]
Expr_TB3 <- Expr_TB3[!is.na(rownames(Expr_TB3)),]
dim(Expr_TB3)

range(Expr_TB3)
Expr_TB3 <- log2(Expr_TB3 + 21)
Expr_TB3 <- t(scale(t(Expr_TB3), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# primary TB and pneumonia vs latent TB

# Pheno1

Pheno_TB3$DiseaseStatus <- as.factor(Pheno_TB3$`disease state:ch1`)
levels(Pheno_TB3$DiseaseStatus) <- c("case", "control", "control")
Pheno_TB3$DiseaseStatus <- factor(Pheno_TB3$DiseaseStatus, levels = c("control", "case"))
table(Pheno_TB3$DiseaseStatus)
all(rownames(Pheno_TB3) == colnames(Expr_TB3))

ClassTBVsLatentTBandHealthy<- Pheno_TB3$DiseaseStatus

####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict

TestingData_TB3 <- t(Expr_TB3)

PredVotes_TB3 <- predict(RF_Comp, newdata = TestingData_TB3, type = "vote")
PredResponse_TB3 <- predict(RF_Comp, TestingData_TB3, type="response")

ROCTest <- roc(ClassTBVsLatentTBandHealthy, PredVotes_TB3[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_TB3 <- evalmod(scores = PredVotes_TB3[,2], labels = ClassTBVsLatentTBandHealthy)
sscurves_TB3
ROC_TB3 <- autoplot(sscurves_TB3, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE62525 (primary TB vs latent TB and healthy)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.566"), size = 4)
PRC_TB3 <- autoplot(sscurves_TB3, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE62525 (primary TB vs latent TB and healthy)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.46"), size = 4)

save(ROC_TB3, PRC_TB3, file = "./Objs/TB3_Curves.rda")

########################################################################################
####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

#################
## Predict
PredVotes_TB3_cerebral <- predict(RF_Cerebral, newdata = TestingData_TB3, type = "vote")
PredResponse_TB3_cerebral <- predict(RF_Cerebral, TestingData_TB3, type="response")

ROCTest_cerebral <- roc(ClassTBVsLatentTBandHealthy, PredVotes_TB3_cerebral[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest_cerebral

# For ROC and PRC curves
sscurves_TB3_cerebral <- evalmod(scores = PredVotes_TB3_cerebral[,2], labels = ClassTBVsLatentTBandHealthy)
sscurves_TB3_cerebral
ROC_TB3_cerebral <- autoplot(sscurves_TB3_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE62525 (primary TB vs latent TB and healthy)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.15"), size = 4)
PRC_TB3_cerebral <- autoplot(sscurves_TB3_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE62525 (primary TB vs latent TB and healthy)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.21"), size = 4)

save(ROC_TB3_cerebral, PRC_TB3_cerebral, file = "./Objs/TB3_Curves_cerebral.rda")

###################
sessionInfo()

