rm(list = ls())

setwd("~/Documents/Research/Projects/Malaria")

renv::activate("~/Documents/Research/Projects/Malaria")

library(GEOquery)
library(precrec)
library(ggplot2)
library(randomForest)
library(pROC)

# HIVandTB <- getGEO("GSE39940", GSEMatrix = T, AnnotGPL = T)
# HIVandTB <- HIVandTB$GSE39940_series_matrix.txt.gz
# 
# save(HIVandTB, file = "./Data/HIVandTB.rda")

load("./Data/HIVandTB.rda")

Expr_HIVandTB <- exprs(HIVandTB)
Pheno_HIVandTB <- pData(HIVandTB)
FeatData_HIVandTB <- fData(HIVandTB)

############################
## Annotation

## Expr_HIVandTB
head(rownames(Expr_HIVandTB))
rownames(Expr_HIVandTB) <- FeatData_HIVandTB$`Gene symbol`
summary(is.na(rownames(Expr_HIVandTB)))
#rownames(Expr_HIVandTB) <- gsub("-","", rownames(Expr_HIVandTB))
#rownames(Expr_HIVandTB) <- gsub("_","",rownames(Expr_HIVandTB))
sel <- which(apply(Expr_HIVandTB, 1, function(x) all(is.finite(x)) ))
Expr_HIVandTB <- Expr_HIVandTB[sel, ]
Expr_HIVandTB <- Expr_HIVandTB[!is.na(rownames(Expr_HIVandTB)),]
dim(Expr_HIVandTB)

range(Expr_HIVandTB)
Expr_HIVandTB <- log2(Expr_HIVandTB + 57)
Expr_HIVandTB <- t(scale(t(Expr_HIVandTB), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# Any Hiv positive or active TB or both will be case , if both are negative > control

# Pheno1

Pheno_HIVandTB$DiseaseStatus <- as.factor(Pheno_HIVandTB$`disease status:ch1`)
Pheno_HIVandTB$DiseaseStatus <- ifelse(Pheno_HIVandTB$`disease status:ch1` %in% c("latent TB infection", "other disease") & Pheno_HIVandTB$`hiv status:ch1` == "HIV negative", "control", "case")
table(Pheno_HIVandTB$DiseaseStatus)
Pheno_HIVandTB$DiseaseStatus <- factor(Pheno_HIVandTB$DiseaseStatus, levels = c("control", "case"))

all(rownames(Pheno_HIVandTB) == colnames(Expr_HIVandTB))

ClassHIVandTBVsHealthy<- Pheno_HIVandTB$DiseaseStatus

####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict

TestingData_HIVandTB <- t(Expr_HIVandTB)

PredVotes_HIVandTB <- predict(RF_Comp, newdata = TestingData_HIVandTB, type = "vote")
PredResponse_HIVandTB <- predict(RF_Comp, TestingData_HIVandTB, type="response")

ROCTest <- roc(ClassHIVandTBVsHealthy, PredVotes_HIVandTB[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_HIVandTB <- evalmod(scores = PredVotes_HIVandTB[,2], labels = ClassHIVandTBVsHealthy)
sscurves_HIVandTB
ROC_HIVandTB <- autoplot(sscurves_HIVandTB, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE39940 (HIV +ve or active TB vs healthy)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.52"), size = 4)
PRC_HIVandTB <- autoplot(sscurves_HIVandTB, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE39940 (HIV +ve or active TB vs healthy)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.57"), size = 4)

save(ROC_HIVandTB, PRC_HIVandTB, file = "./Objs/HIVandTB_Curves.rda")

########################################################################################
####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

#################
## Predict
PredVotes_HIVandTB_cerebral <- predict(RF_Cerebral, newdata = TestingData_HIVandTB, type = "vote")
PredResponse_HIVandTB_cerebral <- predict(RF_Cerebral, TestingData_HIVandTB, type="response")

ROCTest_cerebral <- roc(ClassHIVandTBVsHealthy, PredVotes_HIVandTB_cerebral[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest_cerebral

# For ROC and PRC curves
sscurves_HIVandTB_cerebral <- evalmod(scores = PredVotes_HIVandTB_cerebral[,2], labels = ClassHIVandTBVsHealthy)
sscurves_HIVandTB_cerebral
ROC_HIVandTB_cerebral <- autoplot(sscurves_HIVandTB_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE39940 (HIV +ve or active TB vs healthy)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.38"), size = 4)
PRC_HIVandTB_cerebral <- autoplot(sscurves_HIVandTB_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE39940 (HIV +ve or active TB vs healthy)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.47"), size = 4)

save(ROC_HIVandTB_cerebral, PRC_HIVandTB_cerebral, file = "./Objs/HIVandTB_Curves_cerebral.rda")

#####################
sessionInfo()
