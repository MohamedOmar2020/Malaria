##################################
## Test both the severe and cerebral malaria signatures on many infectious diseases

rm(list = ls())

setwd("~/Documents/Research/Projects/Malaria")

renv::activate("~/Documents/Research/Projects/Malaria")

library(GEOquery)
library(randomForest)
library(pROC)
library(precrec)
library(ggplot2)

# ManyInfections1 <- getGEO("GSE42026", GSEMatrix = T, AnnotGPL = T)
# ManyInfections1 <- ManyInfections1$GSE42026_series_matrix.txt.gz
# 
# save(ManyInfections1, file = "./Data/ManyInfections1.rda")

load("./Data/ManyInfections1.rda")

Expr_ManyInfections1 <- exprs(ManyInfections1)
Pheno_ManyInfections1 <- pData(ManyInfections1)
FeatData_ManyInfections1 <- fData(ManyInfections1)


############################
## Annotation

## Expr_ManyInfections1
head(rownames(Expr_ManyInfections1))
rownames(Expr_ManyInfections1) <- FeatData_ManyInfections1$`Gene symbol`
summary(is.na(rownames(Expr_ManyInfections1)))
#rownames(Expr_ManyInfections1) <- gsub("-","", rownames(Expr_ManyInfections1))
#rownames(Expr_ManyInfections1) <- gsub("_","",rownames(Expr_ManyInfections1))
sel <- which(apply(Expr_ManyInfections1, 1, function(x) all(is.finite(x)) ))
Expr_ManyInfections1 <- Expr_ManyInfections1[sel, ]
Expr_ManyInfections1 <- Expr_ManyInfections1[!is.na(rownames(Expr_ManyInfections1)),]
dim(Expr_ManyInfections1)

range(Expr_ManyInfections1)

Expr_ManyInfections1 <- log2(Expr_ManyInfections1 + 21)
Expr_ManyInfections1 <- t(scale(t(Expr_ManyInfections1), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# Asymptomatic vs severe disease

# Pheno1

Pheno_ManyInfections1$DiseaseStatus <- as.factor(Pheno_ManyInfections1$`infecting pathogen:ch1`)
levels(Pheno_ManyInfections1$DiseaseStatus) <- c("case", "case", "control", "case") 
table(Pheno_ManyInfections1$DiseaseStatus)
Pheno_ManyInfections1$DiseaseStatus <- factor(Pheno_ManyInfections1$DiseaseStatus, levels = c("control", "case"))

all(rownames(Pheno_ManyInfections1) == colnames(Expr_ManyInfections1))

ClassInfectionsVsHealthy<- Pheno_ManyInfections1$DiseaseStatus

####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict 

TestingData_ManyInfections <- t(Expr_ManyInfections1)

PredVotes_ManyInfections <- predict(RF_Comp, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections <- predict(RF_Comp, TestingData_ManyInfections, type="response")

ROCTest <- roc(ClassInfectionsVsHealthy, PredVotes_ManyInfections[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_ManyInfections <- evalmod(scores = PredVotes_ManyInfections[,2], labels = ClassInfectionsVsHealthy)
sscurves_ManyInfections
ROC_ManyInfections1 <- autoplot(sscurves_ManyInfections, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE42026 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.49"), size = 4)
PRC_ManyInfections1 <- autoplot(sscurves_ManyInfections, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE42026 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.68"), size = 4)

save(ROC_ManyInfections1, PRC_ManyInfections1, file = "./Objs/ManyInfections1_Curves.rda")

########################################################################################


####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

#################
## Predict

PredVotes_ManyInfections_cerebral <- predict(RF_Cerebral, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections_cerebral <- predict(RF_Cerebral, TestingData_ManyInfections, type="response")

ROCTest <- roc(ClassInfectionsVsHealthy, PredVotes_ManyInfections_cerebral[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_ManyInfections_cerebral <- evalmod(scores = PredVotes_ManyInfections_cerebral[,2], labels = ClassInfectionsVsHealthy)
sscurves_ManyInfections_cerebral
ROC_ManyInfections1_cerebral <- autoplot(sscurves_ManyInfections_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE42026 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.31"), size = 4)
PRC_ManyInfections1_cerebral <- autoplot(sscurves_ManyInfections_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE42026 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.52"), size = 4)

save(ROC_ManyInfections1_cerebral, PRC_ManyInfections1_cerebral, file = "./Objs/ManyInfections1_Curves_cerebral.rda")

#################################
sessionInfo()