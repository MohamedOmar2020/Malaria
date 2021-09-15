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

# ManyInfections2 <- getGEO("GSE6269", GSEMatrix = T, AnnotGPL = T)
# ManyInfections2 <- ManyInfections2$`GSE6269-GPL96_series_matrix.txt.gz`
# 
# save(ManyInfections2, file = "./Data/ManyInfections2.rda")

load("./Data/ManyInfections2.rda")

Expr_ManyInfections2 <- exprs(ManyInfections2)
Pheno_ManyInfections2 <- pData(ManyInfections2)
FeatData_ManyInfections2 <- fData(ManyInfections2)

############################
## Annotation

## Expr_ManyInfections2
head(rownames(Expr_ManyInfections2))
rownames(Expr_ManyInfections2) <- FeatData_ManyInfections2$`Gene symbol`
summary(is.na(rownames(Expr_ManyInfections2)))
#rownames(Expr_ManyInfections2) <- gsub("-","", rownames(Expr_ManyInfections2))
#rownames(Expr_ManyInfections2) <- gsub("_","",rownames(Expr_ManyInfections2))
sel <- which(apply(Expr_ManyInfections2, 1, function(x) all(is.finite(x)) ))
Expr_ManyInfections2 <- Expr_ManyInfections2[sel, ]
Expr_ManyInfections2 <- Expr_ManyInfections2[!is.na(rownames(Expr_ManyInfections2)),]
dim(Expr_ManyInfections2)

range(Expr_ManyInfections2)

Expr_ManyInfections2 <- log2(Expr_ManyInfections2)
Expr_ManyInfections2 <- t(scale(t(Expr_ManyInfections2), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# Asymptomatic vs severe disease

# Pheno1

Pheno_ManyInfections2$DiseaseStatus <- as.factor(Pheno_ManyInfections2$`Pathogen:ch1`)
levels(Pheno_ManyInfections2$DiseaseStatus) <- c("case", "case", "control", "case", "case", "case") 
table(Pheno_ManyInfections2$DiseaseStatus)
Pheno_ManyInfections2$DiseaseStatus <- factor(Pheno_ManyInfections2$DiseaseStatus, levels = c("control", "case"))

all(rownames(Pheno_ManyInfections2) == colnames(Expr_ManyInfections2))

ClassInfectionsVsHealthy<- Pheno_ManyInfections2$DiseaseStatus

####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict

TestingData_ManyInfections <- t(Expr_ManyInfections2)

PredVotes_ManyInfections <- predict(RF_Comp, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections <- predict(RF_Comp, TestingData_ManyInfections, type="response")

ROCTest <- roc(ClassInfectionsVsHealthy, PredVotes_ManyInfections[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_ManyInfections <- evalmod(scores = PredVotes_ManyInfections[,2], labels = ClassInfectionsVsHealthy)
sscurves_ManyInfections
ROC_ManyInfections2 <- autoplot(sscurves_ManyInfections, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE6269-GPL96 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.35"), size = 4)
PRC_ManyInfections2 <- autoplot(sscurves_ManyInfections, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE6269-GPL96 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.92"), size = 4)

save(ROC_ManyInfections2, PRC_ManyInfections2, file = "./Objs/ManyInfections2_Curves.rda")

########################################################################################


####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

#################
## Predict 

PredVotes_ManyInfections2_cerebral <- predict(RF_Cerebral, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections2_cerebral <- predict(RF_Cerebral, TestingData_ManyInfections, type="response")

# For ROC and PRC curves
sscurves_ManyInfections2_cerebral <- evalmod(scores = PredVotes_ManyInfections2_cerebral[,2], labels = ClassInfectionsVsHealthy)
sscurves_ManyInfections2_cerebral
ROC_ManyInfections2_cerebral <- autoplot(sscurves_ManyInfections2_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE6269-GPL96 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.26"), size = 4)
PRC_ManyInfections2_cerebral <- autoplot(sscurves_ManyInfections2_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE6269-GPL96 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.90"), size = 4)

save(ROC_ManyInfections2_cerebral, PRC_ManyInfections2_cerebral, file = "./Objs/ManyInfections2_Curves_cerebral.rda")

############################
sessionInfo()