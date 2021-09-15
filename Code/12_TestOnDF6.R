##################################
## Test both the severe and cerebral malaria signatures on Dengue fever

rm(list = ls())

setwd("~/Documents/Research/Projects/Malaria")

renv::activate("~/Documents/Research/Projects/Malaria")

library(GEOquery)
library(randomForest)
library(pROC)
library(precrec)
library(ggplot2)

# DengueDataset6 <- getGEO("GSE13052", GSEMatrix = T, AnnotGPL = T)
# DengueDataset6 <- DengueDataset6$GSE13052_series_matrix.txt.gz
# 
# save(DengueDataset6, file = "./Data/DengueDataset6.rda")

load("./Data/DengueDataset6.rda")


Expr_Dengue6 <- exprs(DengueDataset6)
Pheno_Dengue6 <- pData(DengueDataset6)
FeatData_Dengue6 <- fData(DengueDataset6)


############################
## Annotation

## Expr_Dengue6
head(rownames(Expr_Dengue6))
rownames(Expr_Dengue6) <- FeatData_Dengue6$`Gene symbol`
summary(is.na(rownames(Expr_Dengue6)))
#rownames(Expr_Dengue6) <- gsub("-","", rownames(Expr_Dengue6))
#rownames(Expr_Dengue6) <- gsub("_","",rownames(Expr_Dengue6))
sel <- which(apply(Expr_Dengue6, 1, function(x) all(is.finite(x)) ))
Expr_Dengue6 <- Expr_Dengue6[sel, ]
Expr_Dengue6 <- Expr_Dengue6[!is.na(rownames(Expr_Dengue6)),]
dim(Expr_Dengue6)

range(Expr_Dengue6) 
Expr_Dengue6 <- log2(Expr_Dengue6 + 68)
Expr_Dengue6 <- t(scale(t(Expr_Dengue6), center = TRUE, scale = TRUE))


####################################


########################################################################################
#########################################################################################
###### DF vs DHF and DSS

### Modify the phenotype

# Pheno1
Pheno_Dengue6$DiseaseStatus <- as.factor(Pheno_Dengue6$characteristics_ch1)
levels(Pheno_Dengue6$DiseaseStatus) <- c("DSS", "DSS", "DF", "DF") 
table(Pheno_Dengue6$DiseaseStatus)
Pheno_Dengue6$DiseaseStatus <- factor(Pheno_Dengue6$DiseaseStatus, levels = c("DF", "DSS"))

all(rownames(Pheno_Dengue6) == colnames(Expr_Dengue6))

ClassDFvsDSS <- Pheno_Dengue6$DiseaseStatus

####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (DF vs DSS)

TestingData_Dengue <- t(Expr_Dengue6)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDFvsDSS, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DSS"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue6 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDFvsDSS)
sscurves_Dengue6
ROC_Dengue6 <- autoplot(sscurves_Dengue6, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE13052 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.64"), size = 3)
PRC_Dengue6 <- autoplot(sscurves_Dengue6, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE13052 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.566"), size = 3)

save(ROC_Dengue6, PRC_Dengue6, file = "./Objs/Dengue6_Curves.rda")


####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the Dengue dataset (DF vs DSS)
PredVotes_Dengue_cerebral <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue_cerebral <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest_cerebral <- roc(ClassDFvsDSS, PredVotes_Dengue_cerebral[,2], plot = F, print.auc = TRUE, levels = c("DF", "DSS"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest_cerebral

# For ROC and PRC curves
sscurves_Dengue6_cerebral <- evalmod(scores = PredVotes_Dengue_cerebral[,2], labels = ClassDFvsDSS)
sscurves_Dengue6_cerebral
ROC_Dengue6_cerebral <- autoplot(sscurves_Dengue6_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE13052 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.30"), size = 3)
PRC_Dengue6_cerebral <- autoplot(sscurves_Dengue6_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE13052 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.38"), size = 3)

save(ROC_Dengue6_cerebral, PRC_Dengue6_cerebral, file = "./Objs/Dengue6_Curves_cerebral.rda")


