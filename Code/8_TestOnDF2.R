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


#DengueDataset2 <- getGEO("GSE96656", GSEMatrix = T, AnnotGPL = T)
#DengueDataset2 <- DengueDataset2$GSE96656_series_matrix.txt.gz

#save(DengueDataset2, file = "./Data/DengueDataset2.rda")

load("./Data/DengueDataset2.rda")

Expr_Dengue2 <- exprs(DengueDataset2)
Pheno_Dengue2 <- pData(DengueDataset2)
FeatData_Dengue2 <- fData(DengueDataset2)


############################
## Annotation

## Expr_Dengue2
head(rownames(Expr_Dengue2))
rownames(Expr_Dengue2) <- FeatData_Dengue2$ORF
summary(is.na(rownames(Expr_Dengue2)))
#rownames(Expr_Dengue2) <- gsub("-","", rownames(Expr_Dengue2))
#rownames(Expr_Dengue2) <- gsub("_","",rownames(Expr_Dengue2))
sel <- which(apply(Expr_Dengue2, 1, function(x) all(is.finite(x)) ))
Expr_Dengue2 <- Expr_Dengue2[sel, ]
Expr_Dengue2 <- Expr_Dengue2[!is.na(rownames(Expr_Dengue2)),]
dim(Expr_Dengue2)

range(Expr_Dengue2)  

# Already Z-transformed

#Expr_Dengue2 <- t(scale(t(Expr_Dengue2), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# Denge Fever and DHF vs Healthy

# Pheno1
Pheno_Dengue2$DiseaseStatus <- as.factor(Pheno_Dengue2$`disease state:ch2`)
levels(Pheno_Dengue2$DiseaseStatus) <- c("case", "case", "control") 
table(Pheno_Dengue2$DiseaseStatus)
Pheno_Dengue2$DiseaseStatus <- factor(Pheno_Dengue2$DiseaseStatus, levels = c("control", "case"))

ClassDengueVsNormal <- Pheno_Dengue2$DiseaseStatus

####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

## Some features (2) are present in the RF model but not in the expression matrix >> removed them
CommonGns <- intersect(rownames(Expr_Dengue2), rownames(RF_Comp$importance))
RF_Comp$importance <- RF_Comp$importance[CommonGns, ]
#RF_Comp$importanceSD <- RF_Comp$importanceSD[CommonGns, ]
RF_Comp$forest$ncat <- RF_Comp$forest$ncat[CommonGns]

####
## Predict in the Dengue dataset (Dengue vs normal)

TestingData_Dengue <- t(Expr_Dengue2)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROC_DF_severe <- roc(ClassDengueVsNormal, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROC_DF_severe

# For ROC and PRC curves
sscurves_Dengue2 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDengueVsNormal)
sscurves_Dengue2
ROC_Dengue2 <- autoplot(sscurves_Dengue2, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE96656 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = ", round(ROC_DF_severe$auc, 2)), size = 3)
PRC_Dengue2 <- autoplot(sscurves_Dengue2, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE96656 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.67"), size = 3)

save(ROC_Dengue2, PRC_Dengue2, file = "./Objs/Dengue2_Curves.rda")

####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

## Some features (2) are present in the RF model but not in the expression matrix >> removed them
CommonGns <- intersect(rownames(Expr_Dengue2), rownames(RF_Cerebral$importance))
RF_Cerebral$importance <- RF_Cerebral$importance[CommonGns, ]
#RF_Cerebral$importanceSD <- RF_Cerebral$importanceSD[CommonGns, ]
RF_Cerebral$forest$ncat <- RF_Cerebral$forest$ncat[CommonGns]

####
## Predict in the Dengue dataset (Dengue vs normal)

PredVotes_Dengue_cerebral <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue_cerebral <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDengueVsNormal, PredVotes_Dengue_cerebral[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue2_cerebral <- evalmod(scores = PredVotes_Dengue_cerebral[,2], labels = ClassDengueVsNormal)
sscurves_Dengue2_cerebral
ROC_Dengue2_cerebral <- autoplot(sscurves_Dengue2_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE96656 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.92"), size = 3)
PRC_Dengue2_cerebral <- autoplot(sscurves_Dengue2_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE96656 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.97"), size = 3)

save(ROC_Dengue2_cerebral, PRC_Dengue2_cerebral, file = "./Objs/Dengue2_Curves_cerebral.rda")

########################################################################################
#########################################################################################
###### HDF vs DF

### Modify the phenotype
# Remove controls

# Pheno1
Pheno_Dengue2 <- Pheno_Dengue2[!(Pheno_Dengue2$`disease state:ch2` == "Healthy"), ]
Pheno_Dengue2$DiseaseStatus2 <- as.factor(Pheno_Dengue2$`disease state:ch2`)
levels(Pheno_Dengue2$DiseaseStatus2) <- c("DF", "DHF") 
table(Pheno_Dengue2$DiseaseStatus2)

Expr_Dengue2 <- Expr_Dengue2[, colnames(Expr_Dengue2) %in% rownames(Pheno_Dengue2)]
all(rownames(Pheno_Dengue2) == colnames(Expr_Dengue2))

ClassDHFvsDF <- Pheno_Dengue2$DiseaseStatus2

####################################
## Load the complicated malaria signature
#load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (DHF vs DF)

TestingData_Dengue <- t(Expr_Dengue2)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDHFvsDF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue2 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDHFvsDF)
sscurves_Dengue2
ROC_Dengue2_DFvsDHF <- autoplot(sscurves_Dengue2, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE96656 (DF vs DHF)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.37"), size = 5)
PRC_Dengue2_DFvsDHF <- autoplot(sscurves_Dengue2, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE96656 (DF vs DHF)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.71"), size = 5)

save(ROC_Dengue2_DFvsDHF, PRC_Dengue2_DFvsDHF, file = "./Objs/Dengue2_Curves_DFvsDHF.rda")

####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

## Some features (2) are present in the RF model but not in the expression matrix >> removed them
CommonGns <- intersect(rownames(Expr_Dengue2), rownames(RF_Cerebral$importance))
RF_Cerebral$importance <- RF_Cerebral$importance[CommonGns, ]

#################
## Predict in the Dengue dataset (DHF vs DF)

PredVotes_Dengue2_cerebral <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue2_cerebral <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROC_DHF_cerebral <- roc(ClassDHFvsDF, PredVotes_Dengue2_cerebral[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROC_DHF_cerebral

# For ROC and PRC curves
sscurves_DHF2_cerebral <- evalmod(scores = PredVotes_Dengue2_cerebral[,2], labels = ClassDHFvsDF)
sscurves_DHF2_cerebral
ROC_Dengue2_DFvsDHF_cerebral <- autoplot(sscurves_DHF2_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE96656 (DF vs DHF)") + annotate("text", x = .65, y = .25, label = paste("AUC = ", round(ROC_DHF_cerebral$auc, 2)), size = 5)
ROC_Dengue2_DFvsDHF_cerebral <- autoplot(sscurves_DHF2_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE96656 (DF vs DHF)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.23"), size = 5)

save(ROC_Dengue2_DFvsDHF_cerebral, ROC_Dengue2_DFvsDHF_cerebral, file = "./Objs/Dengue2_Curves_DFvsDHF_cerebral.rda")

###############
sessionInfo()