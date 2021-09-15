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

#DengueDataset1 <- getGEO("GSE51808", GSEMatrix = T, AnnotGPL = T)
#DengueDataset1 <- DengueDataset1$GSE51808_series_matrix.txt.gz

#save(DengueDataset1, file = "./Data/DengueDataset1.rda")

load("./Data/DengueDataset1.rda")

Expr_Dengue1 <- exprs(DengueDataset1)
Pheno_Dengue1 <- pData(DengueDataset1)
FeatData_Dengue1 <- fData(DengueDataset1)


############################
## Annotation

## Expr_Dengue1
head(rownames(Expr_Dengue1))
rownames(Expr_Dengue1) <- FeatData_Dengue1$`Gene symbol`
summary(is.na(rownames(Expr_Dengue1)))
sel <- which(apply(Expr_Dengue1, 1, function(x) all(is.finite(x)) ))
Expr_Dengue1 <- Expr_Dengue1[sel, ]
Expr_Dengue1 <- Expr_Dengue1[!is.na(rownames(Expr_Dengue1)),]
dim(Expr_Dengue1)

range(Expr_Dengue1)
 
Expr_Dengue1 <- t(scale(t(Expr_Dengue1), center = TRUE, scale = TRUE))


####################################
### Modify the phenotype

# Control and convalescent VS DF and DHF
# Pheno1
Pheno_Dengue1$DiseaseStatus <- as.factor(Pheno_Dengue1$`status:ch1`)
levels(Pheno_Dengue1$DiseaseStatus) <- c("control", "control", "case", "case") 
table(Pheno_Dengue1$DiseaseStatus)


ClassDengueVsNormal <- Pheno_Dengue1$DiseaseStatus

####################################
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (Dengue vs normal)

TestingData_Dengue <- t(Expr_Dengue1)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDengueVsNormal, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue1 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDengueVsNormal)
sscurves_Dengue1
ROC_Dengue <- autoplot(sscurves_Dengue1, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE51808 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.52"), size = 3)
PRC_Dengue <- autoplot(sscurves_Dengue1, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE51808 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.54"), size = 3)

save(ROC_Dengue, PRC_Dengue, file = "./Objs/Dengue1_Curves.rda")

####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the Dengue dataset (DF vs normal)

PredVotes_Dengue <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDengueVsNormal, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue1_cerebral <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDengueVsNormal)
sscurves_Dengue1_cerebral
ROC_Dengue_cerebral <- autoplot(sscurves_Dengue1_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE51808 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.50"), size = 3)
PRC_Dengue_cerebral <- autoplot(sscurves_Dengue1_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE51808 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.50"), size = 3)

save(ROC_Dengue_cerebral, PRC_Dengue_cerebral, file = "./Objs/Dengue1_Curves_cerebral.rda")

########################################################################################
#########################################################################################
###### HDF vs DF

### Modify the phenotype
# Remove controls

# Pheno1
Pheno_Dengue1 <- Pheno_Dengue1[!(Pheno_Dengue1$`status:ch1` %in% c("convalescent", "control")), ]
Pheno_Dengue1$DiseaseStatus2 <- as.factor(Pheno_Dengue1$`status:ch1`)
#levels(Pheno_Dengue1$DiseaseStatus) <- c("DF", "DHF") 
table(Pheno_Dengue1$DiseaseStatus2)

Expr_Dengue1 <- Expr_Dengue1[, colnames(Expr_Dengue1) %in% rownames(Pheno_Dengue1)]
all(rownames(Pheno_Dengue1) == colnames(Expr_Dengue1))

ClassDHFvsDF <- Pheno_Dengue1$DiseaseStatus2

###############
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

## Predict in the Dengue dataset (DHF vs DF)
TestingData_Dengue <- t(Expr_Dengue1)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDHFvsDF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue1 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDHFvsDF)
sscurves_Dengue1
ROC_Dengue_DFvsDHF <- autoplot(sscurves_Dengue1, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE51808 (DF vs DHF)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.57"), size = 5)
PRC_Dengue_DFvsDHF <- autoplot(sscurves_Dengue1, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE51808 (DF vs DHF)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.39"), size = 5)

save(ROC_Dengue_DFvsDHF, PRC_Dengue_DFvsDHF, file = "./Objs/Dengue1_Curves_DFvsDHF.rda")


###############
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

## Predict in the Dengue dataset (DHF vs DF)

PredVotes_Dengue_cerebral <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue_cerebral <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDHFvsDF, PredVotes_Dengue_cerebral[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue1_cerebral <- evalmod(scores = PredVotes_Dengue_cerebral[,2], labels = ClassDHFvsDF)
sscurves_Dengue1_cerebral
ROC_Dengue_DFvsDHF_cerebral <- autoplot(sscurves_Dengue1_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE51808 (DF vs DHF)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.51"), size = 5)
PRC_Dengue_DFvsDHF_cerebral <- autoplot(sscurves_Dengue1_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE51808 (DF vs DHF)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.35"), size = 5)

save(ROC_Dengue_DFvsDHF_cerebral, PRC_Dengue_DFvsDHF_cerebral, file = "./Objs/Dengue1_Curves_DFvsDHF_cerebral.rda")

##########################
sessionInfo()