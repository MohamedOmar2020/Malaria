
##################################
## Test both the severe and cerebral malaria signatures on West Nile Virus infection

rm(list = ls())

setwd("~/Documents/Research/Projects/Malaria")

renv::activate("~/Documents/Research/Projects/Malaria")

library(GEOquery)
library(randomForest)
library(pROC)
library(precrec)
library(ggplot2)

# WestNileDataset1 <- getGEO("GSE46681", GSEMatrix = T, AnnotGPL = T)
# WestNileDataset1 <- WestNileDataset1$GSE46681_series_matrix.txt.gz
# 
# save(WestNileDataset1, file = "./Data/WestNileDataset1.rda")

load("./Data/WestNileDataset1.rda")

Expr_WestNile1 <- exprs(WestNileDataset1)
Pheno_WestNile1 <- pData(WestNileDataset1)
FeatData_WestNile1 <- fData(WestNileDataset1)


############################
## Annotation

## Expr_WestNile1
head(rownames(Expr_WestNile1))
rownames(Expr_WestNile1) <- FeatData_WestNile1$`Gene symbol`
summary(is.na(rownames(Expr_WestNile1)))
#rownames(Expr_WestNile1) <- gsub("-","", rownames(Expr_WestNile1))
#rownames(Expr_WestNile1) <- gsub("_","",rownames(Expr_WestNile1))
sel <- which(apply(Expr_WestNile1, 1, function(x) all(is.finite(x)) ))
Expr_WestNile1 <- Expr_WestNile1[sel, ]
Expr_WestNile1 <- Expr_WestNile1[!is.na(rownames(Expr_WestNile1)),]
dim(Expr_WestNile1)

range(Expr_WestNile1)
#plot(density(Expr_WestNile1))
#boxplot(Expr_WestNile1)
# X1 <- Expr_WestNile1
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt1 <- genefilter(2^X1,ffun)
# Expr_WestNile1 <- Expr_WestNile1[filt1,]
# 
Expr_WestNile1 <- t(scale(t(Expr_WestNile1), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# Asymptomatic vs severe disease

# Pheno1
Pheno_WestNile1 <- Pheno_WestNile1[Pheno_WestNile1$`cell type:ch1` == "peripheral blood mononuclear cells (PBMCs)", ]
Pheno_WestNile1 <- Pheno_WestNile1[!is.na(Pheno_WestNile1$`disease status:ch1`), ]

Pheno_WestNile1$DiseaseStatus <- as.factor(Pheno_WestNile1$`disease status:ch1`)
levels(Pheno_WestNile1$DiseaseStatus) <- c("control", "case") 
table(Pheno_WestNile1$DiseaseStatus)

Expr_WestNile1 <- Expr_WestNile1[, colnames(Expr_WestNile1) %in% rownames(Pheno_WestNile1)]
all(rownames(Pheno_WestNile1) == colnames(Expr_WestNile1))

ClassSeverVsAsymp <- Pheno_WestNile1$DiseaseStatus

####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict in the WestNile dataset (WestNile vs normal)

TestingData_WestNile <- t(Expr_WestNile1)

PredVotes_WestNile <- predict(RF_Comp, newdata = TestingData_WestNile, type = "vote")
PredResponse_WestNile <- predict(RF_Comp, TestingData_WestNile, type="response")

ROCTest <- roc(ClassSeverVsAsymp, PredVotes_WestNile[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_WestNile1 <- evalmod(scores = PredVotes_WestNile[,2], labels = ClassSeverVsAsymp)
sscurves_WestNile1
ROC_WestNile <- autoplot(sscurves_WestNile1, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE46681 (West Nile)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.55"), size = 4)
PRC_WestNile <- autoplot(sscurves_WestNile1, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE46681 (West Nile)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.466"), size = 4)

save(ROC_WestNile, PRC_WestNile, file = "./Objs/WestNile1_Curves.rda")

########################################################################################


####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the WestNile dataset (DHF vs DF)
PredVotes_WestNile_cerebral <- predict(RF_Cerebral, newdata = TestingData_WestNile, type = "vote")
PredResponse_WestNile_cerebral <- predict(RF_Cerebral, TestingData_WestNile, type="response")

ROCTest_cerebral <- roc(ClassSeverVsAsymp, PredVotes_WestNile_cerebral[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest_cerebral

# For ROC and PRC curves
sscurves_WestNile1_cerebral <- evalmod(scores = PredVotes_WestNile_cerebral[,2], labels = ClassSeverVsAsymp)
sscurves_WestNile1_cerebral
ROC_WestNile_cerebral <- autoplot(sscurves_WestNile1_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE46681 (West Nile)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.40"), size = 4)
PRC_WestNile_cerebral <- autoplot(sscurves_WestNile1_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE46681 (West Nile)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.44"), size = 4)

save(ROC_WestNile_cerebral, PRC_WestNile_cerebral, file = "./Objs/WestNile1_Curves_cerebral.rda")

#####################
sessionInfo()
