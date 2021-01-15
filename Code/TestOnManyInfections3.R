rm(list = ls())

library(GEOquery)

# ManyInfections3 <- getGEO("GSE63990", GSEMatrix = T, AnnotGPL = T)
# ManyInfections3 <- ManyInfections3$GSE63990_series_matrix.txt.gz
# 
# save(ManyInfections3, file = "./Data/ManyInfections3.rda")

load("./Data/ManyInfections3.rda")

Expr_ManyInfections3 <- exprs(ManyInfections3)
Pheno_ManyInfections3 <- pData(ManyInfections3)
FeatData_ManyInfections3 <- fData(ManyInfections3)

############################
## Annotation

## Expr_ManyInfections3
head(rownames(Expr_ManyInfections3))
rownames(Expr_ManyInfections3) <- FeatData_ManyInfections3$`Gene symbol`
summary(is.na(rownames(Expr_ManyInfections3)))
#rownames(Expr_ManyInfections3) <- gsub("-","", rownames(Expr_ManyInfections3))
#rownames(Expr_ManyInfections3) <- gsub("_","",rownames(Expr_ManyInfections3))
sel <- which(apply(Expr_ManyInfections3, 1, function(x) all(is.finite(x)) ))
Expr_ManyInfections3 <- Expr_ManyInfections3[sel, ]
Expr_ManyInfections3 <- Expr_ManyInfections3[!is.na(rownames(Expr_ManyInfections3)),]
dim(Expr_ManyInfections3)

range(Expr_ManyInfections3)

Expr_ManyInfections3 <- t(scale(t(Expr_ManyInfections3), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# Asymptomatic vs severe disease

# Pheno1

Pheno_ManyInfections3$DiseaseStatus <- as.factor(Pheno_ManyInfections3$`infection_status:ch1`)
levels(Pheno_ManyInfections3$DiseaseStatus) <- c("case", "control", "case") 
table(Pheno_ManyInfections3$DiseaseStatus)
Pheno_ManyInfections3$DiseaseStatus <- factor(Pheno_ManyInfections3$DiseaseStatus, levels = c("control", "case"))

all(rownames(Pheno_ManyInfections3) == colnames(Expr_ManyInfections3))

ClassInfectionsVsHealthy<- Pheno_ManyInfections3$DiseaseStatus

####################################
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the WestNile dataset (WestNile vs normal)

TestingData_ManyInfections <- t(Expr_ManyInfections3)

PredVotes_ManyInfections <- predict(RF_Comp, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections <- predict(RF_Comp, TestingData_ManyInfections, type="response")

ROCTest <- roc(ClassInfectionsVsHealthy, PredVotes_ManyInfections[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_ManyInfections <- evalmod(scores = PredVotes_ManyInfections[,2], labels = ClassInfectionsVsHealthy)
sscurves_ManyInfections
ROC_ManyInfections3 <- autoplot(sscurves_ManyInfections, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE63990 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.39"), size = 5)
PRC_ManyInfections3 <- autoplot(sscurves_ManyInfections, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE63990 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.62"), size = 5)

save(ROC_ManyInfections3, PRC_ManyInfections3, file = "./Objs/ManyInfections3_Curves.rda")

########################################################################################


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the WestNile dataset (DHF vs DF)

#TestingData_ManyInfections <- t(Expr_ManyInfections3)

PredVotes_ManyInfections <- predict(RF_Cerebral, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections <- predict(RF_Cerebral, TestingData_ManyInfections, type="response")

ROCTest <- roc(ClassInfectionsVsHealthy, PredVotes_ManyInfections[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
