rm(list = ls())

library(GEOquery)

# TB2 <- getGEO("GSE73408", GSEMatrix = T, AnnotGPL = T)
# TB2 <- TB2$GSE73408_series_matrix.txt.gz
# 
# save(TB2, file = "./Data/TB2.rda")

load("./Data/TB2.rda")

Expr_TB2 <- exprs(TB2)
Pheno_TB2 <- pData(TB2)
FeatData_TB2 <- fData(TB2)

############################
## Annotation

## Expr_TB2
head(rownames(Expr_TB2))
rownames(Expr_TB2) <- FeatData_TB2$`Gene symbol`
summary(is.na(rownames(Expr_TB2)))
#rownames(Expr_TB2) <- gsub("-","", rownames(Expr_TB2))
#rownames(Expr_TB2) <- gsub("_","",rownames(Expr_TB2))
sel <- which(apply(Expr_TB2, 1, function(x) all(is.finite(x)) ))
Expr_TB2 <- Expr_TB2[sel, ]
Expr_TB2 <- Expr_TB2[!is.na(rownames(Expr_TB2)),]
dim(Expr_TB2)

range(Expr_TB2)
Expr_TB2 <- t(scale(t(Expr_TB2), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# primary TB and pneumonia vs latent TB

# Pheno1

Pheno_TB2$DiseaseStatus <- as.factor(Pheno_TB2$`clinical group:ch1`)
levels(Pheno_TB2$DiseaseStatus) <- c("control", "case", "case")
#Pheno_TB2$DiseaseStatus <- factor(Pheno_TB2$DiseaseStatus, levels = c("control", "case"))
table(Pheno_TB2$DiseaseStatus)
all(rownames(Pheno_TB2) == colnames(Expr_TB2))

ClassTBandPNAVsLatentTB<- Pheno_TB2$DiseaseStatus

####################################
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the WestNile dataset (WestNile vs normal)

TestingData_TB2 <- t(Expr_TB2)

PredVotes_TB2 <- predict(RF_Comp, newdata = TestingData_TB2, type = "vote")
PredResponse_TB2 <- predict(RF_Comp, TestingData_TB2, type="response")

ROCTest <- roc(ClassTBandPNAVsLatentTB, PredVotes_TB2[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_TB2 <- evalmod(scores = PredVotes_TB2[,2], labels = ClassTBandPNAVsLatentTB)
sscurves_TB2
ROC_TB2 <- autoplot(sscurves_TB2, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE73408 (primary TB and pneumonia vs latent TB)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.55"), size = 4)
PRC_TB2 <- autoplot(sscurves_TB2, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE73408 (primary TB and pneumonia vs latent TB)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.73"), size = 4)

save(ROC_TB2, PRC_TB2, file = "./Objs/TB2_Curves.rda")

########################################################################################


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the WestNile dataset (DHF vs DF)
RF_Cerebral$importance <- RF_Cerebral$importance[rownames(RF_Cerebral$importance) %in% rownames(Expr_TB2), ]
#TestingData_ManyInfections <- t(Expr_TB2)

PredVotes_TB2 <- predict(RF_Cerebral, newdata = TestingData_TB2, type = "vote")
PredResponse_TB2 <- predict(RF_Cerebral, TestingData_TB2, type="response")

ROCTest <- roc(ClassTBandPNAVsLatentTB, PredVotes_TB2[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
