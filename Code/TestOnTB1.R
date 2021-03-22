rm(list = ls())

library(GEOquery)

# TB1 <- getGEO("GSE19444", GSEMatrix = T, AnnotGPL = T)
# TB1 <- TB1$GSE19444_series_matrix.txt.gz
# 
# save(TB1, file = "./Data/TB1.rda")
# 
load("./Data/TB1.rda")

Expr_TB1 <- exprs(TB1)
Pheno_TB1 <- pData(TB1)
FeatData_TB1 <- fData(TB1)

############################
## Annotation

## Expr_TB1
head(rownames(Expr_TB1))
rownames(Expr_TB1) <- FeatData_TB1$`Gene symbol`
summary(is.na(rownames(Expr_TB1)))
#rownames(Expr_TB1) <- gsub("-","", rownames(Expr_TB1))
#rownames(Expr_TB1) <- gsub("_","",rownames(Expr_TB1))
sel <- which(apply(Expr_TB1, 1, function(x) all(is.finite(x)) ))
Expr_TB1 <- Expr_TB1[sel, ]
Expr_TB1 <- Expr_TB1[!is.na(rownames(Expr_TB1)),]
dim(Expr_TB1)

range(Expr_TB1)
Expr_TB1 <- log2(Expr_TB1 + 25)
Expr_TB1 <- t(scale(t(Expr_TB1), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# primary TB vs latent TB and healthy

# Pheno1

Pheno_TB1$DiseaseStatus <- as.factor(Pheno_TB1$`illness:ch1`)
levels(Pheno_TB1$DiseaseStatus) <- c("control", "control", "case")
Pheno_TB1$DiseaseStatus <- factor(Pheno_TB1$DiseaseStatus, levels = c("control", "case"))
table(Pheno_TB1$DiseaseStatus)
all(rownames(Pheno_TB1) == colnames(Expr_TB1))

ClassTBVsHealthy<- Pheno_TB1$DiseaseStatus

####################################
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the WestNile dataset (WestNile vs normal)

TestingData_TB1 <- t(Expr_TB1)

PredVotes_TB1 <- predict(RF_Comp, newdata = TestingData_TB1, type = "vote")
PredResponse_TB1 <- predict(RF_Comp, TestingData_TB1, type="response")

ROCTest <- roc(ClassTBVsHealthy, PredVotes_TB1[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_TB1 <- evalmod(scores = PredVotes_TB1[,2], labels = ClassTBVsHealthy)
sscurves_TB1
ROC_TB1 <- autoplot(sscurves_TB1, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE19444 (primary TB vs latent TB and healthy)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.32"), size = 4)
PRC_TB1 <- autoplot(sscurves_TB1, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE19444 (primary TB vs latent TB and healthy)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.286"), size = 4)

save(ROC_TB1, PRC_TB1, file = "./Objs/TB1_Curves.rda")

########################################################################################


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the WestNile dataset (DHF vs DF)

#TestingData_ManyInfections <- t(Expr_TB1)

PredVotes_TB1 <- predict(RF_Cerebral, newdata = TestingData_TB1, type = "vote")
PredResponse_TB1 <- predict(RF_Cerebral, TestingData_TB1, type="response")

ROCTest <- roc(ClassTBVsHealthy, PredVotes_TB1[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
