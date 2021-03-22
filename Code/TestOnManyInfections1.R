
rm(list = ls())

library(GEOquery)

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
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the WestNile dataset (WestNile vs normal)

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
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the WestNile dataset (DHF vs DF)

#TestingData_ManyInfections <- t(Expr_ManyInfections1)

PredVotes_ManyInfections <- predict(RF_Cerebral, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections <- predict(RF_Cerebral, TestingData_ManyInfections, type="response")

ROCTest <- roc(ClassInfectionsVsHealthy, PredVotes_ManyInfections[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
