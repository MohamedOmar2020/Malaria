rm(list = ls())

library(GEOquery)

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
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the WestNile dataset (WestNile vs normal)

TestingData_ManyInfections <- t(Expr_ManyInfections2)

PredVotes_ManyInfections <- predict(RF_Comp, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections <- predict(RF_Comp, TestingData_ManyInfections, type="response")

ROCTest <- roc(ClassInfectionsVsHealthy, PredVotes_ManyInfections[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_ManyInfections <- evalmod(scores = PredVotes_ManyInfections[,2], labels = ClassInfectionsVsHealthy)
sscurves_ManyInfections
ROC_ManyInfections2 <- autoplot(sscurves_ManyInfections, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE6269-GPL96 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.35"), size = 5)
PRC_ManyInfections2 <- autoplot(sscurves_ManyInfections, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE6269-GPL96 (Multiple infections vs control)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.92"), size = 5)

save(ROC_ManyInfections2, PRC_ManyInfections2, file = "./Objs/ManyInfections2_Curves.rda")

########################################################################################


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the WestNile dataset (DHF vs DF)

#TestingData_ManyInfections <- t(Expr_ManyInfections2)

PredVotes_ManyInfections <- predict(RF_Cerebral, newdata = TestingData_ManyInfections, type = "vote")
PredResponse_ManyInfections <- predict(RF_Cerebral, TestingData_ManyInfections, type="response")

ROCTest <- roc(ClassInfectionsVsHealthy, PredVotes_ManyInfections[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
