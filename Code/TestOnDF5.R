rm(list = ls())

library(GEOquery)

# DengueDataset5 <- getGEO("GSE17924", GSEMatrix = T, AnnotGPL = T)
# DengueDataset5 <- DengueDataset5$GSE17924_series_matrix.txt.gz
# 
# save(DengueDataset5, file = "./Data/DengueDataset5.rda")

load("./Data/DengueDataset5.rda")


Expr_Dengue5 <- exprs(DengueDataset5)
Pheno_Dengue5 <- pData(DengueDataset5)
FeatData_Dengue5 <- fData(DengueDataset5)


############################
## Annotation

## Expr_Dengue5
head(rownames(Expr_Dengue5))
rownames(Expr_Dengue5) <- FeatData_Dengue5$`Gene symbol`
summary(is.na(rownames(Expr_Dengue5)))
#rownames(Expr_Dengue5) <- gsub("-","", rownames(Expr_Dengue5))
#rownames(Expr_Dengue5) <- gsub("_","",rownames(Expr_Dengue5))
sel <- which(apply(Expr_Dengue5, 1, function(x) all(is.finite(x)) ))
Expr_Dengue5 <- Expr_Dengue5[sel, ]
Expr_Dengue5 <- Expr_Dengue5[!is.na(rownames(Expr_Dengue5)),]
dim(Expr_Dengue5)

range(Expr_Dengue5) 

Expr_Dengue5 <- t(scale(t(Expr_Dengue5), center = TRUE, scale = TRUE))


####################################


########################################################################################
#########################################################################################
###### DF vs DHF and DSS

### Modify the phenotype

# Pheno1
Pheno_Dengue5$DiseaseStatus <- as.factor(Pheno_Dengue5$`severity:ch1`)
levels(Pheno_Dengue5$DiseaseStatus) <- c("DF", "ComplicatedDF", "ComplicatedDF") 
table(Pheno_Dengue5$DiseaseStatus)

all(rownames(Pheno_Dengue5) == colnames(Expr_Dengue5))

ClassDFvsComDengue <- Pheno_Dengue5$DiseaseStatus

####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict

TestingData_Dengue <- t(Expr_Dengue5)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDFvsComDengue, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "ComplicatedDF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue5 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDFvsComDengue)
sscurves_Dengue5
ROC_Dengue5 <- autoplot(sscurves_Dengue5, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE17924 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.52"), size = 3)
PRC_Dengue5 <- autoplot(sscurves_Dengue5, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE17924 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.755"), size = 3)

save(ROC_Dengue5, PRC_Dengue5, file = "./Objs/Dengue5_Curves.rda")


####################################
## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

#################
## Predict
PredVotes_Dengue_cerebral <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue_cerebral <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest_cerebral <- roc(ClassDFvsComDengue, PredVotes_Dengue_cerebral[,2], plot = F, print.auc = TRUE, levels = c("DF", "ComplicatedDF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest_cerebral

# For ROC and PRC curves
sscurves_Dengue5_cerebral <- evalmod(scores = PredVotes_Dengue_cerebral[,2], labels = ClassDFvsComDengue)
sscurves_Dengue5_cerebral
ROC_Dengue5_cerebral <- autoplot(sscurves_Dengue5_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE17924 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.66"), size = 3)
PRC_Dengue5_cerebral <- autoplot(sscurves_Dengue5_cerebral, curvetype = c("PRC")) + labs(title = "PRC curve of the cerebral malaria signature in GSE17924 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.78"), size = 3)

save(ROC_Dengue5_cerebral, PRC_Dengue5_cerebral, file = "./Objs/Dengue5_Curves_cerebral.rda")

