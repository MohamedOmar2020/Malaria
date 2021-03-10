

rm(list = ls())

library(GEOquery)
library(pd.hugene.2.1.st)
library(affycoretools)
library(pROC)
library(caret)
library(randomForest)



ExternalDataset <- getGEO("GSE35859", GSEMatrix = T, AnnotGPL = T)
ExternalDataset <- ExternalDataset$GSE35859_series_matrix.txt.gz
save(ExternalDataset, file = "./Objs/ExternalDataset.rda")  

load("./Objs/ExternalDataset.rda")


Expr_Test <- exprs(ExternalDataset)
Pheno_Test <- pData(ExternalDataset)
FeatData_Test <- fData(ExternalDataset)


############################
## Annotation

## Expr_Cerebral
head(rownames(Expr_Test))
rownames(Expr_Test) <- FeatData_Test$`Composite Element Database Entry[Gene Symbol]`
summary(is.na(rownames(Expr_Test)))
#rownames(Expr_Cerebral) <- gsub("-","", rownames(Expr_Cerebral))
#rownames(Expr_Cerebral) <- gsub("_","",rownames(Expr_Cerebral))
sel <- which(apply(Expr_Test, 1, function(x) all(is.finite(x)) ))
Expr_Test <- Expr_Test[sel, ]
Expr_Test <- Expr_Test[!is.na(rownames(Expr_Test)),]
dim(Expr_Test)

range(Expr_Test)
#Expr_Test <- t(scale(t(Expr_Test), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

Pheno_Test <- Pheno_Test[!(Pheno_Test$`disease state:ch2` == "Healthy"), ]
Pheno_Test$DiseaseStatus <- as.factor(Pheno_Test$`disease state:ch2`)
levels(Pheno_Test$DiseaseStatus) <- c("Complicated", "unComplicated") 
table(Pheno_Test$DiseaseStatus)
Pheno_Test$DiseaseStatus <- factor(Pheno_Test$DiseaseStatus, levels = c("unComplicated", "Complicated"))

Expr_Test <- Expr_Test[, colnames(Expr_Test) %in% rownames(Pheno_Test)]
all(rownames(Pheno_Test) == colnames(Expr_Test))


ClassComplicatedVsUnComplicated <- Pheno_Test$DiseaseStatus


load("./Objs/RF_Comp.rda")
####
Sel <- rownames(RF_Comp$importance)
Sel <- intersect(Sel, rownames(Expr_Test))  
RF_Comp$importance <- RF_Comp$importance[Sel, ]

  
usedTestMat_Filt <- Expr_Test[Sel, ]

TestingData_Filt <- t(usedTestMat_Filt)

#######################
PredVotes_Test <- predict(RF_Comp, newdata = TestingData_Filt, type = "vote")
PredResponse_Test <- predict(RF_Comp, TestingData_Filt, type="response")

ROCTest <- roc(ClassComplicatedVsUnComplicated, PredVotes_Test[,2], plot = F, print.auc = TRUE, levels = c("unComplicated", "Complicated"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

### Resubstitution peRF_Compormance in the Test set
ConfusionTest <- confusionMatrix(PredResponse_Test, ClassComplicatedVsUnComplicated, positive = "Complicated", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = PredResponse_Test, actuals = ClassComplicatedVsUnComplicated)
MCC_Test

# For ROC and PRC curves
sscurves_Test_Comp <- evalmod(scores = PredVotes_Test[,2], labels = ClassComplicatedVsUnComplicated)
sscurves_Test_Comp
ROC_Test_Comp <- autoplot(sscurves_Test_Comp, curvetype = c("ROC")) + labs(title = "ROC curve 1st testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.85"), size = 5)
PRC_Test_Comp <- autoplot(sscurves_Test_Comp, curvetype = c("PRC")) + labs(title = "PRC curve 1st testing dataset") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.95"), size = 5)





