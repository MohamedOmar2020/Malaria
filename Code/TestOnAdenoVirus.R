


rm(list = ls())

library(GEOquery)

AdenovirusDataset1 <- getGEO("GSE40396", GSEMatrix = T, AnnotGPL = T)
AdenovirusDataset1 <- AdenovirusDataset1$GSE40396_series_matrix.txt.gz

save(AdenovirusDataset1, file = "./Data/AdenovirusDataset1.rda")

load("./Data/AdenovirusDataset1.rda")

Expr_Adenovirus1 <- exprs(AdenovirusDataset1)
Pheno_Adenovirus1 <- pData(AdenovirusDataset1)
FeatData_Adenovirus1 <- fData(AdenovirusDataset1)


############################
## Annotation

## Expr_Adenovirus1
head(rownames(Expr_Adenovirus1))
rownames(Expr_Adenovirus1) <- FeatData_Adenovirus1$`Gene symbol`
summary(is.na(rownames(Expr_Adenovirus1)))
#rownames(Expr_Adenovirus1) <- gsub("-","", rownames(Expr_Adenovirus1))
#rownames(Expr_Adenovirus1) <- gsub("_","",rownames(Expr_Adenovirus1))
sel <- which(apply(Expr_Adenovirus1, 1, function(x) all(is.finite(x)) ))
Expr_Adenovirus1 <- Expr_Adenovirus1[sel, ]
Expr_Adenovirus1 <- Expr_Adenovirus1[!is.na(rownames(Expr_Adenovirus1)),]
dim(Expr_Adenovirus1)

range(Expr_Adenovirus1)
 
Expr_Adenovirus1 <- t(scale(t(Expr_Adenovirus1), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# Many infections vs normal
# Pheno1
Pheno_Adenovirus1$DiseaseStatus <- as.factor(Pheno_Adenovirus1$`pathogen:ch1`)
levels(Pheno_Adenovirus1$DiseaseStatus) <- c("case", "case", "case", "case", "case", "case", "case", "control", "case", "case") 
table(Pheno_Adenovirus1$DiseaseStatus)
Pheno_Adenovirus1$DiseaseStatus <- factor(Pheno_Adenovirus1$DiseaseStatus, levels = c("control", "case"))

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Adenovirus1)]
#all(rownames(Pheno_Adenovirus1) == colnames(expr1))

ClassInfectionVsNormal <- Pheno_Adenovirus1$DiseaseStatus

####################################
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the Adenovirus dataset (Adenovirus vs normal)

TestingData_Adenovirus <- t(Expr_Adenovirus1)

PredVotes_Adenovirus <- predict(RF_Comp, newdata = TestingData_Adenovirus, type = "vote")
PredResponse_Adenovirus <- predict(RF_Comp, TestingData_Adenovirus, type="response")

ROCTest <- roc(ClassInfectionVsNormal, PredVotes_Adenovirus[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Adenovirus1 <- evalmod(scores = PredVotes_Adenovirus[,2], labels = ClassInfectionVsNormal)
sscurves_Adenovirus1
ROC_Adenovirus <- autoplot(sscurves_Adenovirus1, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE40396 (Many infections)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.60"), size = 5)
PRC_Adenovirus <- autoplot(sscurves_Adenovirus1, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE40396 (Many infections)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.736"), size = 5)

save(ROC_Adenovirus, PRC_Adenovirus, file = "./Objs/Adenovirus1_Curves.rda")

########################################################################################
#########################################################################################
###### AdenoVirus vs normal

### Modify the phenotype
# Remove controls

# Pheno1
Pheno_Adenovirus1 <- Pheno_Adenovirus1[Pheno_Adenovirus1$`pathogen:ch1` %in%  c("Adenovirus", "None"), ]
Pheno_Adenovirus1$DiseaseStatus2 <- as.factor(Pheno_Adenovirus1$`pathogen:ch1`)
levels(Pheno_Adenovirus1$DiseaseStatus2) <- c("AdenoVirus", "control") 
table(Pheno_Adenovirus1$DiseaseStatus2)
Pheno_Adenovirus1$DiseaseStatus2 <- factor(Pheno_Adenovirus1$DiseaseStatus2, levels = c("control", "AdenoVirus"))

Expr_Adenovirus1 <- Expr_Adenovirus1[, colnames(Expr_Adenovirus1) %in% rownames(Pheno_Adenovirus1)]
all(rownames(Pheno_Adenovirus1) == colnames(Expr_Adenovirus1))

ClassAdenoVsNormal <- Pheno_Adenovirus1$DiseaseStatus2

####################################
## Load the complicated malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict in the Adenovirus dataset (DHF vs DF)

TestingData_Adenovirus <- t(Expr_Adenovirus1)

PredVotes_Adenovirus <- predict(RF_Comp, newdata = TestingData_Adenovirus, type = "vote")
PredResponse_Adenovirus <- predict(RF_Comp, TestingData_Adenovirus, type="response")

ROCTest <- roc(ClassAdenoVsNormal, PredVotes_Adenovirus[,2], plot = F, print.auc = TRUE, levels = c("control", "AdenoVirus"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the Adenovirus dataset (DHF vs DF)

#TestingData_Adenovirus <- t(Expr_Adenovirus1)

PredVotes_Adenovirus <- predict(RF_Cerebral, newdata = TestingData_Adenovirus, type = "vote")
PredResponse_Adenovirus <- predict(RF_Cerebral, TestingData_Adenovirus, type="response")

ROCTest <- roc(ClassAdenoVsNormal, PredVotes_Adenovirus[,2], plot = F, print.auc = TRUE, levels = c("control", "AdenoVirus"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
