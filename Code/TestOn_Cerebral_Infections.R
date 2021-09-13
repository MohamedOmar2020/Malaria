

rm(list = ls())

library(GEOquery)
library(randomForest)
library(RRF)
library(precrec)
library(ggplot2)

MeningitisDataset2 <- getGEO("GSE80496", GSEMatrix = T, AnnotGPL = T)
MeningitisDataset2 <- MeningitisDataset2$GSE80496_series_matrix.txt.gz

MeningitisDataset4 <- getGEO("GSE40586", GSEMatrix = T, AnnotGPL = T)
MeningitisDataset4 <- MeningitisDataset4$GSE40586_series_matrix.txt.gz

  
save(MeningitisDataset2, MeningitisDataset4, file = "./Data/MeningitisDatasets.rda")
load("./Data/MeningitisDatasets.rda")

############################################
## Get the expression matrices
Expr_Meningitis2 <- exprs(MeningitisDataset2)
Expr_Meningitis4 <- exprs(MeningitisDataset4)

############################################
## Get the phenotype
Pheno_Meningitis2 <- pData(MeningitisDataset2)
Pheno_Meningitis4 <- pData(MeningitisDataset4)


###########################################
## Get the feature data for annotation
FeatData_Meningitis2 <- fData(MeningitisDataset2)
FeatData_Meningitis4 <- fData(MeningitisDataset4)

############################
## Annotation
##############################################
## Expr_Meningitis2
head(rownames(Expr_Meningitis2))
rownames(Expr_Meningitis2) <- FeatData_Meningitis2$`Gene symbol`
summary(is.na(rownames(Expr_Meningitis2)))
sel <- which(apply(Expr_Meningitis2, 1, function(x) all(is.finite(x)) ))
Expr_Meningitis2 <- Expr_Meningitis2[sel, ]
Expr_Meningitis2 <- Expr_Meningitis2[!is.na(rownames(Expr_Meningitis2)),]
dim(Expr_Meningitis2)
range(Expr_Meningitis2)
Expr_Meningitis2 <- t(scale(t(Expr_Meningitis2), center = TRUE, scale = TRUE))

##############################################
## Expr_Meningitis4
head(rownames(Expr_Meningitis4))
rownames(Expr_Meningitis4) <- FeatData_Meningitis4$`Gene symbol`
summary(is.na(rownames(Expr_Meningitis4)))
sel <- which(apply(Expr_Meningitis4, 1, function(x) all(is.finite(x)) ))
Expr_Meningitis4 <- Expr_Meningitis4[sel, ]
Expr_Meningitis4 <- Expr_Meningitis4[!is.na(rownames(Expr_Meningitis4)),]
dim(Expr_Meningitis4)
range(Expr_Meningitis4)
Expr_Meningitis4 <- t(scale(t(Expr_Meningitis4), center = TRUE, scale = TRUE))

####################################
### Modify the phenotype
###################################
# Pheno2
table(Pheno_Meningitis2$`disease state:ch1`)
Pheno_Meningitis2$DiseaseStatus <- as.factor(Pheno_Meningitis2$`disease state:ch1`)
levels(Pheno_Meningitis2$DiseaseStatus) <- c("control", "case", "case") 
table(Pheno_Meningitis2$DiseaseStatus)

ClassMeningitisVsNormal <- Pheno_Meningitis2$DiseaseStatus

###################################
# Pheno4
table(Pheno_Meningitis4$`sample group:ch1`)
Pheno_Meningitis4$DiseaseStatus <- as.factor(Pheno_Meningitis4$`sample group:ch1`)
levels(Pheno_Meningitis4$DiseaseStatus) <- c("case", "control")
Pheno_Meningitis4$DiseaseStatus <- factor(Pheno_Meningitis4$DiseaseStatus, levels = c('control', 'case'))
table(Pheno_Meningitis4$DiseaseStatus)

ClassMeningitisVsNormal4 <- Pheno_Meningitis4$DiseaseStatus


####################################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")
#################
## Predict in the 2nd cerebral dataset meningitis vs normal) using the two malaria signature
sel_severe <- intersect(rownames(RF_Comp$importance), rownames(Expr_Meningitis2))
sel_cerebral <- intersect(rownames(RF_Cerebral$importance), rownames(Expr_Meningitis2))

# Subset the severe signature
RF_Comp$importance <- RF_Comp$importance[sel_severe, ]
RF_Comp$importanceSD <- RF_Comp$importanceSD[sel_severe, ]
RF_Comp$forest$ncat <- RF_Comp$forest$ncat[sel_severe]

# Subset the cerebral signature
RF_Cerebral$importance <- RF_Cerebral$importance[sel_cerebral, ]
RF_Cerebral$importanceSD <- RF_Cerebral$importanceSD[sel_cerebral, ]
RF_Cerebral$forest$ncat <- RF_Cerebral$forest$ncat[sel_cerebral]

# transpose the expression matrix
TestingData_Meningitis2 <- t(Expr_Meningitis2)

PredVotes_Meningitis2_severe <- predict(RF_Comp, newdata = TestingData_Meningitis2, type = "vote")
PredVotes_Meningitis2_cerebral <- predict(RF_Cerebral, newdata = TestingData_Meningitis2, type = "vote")

# For ROC and PRC curves
# Severe signauture:
sscurves_Meningitis2_severe <- evalmod(scores = PredVotes_Meningitis2_severe[,2], labels = ClassMeningitisVsNormal)
sscurves_Meningitis2_severe
ROC_Meningitis2_severe <- autoplot(sscurves_Meningitis2_severe, curvetype = c("ROC")) + labs(title = "GSE80496 (meningitis vs healthy controls)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.53"), size = 3)

# cerebral signauture:
sscurves_Meningitis2_cerebral <- evalmod(scores = PredVotes_Meningitis2_cerebral[,2], labels = ClassMeningitisVsNormal)
sscurves_Meningitis2_cerebral
ROC_Meningitis2_cerebral <- autoplot(sscurves_Meningitis2_cerebral, curvetype = c("ROC")) + labs(title = "GSE80496 (meningitis vs healthy controls)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.12"), size = 3)

####################################
#################
## Load the severe malaria signature
load("./Objs/RF_Comp.rda")

## Load the cerebral malaria signature
load("./Objs/RF_Cerebral.rda")

# transpose the expression matrix
TestingData_Meningitis4 <- t(Expr_Meningitis4)

PredVotes_Meningitis4_severe <- predict(RF_Comp, newdata = TestingData_Meningitis4, type = "vote")
PredVotes_Meningitis4_cerebral <- predict(RF_Cerebral, newdata = TestingData_Meningitis4, type = "vote")

# For ROC and PRC curves
# severe signature:
sscurves_Meningitis4_severe <- evalmod(scores = PredVotes_Meningitis4_severe[,2], labels = ClassMeningitisVsNormal4)
sscurves_Meningitis4_severe
ROC_Meningitis4_severe <- autoplot(sscurves_Meningitis4_severe, curvetype = c("ROC")) + labs(title = "GSE40586 (meningitis vs healthy controls)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.75"), size = 3)

# cerebral signature:
sscurves_Meningitis4_cerebral <- evalmod(scores = PredVotes_Meningitis4_cerebral[,2], labels = ClassMeningitisVsNormal4)
sscurves_Meningitis4_cerebral
ROC_Meningitis4_cerebral <- autoplot(sscurves_Meningitis4_cerebral, curvetype = c("ROC")) + labs(title = "GSE40586 (meningitis vs healthy controls)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.30"), size = 3)

#################################################
### Save
save(ROC_Meningitis2_severe, ROC_Meningitis4_severe, file = "./Objs/Meningitis_Curves_severe.rda")
save(ROC_Meningitis2_cerebral, ROC_Meningitis4_cerebral, file = "./Objs/Meningitis_Curves_cerebral.rda")






