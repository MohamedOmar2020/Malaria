

rm(list = ls())

library(GEOquery)


MeningitisDataset1 <- getGEO("GSE83892", GSEMatrix = T, AnnotGPL = T)
MeningitisDataset1 <- MeningitisDataset1$GSE83892_series_matrix.txt.gz

MeningitisDataset2 <- getGEO("GSE80496", GSEMatrix = T, AnnotGPL = T)
MeningitisDataset2 <- MeningitisDataset2$GSE80496_series_matrix.txt.gz

MeningitisDataset3 <- getGEO("GSE72829", GSEMatrix = T, AnnotGPL = T)
MeningitisDataset3_1 <- MeningitisDataset3$`GSE72829-GPL10558_series_matrix.txt.gz`
MeningitisDataset3_2 <- MeningitisDataset3$`GSE72829-GPL6883_series_matrix.txt.gz`
MeningitisDataset3_3 <- MeningitisDataset3$`GSE72829-GPL6947_series_matrix.txt.gz`

MeningitisDataset4 <- getGEO("GSE40586", GSEMatrix = T, AnnotGPL = T)
MeningitisDataset4 <- MeningitisDataset4$GSE40586_series_matrix.txt.gz

  
save(MeningitisDataset1, MeningitisDataset2, MeningitisDataset3_1, MeningitisDataset3_2, MeningitisDataset3_3, MeningitisDataset4, file = "./Data/MeningitisDatasets.rda")
load("./Data/MeningitisDatasets.rda")

############################################
## Get the expression matrices
Expr_Meningitis1 <- exprs(MeningitisDataset1)
Expr_Meningitis2 <- exprs(MeningitisDataset2)
Expr_Meningitis3_1 <- exprs(MeningitisDataset3_1)
Expr_Meningitis3_2 <- exprs(MeningitisDataset3_2)
Expr_Meningitis3_3 <- exprs(MeningitisDataset3_3)
Expr_Meningitis4 <- exprs(MeningitisDataset4)

############################################
## Get the phenotype
Pheno_Meningitis1 <- pData(MeningitisDataset1)
Pheno_Meningitis2 <- pData(MeningitisDataset2)
Pheno_Meningitis3_1 <- pData(MeningitisDataset3_1)
Pheno_Meningitis3_2 <- pData(MeningitisDataset3_2)
Pheno_Meningitis3_3 <- pData(MeningitisDataset3_3)
Pheno_Meningitis4 <- pData(MeningitisDataset4)


###########################################
## Get the feature data for annotation
FeatData_Meningitis1 <- fData(MeningitisDataset1)
FeatData_Meningitis2 <- fData(MeningitisDataset2)
FeatData_Meningitis3_1 <- fData(MeningitisDataset3_1)
FeatData_Meningitis3_2 <- fData(MeningitisDataset3_2)
FeatData_Meningitis3_3 <- fData(MeningitisDataset3_3)
FeatData_Meningitis4 <- fData(MeningitisDataset4)

############################
## Annotation
## Expr_Meningitis1
head(rownames(Expr_Meningitis1))
rownames(Expr_Meningitis1) <- FeatData_Meningitis1$`Gene symbol`
summary(is.na(rownames(Expr_Meningitis1)))
#sel <- which(apply(Expr_Meningitis1, 1, function(x) all(is.finite(x)) ))
#Expr_Meningitis1 <- Expr_Meningitis1[sel, ]
Expr_Meningitis1 <- Expr_Meningitis1[!is.na(rownames(Expr_Meningitis1)), ]
dim(Expr_Meningitis1)
range(Expr_Meningitis1)
Expr_Meningitis1 <- t(scale(t(Expr_Meningitis1), center = TRUE, scale = TRUE))

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
## Expr_Meningitis3_1
head(rownames(Expr_Meningitis3_1))
rownames(Expr_Meningitis3_1) <- FeatData_Meningitis3_1$`Gene symbol`
summary(is.na(rownames(Expr_Meningitis3_1)))
sel <- which(apply(Expr_Meningitis3_1, 1, function(x) all(is.finite(x)) ))
Expr_Meningitis3_1 <- Expr_Meningitis3_1[sel, ]
Expr_Meningitis3_1 <- Expr_Meningitis3_1[!is.na(rownames(Expr_Meningitis3_1)),]
dim(Expr_Meningitis3_1)
range(Expr_Meningitis3_1)
Expr_Meningitis3_1 <- log2(Expr_Meningitis3_1+55)
Expr_Meningitis3_1 <- t(scale(t(Expr_Meningitis3_1), center = TRUE, scale = TRUE))

## Expr_Meningitis3_2
head(rownames(Expr_Meningitis3_2))
rownames(Expr_Meningitis3_2) <- FeatData_Meningitis3_2$`Gene symbol`
summary(is.na(rownames(Expr_Meningitis3_2)))
sel <- which(apply(Expr_Meningitis3_2, 1, function(x) all(is.finite(x)) ))
Expr_Meningitis3_2 <- Expr_Meningitis3_2[sel, ]
Expr_Meningitis3_2 <- Expr_Meningitis3_2[!is.na(rownames(Expr_Meningitis3_2)),]
dim(Expr_Meningitis3_2)
range(Expr_Meningitis3_2)
Expr_Meningitis3_2 <- t(scale(t(Expr_Meningitis3_2), center = TRUE, scale = TRUE))

## Expr_Meningitis3_3
head(rownames(Expr_Meningitis3_3))
rownames(Expr_Meningitis3_3) <- FeatData_Meningitis3_3$`Gene symbol`
summary(is.na(rownames(Expr_Meningitis3_3)))
sel <- which(apply(Expr_Meningitis3_3, 1, function(x) all(is.finite(x)) ))
Expr_Meningitis3_3 <- Expr_Meningitis3_3[sel, ]
Expr_Meningitis3_3 <- Expr_Meningitis3_3[!is.na(rownames(Expr_Meningitis3_3)),]
dim(Expr_Meningitis3_3)
range(Expr_Meningitis3_3)
Expr_Meningitis3_3 <- t(scale(t(Expr_Meningitis3_3), center = TRUE, scale = TRUE))

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

# Pheno1
table(Pheno_Meningitis1$`group:ch1`)
Pheno_Meningitis1$DiseaseStatus <- as.factor(Pheno_Meningitis1$`group:ch1`)
levels(Pheno_Meningitis1$DiseaseStatus) <- c("control", "case", "case") 
table(Pheno_Meningitis1$DiseaseStatus)

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Dengue1)]
#all(rownames(Pheno_Dengue1) == colnames(expr1))

ClassTBMVsNormal <- Pheno_Meningitis1$DiseaseStatus

###################################
# Pheno2
table(Pheno_Meningitis2$`disease state:ch1`)
Pheno_Meningitis2$DiseaseStatus <- as.factor(Pheno_Meningitis2$`disease state:ch1`)
levels(Pheno_Meningitis2$DiseaseStatus) <- c("control", "case", "case") 
table(Pheno_Meningitis2$DiseaseStatus)

ClassMeningitisVsNormal <- Pheno_Meningitis2$DiseaseStatus

###################################
# Pheno3_1
table(Pheno_Meningitis3_1)
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
## Predict in the 1st cerebral dataset (TB-meningitis vs normal) using the two malaria signature
TestingData_Meningitis1 <- t(Expr_Meningitis1)

PredVotes_Meningitis1_severe <- predict(RF_Comp, newdata = TestingData_Meningitis1, type = "vote")
PredVotes_Meningitis1_cerebral <- predict(RF_Cerebral, newdata = TestingData_Meningitis1, type = "vote")

ROCTest <- roc(ClassDengueVsNormal, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue1 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDengueVsNormal)
sscurves_Dengue1
ROC_Dengue <- autoplot(sscurves_Dengue1, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE51808 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.52"), size = 3)
PRC_Dengue <- autoplot(sscurves_Dengue1, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE51808 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.54"), size = 3)

save(ROC_Dengue, PRC_Dengue, file = "./Objs/Dengue1_Curves.rda")

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

ROC_meningitis2_severe <- roc(ClassMeningitisVsNormal, PredVotes_Meningitis2_severe[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROC_meningitis2_severe

ROC_meningitis2_cerebral <- roc(ClassMeningitisVsNormal, PredVotes_Meningitis2_cerebral[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROC_meningitis2_cerebral

# For ROC and PRC curves
# Severe signauture:
sscurves_Meningitis2_severe <- evalmod(scores = PredVotes_Meningitis2_severe[,2], labels = ClassMeningitisVsNormal)
sscurves_Meningitis2_severe
ROC_Meningitis2_severe <- autoplot(sscurves_Meningitis2_severe, curvetype = c("ROC")) + labs(title = "ROC curve of the severe malaria signature in GSE80496 (Meningococal disease vs healthy controls)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.53"), size = 3)

# cerebral signauture:
sscurves_Meningitis2_cerebral <- evalmod(scores = PredVotes_Meningitis2_cerebral[,2], labels = ClassMeningitisVsNormal)
sscurves_Meningitis2_cerebral
ROC_meningitis2_cerebral <- autoplot(sscurves_Meningitis2_cerebral, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE80496 (Meningococal disease vs healthy controls)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.53"), size = 3)

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

ROC_meningitis4_severe <- roc(ClassMeningitisVsNormal4, PredVotes_Meningitis4_severe[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROC_meningitis4_severe

ROC_meningitis4_cerebral <- roc(ClassMeningitisVsNormal4, PredVotes_Meningitis4_cerebral[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROC_meningitis4_cerebral

# For ROC and PRC curves
# severe signature:
sscurves_Meningitis4_severe <- evalmod(scores = PredVotes_Meningitis4_severe[,2], labels = ClassMeningitisVsNormal4)
sscurves_Meningitis4_severe
ROC_Meningitis4_severe <- autoplot(sscurves_Meningitis4_severe, curvetype = c("ROC")) + labs(title = "ROC curve of the severe malaria signature in GSE40586 (Meningococal disease vs healthy controls)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.75"), size = 3)

# cerebral signature:
sscurves_Meningitis4_cerebral <- evalmod(scores = PredVotes_Meningitis4_cerebral[,2], labels = ClassMeningitisVsNormal4)
sscurves_Meningitis4_cerebral
ROC_Meningitis4_cerebral <- autoplot(sscurves_Meningitis4_severe, curvetype = c("ROC")) + labs(title = "ROC curve of the cerebral malaria signature in GSE40586 (Meningococal disease vs healthy controls)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.30"), size = 3)








