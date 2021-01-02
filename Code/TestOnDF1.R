


rm(list = ls())

library(GEOquery)

DengueDataset1 <- getGEO("GSE51808", GSEMatrix = T, AnnotGPL = T)
DengueDataset1 <- DengueDataset1$GSE51808_series_matrix.txt.gz


Expr_Dengue1 <- exprs(DengueDataset1)
Pheno_Dengue1 <- pData(DengueDataset1)
FeatData_Dengue1 <- fData(DengueDataset1)


############################
## Annotation

## Expr_Dengue1
head(rownames(Expr_Dengue1))
rownames(Expr_Dengue1) <- FeatData_Dengue1$`Gene symbol`
summary(is.na(rownames(Expr_Dengue1)))
#rownames(Expr_Dengue1) <- gsub("-","", rownames(Expr_Dengue1))
#rownames(Expr_Dengue1) <- gsub("_","",rownames(Expr_Dengue1))
sel <- which(apply(Expr_Dengue1, 1, function(x) all(is.finite(x)) ))
Expr_Dengue1 <- Expr_Dengue1[sel, ]
Expr_Dengue1 <- Expr_Dengue1[!is.na(rownames(Expr_Dengue1)),]
dim(Expr_Dengue1)

range(Expr_Dengue1)
#plot(density(Expr_Dengue1))
#boxplot(Expr_Dengue1)
# X1 <- Expr_Dengue1
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt1 <- genefilter(2^X1,ffun)
# Expr_Dengue1 <- Expr_Dengue1[filt1,]
# 
Expr_Dengue1 <- t(scale(t(Expr_Dengue1), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype
# Remove controls

# Pheno1
Pheno_Dengue1$DiseaseStatus <- as.factor(Pheno_Dengue1$`status:ch1`)
levels(Pheno_Dengue1$DiseaseStatus) <- c("control", "control", "case", "case") 
table(Pheno_Dengue1$DiseaseStatus)

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Dengue1)]
#all(rownames(Pheno_Dengue1) == colnames(expr1))

ClassDengueVsNormal <- Pheno_Dengue1$DiseaseStatus

####################################
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (Dengue vs normal)

TestingData_Dengue <- t(Expr_Dengue1)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDengueVsNormal, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest


########################################################################################
#########################################################################################
###### HDF vs DF

### Modify the phenotype
# Remove controls

# Pheno1
Pheno_Dengue1 <- Pheno_Dengue1[!(Pheno_Dengue1$`status:ch1` %in% c("convalescent", "control")), ]
Pheno_Dengue1$DiseaseStatus2 <- as.factor(Pheno_Dengue1$`status:ch1`)
#levels(Pheno_Dengue1$DiseaseStatus) <- c("DF", "DHF") 
table(Pheno_Dengue1$DiseaseStatus2)

Expr_Dengue1 <- Expr_Dengue1[, colnames(Expr_Dengue1) %in% rownames(Pheno_Dengue1)]
all(rownames(Pheno_Dengue1) == colnames(Expr_Dengue1))

ClassDHFvsDF <- Pheno_Dengue1$DiseaseStatus2

####################################
## Load the complicated malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (DHF vs DF)

TestingData_Dengue <- t(Expr_Dengue1)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDHFvsDF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the Dengue dataset (DHF vs DF)

#TestingData_Dengue <- t(Expr_Dengue1)

PredVotes_Dengue <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDHFvsDF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
