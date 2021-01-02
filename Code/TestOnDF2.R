


rm(list = ls())

library(GEOquery)

DengueDataset2 <- getGEO("GSE96656", GSEMatrix = T, AnnotGPL = T)
DengueDataset2 <- DengueDataset2$GSE96656_series_matrix.txt.gz


Expr_Dengue2 <- exprs(DengueDataset2)
Pheno_Dengue2 <- pData(DengueDataset2)
FeatData_Dengue2 <- fData(DengueDataset2)


############################
## Annotation

## Expr_Dengue2
head(rownames(Expr_Dengue2))
rownames(Expr_Dengue2) <- FeatData_Dengue2$ORF
summary(is.na(rownames(Expr_Dengue2)))
#rownames(Expr_Dengue2) <- gsub("-","", rownames(Expr_Dengue2))
#rownames(Expr_Dengue2) <- gsub("_","",rownames(Expr_Dengue2))
sel <- which(apply(Expr_Dengue2, 1, function(x) all(is.finite(x)) ))
Expr_Dengue2 <- Expr_Dengue2[sel, ]
Expr_Dengue2 <- Expr_Dengue2[!is.na(rownames(Expr_Dengue2)),]
dim(Expr_Dengue2)

range(Expr_Dengue2)
#plot(density(Expr_Dengue2))
#boxplot(Expr_Dengue2)
# X1 <- Expr_Dengue2
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt1 <- genefilter(2^X1,ffun)
# Expr_Dengue2 <- Expr_Dengue2[filt1,]
# 
#Expr_Dengue2 <- t(scale(t(Expr_Dengue2), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype
# Remove controls

# Pheno1
Pheno_Dengue2$DiseaseStatus <- as.factor(Pheno_Dengue2$`disease state:ch2`)
levels(Pheno_Dengue2$DiseaseStatus) <- c("case", "case", "control") 
table(Pheno_Dengue2$DiseaseStatus)

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Dengue2)]
#all(rownames(Pheno_Dengue2) == colnames(expr1))

ClassDengueVsNormal <- Pheno_Dengue2$DiseaseStatus

####################################
## Load the model
load("./Objs/RF_Comp.rda")

## Some features (2) are present in the RF model but not in the expression matrix >> removed them
CommonGns <- intersect(rownames(Expr_Dengue2), rownames(RF_Comp$importance))
RF_Comp$importance <- RF_Comp$importance[CommonGns, ]

#################
## Predict in the Dengue dataset (Dengue vs normal)

TestingData_Dengue <- t(Expr_Dengue2)

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
Pheno_Dengue2 <- Pheno_Dengue2[!(Pheno_Dengue2$`disease state:ch2` == "Healthy"), ]
Pheno_Dengue2$DiseaseStatus2 <- as.factor(Pheno_Dengue2$`disease state:ch2`)
levels(Pheno_Dengue2$DiseaseStatus2) <- c("DF", "DHF") 
table(Pheno_Dengue2$DiseaseStatus2)

Expr_Dengue2 <- Expr_Dengue2[, colnames(Expr_Dengue2) %in% rownames(Pheno_Dengue2)]
all(rownames(Pheno_Dengue2) == colnames(Expr_Dengue2))

ClassDHFvsDF <- Pheno_Dengue2$DiseaseStatus2

####################################
## Load the complicated malaria signature
#load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (DHF vs DF)

TestingData_Dengue <- t(Expr_Dengue2)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDHFvsDF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

## Some features (2) are present in the RF model but not in the expression matrix >> removed them
CommonGns <- intersect(rownames(Expr_Dengue2), rownames(RF_Cerebral$importance))
RF_Cerebral$importance <- RF_Cerebral$importance[CommonGns, ]


#################
## Predict in the Dengue dataset (DHF vs DF)

#TestingData_Dengue <- t(Expr_Dengue2)

PredVotes_Dengue <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDHFvsDF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
