


rm(list = ls())

library(GEOquery)

#DengueDataset4 <- getGEO("GSE18090", GSEMatrix = T, AnnotGPL = T)
#DengueDataset4 <- DengueDataset4$GSE18090_series_matrix.txt.gz
 
#save(DengueDataset4, file = "./Data/DengueDataset4.rda")

load("./Data/DengueDataset4.rda")

Expr_Dengue4 <- exprs(DengueDataset4)
Pheno_Dengue4 <- pData(DengueDataset4)
FeatData_Dengue4 <- fData(DengueDataset4)


############################
## Annotation

## Expr_Dengue4
head(rownames(Expr_Dengue4))
rownames(Expr_Dengue4) <- FeatData_Dengue4$`Gene symbol`
summary(is.na(rownames(Expr_Dengue4)))
#rownames(Expr_Dengue4) <- gsub("-","", rownames(Expr_Dengue4))
#rownames(Expr_Dengue4) <- gsub("_","",rownames(Expr_Dengue4))
sel <- which(apply(Expr_Dengue4, 1, function(x) all(is.finite(x)) ))
Expr_Dengue4 <- Expr_Dengue4[sel, ]
Expr_Dengue4 <- Expr_Dengue4[!is.na(rownames(Expr_Dengue4)),]
dim(Expr_Dengue4)

range(Expr_Dengue4) 
Expr_Dengue4 <- log2(Expr_Dengue4)
Expr_Dengue4 <- t(scale(t(Expr_Dengue4), center = TRUE, scale = TRUE))


####################################


########################################################################################
#########################################################################################
###### DF and DHF vs non dengue

### Modify the phenotype

# Pheno1
Pheno_Dengue4$DiseaseStatus <- as.factor(Pheno_Dengue4$source_name_ch1)
levels(Pheno_Dengue4$DiseaseStatus) <- c("case", "case", "control") 
table(Pheno_Dengue4$DiseaseStatus)
Pheno_Dengue4$DiseaseStatus <- factor(Pheno_Dengue4$DiseaseStatus, levels = c("control", "case"))
#Expr_Dengue4 <- Expr_Dengue4[, colnames(Expr_Dengue4) %in% rownames(Pheno_Dengue4)]
all(rownames(Pheno_Dengue4) == colnames(Expr_Dengue4))

ClassDFvsnormal <- Pheno_Dengue4$DiseaseStatus

####################################
## Load the complicated malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (DHF vs DF)

TestingData_Dengue <- t(Expr_Dengue4)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDFvsnormal, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue4 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDFvsnormal)
sscurves_Dengue4
ROC_Dengue4 <- autoplot(sscurves_Dengue4, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE18090 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.43"), size = 5)
PRC_Dengue4 <- autoplot(sscurves_Dengue4, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE18090 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.725"), size = 5)

save(ROC_Dengue4, PRC_Dengue4, file = "./Objs/Dengue4_Curves.rda")


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the Dengue dataset (DSS vs DF)

#TestingData_Dengue <- t(Expr_Dengue4)

PredVotes_Dengue <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDFvsnormal, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest


######################################################################################
###### DF and DHF vs non dengue

### Modify the phenotype

# Pheno1
Pheno_Dengue4 <- Pheno_Dengue4[!(Pheno_Dengue4$source_name_ch1 == "PBMCs from ND patient"), ]
Pheno_Dengue4$DiseaseStatus <- as.factor(Pheno_Dengue4$source_name_ch1)
levels(Pheno_Dengue4$DiseaseStatus) <- c("DF", "DHF") 
table(Pheno_Dengue4$DiseaseStatus)

Expr_Dengue4 <- Expr_Dengue4[, colnames(Expr_Dengue4) %in% rownames(Pheno_Dengue4)]
all(rownames(Pheno_Dengue4) == colnames(Expr_Dengue4))

ClassDFvsDHF <- Pheno_Dengue4$DiseaseStatus

####################################
## Load the complicated malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (DF vs DHF)

TestingData_Dengue <- t(Expr_Dengue4)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDFvsDHF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
# sscurves_Dengue4 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDFvsDHF)
# sscurves_Dengue4
# ROC_Dengue4 <- autoplot(sscurves_Dengue4, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE18090 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.43"), size = 5)
# PRC_Dengue4 <- autoplot(sscurves_Dengue4, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE18090 (Dengue fever)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.725"), size = 5)
# 
# save(ROC_Dengue4, PRC_Dengue4, file = "./Objs/Dengue4_Curves.rda")


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the Dengue dataset (DSS vs DF)

#TestingData_Dengue <- t(Expr_Dengue4)

PredVotes_Dengue <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDFvsDHF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DHF"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest


