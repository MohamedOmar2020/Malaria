

rm(list = ls())

library(GEOquery)

#DengueDataset3 <- getGEO("GSE25001", GSEMatrix = T, AnnotGPL = T)
#DengueDataset3 <- DengueDataset3$GSE25001_series_matrix.txt.gz

#save(DengueDataset3, file = "./Data/DengueDataset3.rda")

load("./Data/DengueDataset3.rda")

Expr_Dengue3 <- exprs(DengueDataset3)
Pheno_Dengue3 <- pData(DengueDataset3)
FeatData_Dengue3 <- fData(DengueDataset3)


############################
## Annotation

## Expr_Dengue3
head(rownames(Expr_Dengue3))
rownames(Expr_Dengue3) <- FeatData_Dengue3$`Gene symbol`
summary(is.na(rownames(Expr_Dengue3)))
#rownames(Expr_Dengue3) <- gsub("-","", rownames(Expr_Dengue3))
#rownames(Expr_Dengue3) <- gsub("_","",rownames(Expr_Dengue3))
sel <- which(apply(Expr_Dengue3, 1, function(x) all(is.finite(x)) ))
Expr_Dengue3 <- Expr_Dengue3[sel, ]
Expr_Dengue3 <- Expr_Dengue3[!is.na(rownames(Expr_Dengue3)),]
dim(Expr_Dengue3)

range(Expr_Dengue3)
Expr_Dengue3 <- log2(Expr_Dengue3 + 65)
#plot(density(Expr_Dengue3))
#boxplot(Expr_Dengue3)
# X1 <- Expr_Dengue3
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt1 <- genefilter(2^X1,ffun)
# Expr_Dengue3 <- Expr_Dengue3[filt1,]
# 
Expr_Dengue3 <- t(scale(t(Expr_Dengue3), center = TRUE, scale = TRUE))


####################################


########################################################################################
#########################################################################################
###### DSS vs DF

### Modify the phenotype

# Pheno1
Pheno_Dengue3$DiseaseStatus2 <- as.factor(Pheno_Dengue3$`disease state:ch1`)
levels(Pheno_Dengue3$DiseaseStatus2) <- c("DSS", "DF") 
table(Pheno_Dengue3$DiseaseStatus2)
Pheno_Dengue3$DiseaseStatus2 <- factor(Pheno_Dengue3$DiseaseStatus2, levels = c("DF", "DSS"))
#Expr_Dengue3 <- Expr_Dengue3[, colnames(Expr_Dengue3) %in% rownames(Pheno_Dengue3)]
all(rownames(Pheno_Dengue3) == colnames(Expr_Dengue3))

ClassDSSvsDF <- Pheno_Dengue3$DiseaseStatus2

####################################
## Load the complicated malaria signature
load("./Objs/RF_Comp.rda")

#################
## Predict in the Dengue dataset (DHF vs DF)

TestingData_Dengue <- t(Expr_Dengue3)

PredVotes_Dengue <- predict(RF_Comp, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Comp, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDSSvsDF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DSS"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_Dengue3 <- evalmod(scores = PredVotes_Dengue[,2], labels = ClassDSSvsDF)
sscurves_Dengue3
ROC_Dengue3 <- autoplot(sscurves_Dengue3, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE25001 (DF vs DSS)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.48"), size = 5)
PRC_Dengue3 <- autoplot(sscurves_Dengue3, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE25001 (DF vs DSS)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.32"), size = 5)

save(ROC_Dengue3, PRC_Dengue3, file = "./Objs/Dengue3_Curves.rda")


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the Dengue dataset (DSS vs DF)

#TestingData_Dengue <- t(Expr_Dengue3)

PredVotes_Dengue <- predict(RF_Cerebral, newdata = TestingData_Dengue, type = "vote")
PredResponse_Dengue <- predict(RF_Cerebral, TestingData_Dengue, type="response")

ROCTest <- roc(ClassDSSvsDF, PredVotes_Dengue[,2], plot = F, print.auc = TRUE, levels = c("DF", "DSS"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
