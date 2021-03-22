rm(list = ls())

library(GEOquery)

# TB4 <- getGEO("GSE83456", GSEMatrix = T, AnnotGPL = T)
# TB4 <- TB4$GSE83456_series_matrix.txt.gz
# 
# save(TB4, file = "./Data/TB4.rda")
# 
load("./Data/TB4.rda")

Expr_TB4 <- exprs(TB4)
Pheno_TB4 <- pData(TB4)
FeatData_TB4 <- fData(TB4)

############################
## Annotation

## Expr_TB4
head(rownames(Expr_TB4))
rownames(Expr_TB4) <- FeatData_TB4$`Gene symbol`
summary(is.na(rownames(Expr_TB4)))
#rownames(Expr_TB4) <- gsub("-","", rownames(Expr_TB4))
#rownames(Expr_TB4) <- gsub("_","",rownames(Expr_TB4))
sel <- which(apply(Expr_TB4, 1, function(x) all(is.finite(x)) ))
Expr_TB4 <- Expr_TB4[sel, ]
Expr_TB4 <- Expr_TB4[!is.na(rownames(Expr_TB4)),]
dim(Expr_TB4)

range(Expr_TB4)
#Expr_TB4 <- t(scale(t(Expr_TB4), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# TB vs Healthy

# Pheno1

# Remove sarcoidosis
Pheno_TB4 <- Pheno_TB4[!(Pheno_TB4$`disease state:ch1`) == "Sarcoid", ]
Pheno_TB4$DiseaseStatus <- as.factor(Pheno_TB4$`disease state:ch1`)
levels(Pheno_TB4$DiseaseStatus) <- c("case", "control", "case")
Pheno_TB4$DiseaseStatus <- factor(Pheno_TB4$DiseaseStatus, levels = c("control", "case"))
table(Pheno_TB4$DiseaseStatus)


Expr_TB4 <- Expr_TB4[, colnames(Expr_TB4) %in% rownames(Pheno_TB4)]
all(rownames(Pheno_TB4) == colnames(Expr_TB4))

ClassTBVsHealthy<- Pheno_TB4$DiseaseStatus

####################################
## Load the model
load("./Objs/RF_Comp.rda")

#################
## Predict in the WestNile dataset (WestNile vs normal)

TestingData_TB4 <- t(Expr_TB4)

PredVotes_TB4 <- predict(RF_Comp, newdata = TestingData_TB4, type = "vote")
PredResponse_TB4 <- predict(RF_Comp, TestingData_TB4, type="response")

ROCTest <- roc(ClassTBVsHealthy, PredVotes_TB4[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# For ROC and PRC curves
sscurves_TB4 <- evalmod(scores = PredVotes_TB4[,2], labels = ClassTBVsHealthy)
sscurves_TB4
ROC_TB4 <- autoplot(sscurves_TB4, curvetype = c("ROC")) + labs(title = "ROC curve of the complicated malaria signature in GSE83456 (pulomnary and extrapulmonary TB vs healthy)") + annotate("text", x = .65, y = .25, label = paste("AUC = 0.50"), size = 4)
PRC_TB4 <- autoplot(sscurves_TB4, curvetype = c("PRC")) + labs(title = "PRC curve of the complicated malaria signature in GSE83456 (pulomnary and extrapulmonary TB vs healthy)") + annotate("text", x = .65, y = .25, label = paste("AUPRC = 0.61"), size = 4)

save(ROC_TB4, PRC_TB4, file = "./Objs/TB4_Curves.rda")

########################################################################################


####################################
## Load the cerebral malaria
load("./Objs/RF_Cerebral.rda")

#################
## Predict in the WestNile dataset (DHF vs DF)

#TestingData_ManyInfections <- t(Expr_TB4)

PredVotes_TB4 <- predict(RF_Cerebral, newdata = TestingData_TB4, type = "vote")
PredResponse_TB4 <- predict(RF_Cerebral, TestingData_TB4, type="response")

ROCTest <- roc(ClassTBVsHealthy, PredVotes_TB4[,2], plot = F, print.auc = TRUE, levels = c("control", "case"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest
