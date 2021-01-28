

rm(list = ls())

library(GEOquery)

#GSE7586 <- getGEO("GSE7586", GSEMatrix = T, AnnotGPL = T)
#GSE7586 <- GSE7586$GSE7586_series_matrix.txt.gz

#save(GSE7586, file = "./Data/PlacentalMalariaDataset.rda")

load("./Data/PlacentalMalariaDataset.rda")

Expr_Test2 <- exprs(GSE7586)
Pheno_Test2 <- pData(GSE7586)
FeatData_Test2 <- fData(GSE7586)


############################
## Annotation

## Expr_Cerebral
head(rownames(Expr_Test2))
rownames(Expr_Test2) <- FeatData_Test2$`Gene symbol`
summary(is.na(rownames(Expr_Test2)))
#rownames(Expr_Cerebral) <- gsub("-","", rownames(Expr_Cerebral))
#rownames(Expr_Cerebral) <- gsub("_","",rownames(Expr_Cerebral))
sel <- which(apply(Expr_Test2, 1, function(x) all(is.finite(x)) ))
Expr_Test2 <- Expr_Test2[sel, ]
Expr_Test2 <- Expr_Test2[!is.na(rownames(Expr_Test2)),]
dim(Expr_Test2)

range(Expr_Test2)
#plot(density(Expr_Cerebral))
#boxplot(Expr_Cerebral)
# X1 <- Expr_Cerebral
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt1 <- genefilter(2^X1,ffun)
# Expr_Cerebral <- Expr_Cerebral[filt1,]
# 
Expr_Test2 <- t(scale(t(Expr_Test2), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype
# Placental malaria +ve vs PM-ve

Pheno_Test2$`PM Status:ch1`[is.na(Pheno_Test2$`PM Status:ch1`)] <- "uninfected"
Pheno_Test2$DiseaseStatus <- as.factor(Pheno_Test2$`PM Status:ch1`)
levels(Pheno_Test2$DiseaseStatus) <- c("Complicated", "unComplicated") 
table(Pheno_Test2$DiseaseStatus)
Pheno_Test2$DiseaseStatus <- factor(Pheno_Test2$DiseaseStatus, levels = c("unComplicated", "Complicated"))

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Cerebral)]
#all(rownames(Pheno_Cerebral) == colnames(expr1))


ClassComplicatedVSunComplicated <- Pheno_Test2$DiseaseStatus


###################################
####################################

### Modify the phenotype
# Placental malaria +ve vs PM-ve

Pheno_Test2$`PM Status:ch1`[is.na(Pheno_Test2$`PM Status:ch1`)] <- "uninfected"
Pheno_Test2$InflammationStatus <- as.factor(Pheno_Test2$`Inflammation:ch1`)
levels(Pheno_Test2$InflammationStatus) <- c("No", "Yes") 
table(Pheno_Test2$InflammationStatus)
#Pheno_Test2$InflammationStatus <- factor(Pheno_Test2$DiseaseStatus, levels = c("unComplicated", "Complicated"))

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Cerebral)]
#all(rownames(Pheno_Cerebral) == colnames(expr1))


ClassInflammation <- Pheno_Test2$InflammationStatus

save(Expr_Test2, ClassComplicatedVSunComplicated, ClassInflammation, file = "./Objs/PlacentalMalaria.rda")






