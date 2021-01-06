

rm(list = ls())

library(GEOquery)

GSE116306 <- getGEO("GSE116306", GSEMatrix = T, AnnotGPL = T)
GSE116306 <- GSE116306$GSE116306_series_matrix.txt.gz

save(GSE116306, file = "./Data/ExtraMalariaDataset.rda")

load("./Data/ExtraMalariaDataset.rda")

Expr_Test2 <- exprs(GSE116306)
Pheno_Test2 <- pData(GSE116306)
FeatData_Test2 <- fData(GSE116306)


############################
## Annotation

## Expr_Cerebral
head(rownames(Expr_Test2))
rownames(Expr_Test2) <- FeatData_Test2$GENE_SYMBOL
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

# Cerebral vs non-cerebral
# Pheno1
Pheno_Test2$DiseaseStatus <- as.factor(Pheno_Test2$`disease state:ch1`)
levels(Pheno_Test2$DiseaseStatus) <- c("cerebral", "nonCerebral", "nonCerebral") 
table(Pheno_Test2$DiseaseStatus)
Pheno_Test2$DiseaseStatus <- factor(Pheno_Test2$DiseaseStatus, levels = c("nonCerebral", "cerebral"))

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Cerebral)]
#all(rownames(Pheno_Cerebral) == colnames(expr1))

ClassCerebralVsNonCerebral <- Pheno_Test2$DiseaseStatus

save(Expr_Test2, ClassCerebralVsNonCerebral, file = "./Objs/CerebralExtraValidation.rda")


##################################
# Cerebral vs non-cerebral
# Pheno1
Pheno_Test2$DiseaseStatus <- as.factor(Pheno_Test2$`disease state:ch1`)
levels(Pheno_Test2$DiseaseStatus) <- c("Complicated", "unComplicated", "Complicated") 
table(Pheno_Test2$DiseaseStatus)
Pheno_Test2$DiseaseStatus <- factor(Pheno_Test2$DiseaseStatus, levels = c("unComplicated", "Complicated"))

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Cerebral)]
#all(rownames(Pheno_Cerebral) == colnames(expr1))

ClassComplicatedVSunComplicated <- Pheno_Test2$DiseaseStatus

save(Expr_Test2, ClassComplicatedVSunComplicated, file = "./Objs/CompExtraValidation.rda")

