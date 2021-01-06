


library(GEOquery)

CerebralDataset <- getGEO("GSE116306", GSEMatrix = T, AnnotGPL = T)
CerebralDataset <- CerebralDataset$GSE116306_series_matrix.txt.gz

save(CerebralDataset, file = "./Data/CerebralDataset.rda")

Expr_Cerebral <- exprs(CerebralDataset)
Pheno_Cerebral <- pData(CerebralDataset)
FeatData_Cerebral <- fData(CerebralDataset)


############################
## Annotation

## Expr_Cerebral
head(rownames(Expr_Cerebral))
rownames(Expr_Cerebral) <- FeatData_Cerebral$GENE_SYMBOL
summary(is.na(rownames(Expr_Cerebral)))
#rownames(Expr_Cerebral) <- gsub("-","", rownames(Expr_Cerebral))
#rownames(Expr_Cerebral) <- gsub("_","",rownames(Expr_Cerebral))
sel <- which(apply(Expr_Cerebral, 1, function(x) all(is.finite(x)) ))
Expr_Cerebral <- Expr_Cerebral[sel, ]
Expr_Cerebral <- Expr_Cerebral[!is.na(rownames(Expr_Cerebral)),]
dim(Expr_Cerebral)

range(Expr_Cerebral)
#plot(density(Expr_Cerebral))
#boxplot(Expr_Cerebral)
# X1 <- Expr_Cerebral
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt1 <- genefilter(2^X1,ffun)
# Expr_Cerebral <- Expr_Cerebral[filt1,]
# 
Expr_Cerebral <- t(scale(t(Expr_Cerebral), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype

# Control and convalescent VS DF and DHF
# Pheno1
Pheno_Cerebral$DiseaseStatus <- as.factor(Pheno_Cerebral$`disease state:ch1`)
levels(Pheno_Cerebral$DiseaseStatus) <- c("cerebral", "nonCerebral", "nonCerebral") 
table(Pheno_Cerebral$DiseaseStatus)
Pheno_Cerebral$DiseaseStatus <- factor(Pheno_Cerebral$DiseaseStatus, levels = c("nonCerebral", "cerebral"))

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Cerebral)]
#all(rownames(Pheno_Cerebral) == colnames(expr1))

ClassCerebralVsNonCerebral <- Pheno_Cerebral$DiseaseStatus

save(Expr_Cerebral, ClassCerebralVsNonCerebral, file = "./Objs/CerebralExtraValidation.rda")
