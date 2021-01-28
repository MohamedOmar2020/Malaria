

rm(list = ls())

library(GEOquery)

GSE35861 <- getGEO("GSE35861", GSEMatrix = T, AnnotGPL = T)
GSE35861 <- GSE35861$GSE35861_series_matrix.txt.gz

save(GSE35861, file = "./Data/ExtraMalariaDataset.rda")

load("./Data/ExtraMalariaDataset.rda")

Expr_Test3 <- exprs(GSE35861)
Pheno_Test3 <- pData(GSE35861)
FeatData_Test3 <- fData(GSE35861)


############################
## Annotation

## Expr_Cerebral
head(rownames(Expr_Test3))
rownames(Expr_Test3) <- FeatData_Test3$
summary(is.na(rownames(Expr_Test3)))
#rownames(Expr_Cerebral) <- gsub("-","", rownames(Expr_Cerebral))
#rownames(Expr_Cerebral) <- gsub("_","",rownames(Expr_Cerebral))
sel <- which(apply(Expr_Test3, 1, function(x) all(is.finite(x)) ))
Expr_Test3 <- Expr_Test3[sel, ]
Expr_Test3 <- Expr_Test3[!is.na(rownames(Expr_Test3)),]
dim(Expr_Test3)

range(Expr_Test3)
Expr_Test3 <- log2(Expr_Test3)
Expr_Test3 <- t(scale(t(Expr_Test3), center = TRUE, scale = TRUE))


####################################

### Modify the phenotype
# Placental malaria +ve vs PM-ve

Pheno_Test3$DiseaseStatus <- as.factor(Pheno_Test3$`malaria status:ch1`)
levels(Pheno_Test3$DiseaseStatus) <- c("unComplicated", "Complicated") 
table(Pheno_Test3$DiseaseStatus)
#Pheno_Test3$DiseaseStatus <- factor(Pheno_Test3$DiseaseStatus, levels = c("unComplicated", "Complicated"))

#expr1 <- expr1[, colnames(expr1) %in% rownames(Pheno_Cerebral)]
#all(rownames(Pheno_Cerebral) == colnames(expr1))


ClassComplicatedVSunComplicated2 <- Pheno_Test3$DiseaseStatus


###################################
####################################


save(Expr_Test3, ClassComplicatedVSunComplicated2, file = "./Objs/ExtraMalaria.rda")






