rm(list = ls())

library(GEOquery)

GSE1124 <- getGEO("GSE1124", GSEMatrix = T, AnnotGPL = T)
GSE1124_1 <- GSE1124$`GSE1124-GPL96_series_matrix.txt.gz`
GSE1124_2 <- GSE1124$`GSE1124-GPL97_series_matrix.txt.gz`
save(GSE1124_1, GSE1124_2, file = "./Data/ExtraMalariaDataset2.rda")

load("./Data/ExtraMalariaDataset2.rda")

Expr_Test3_1 <- exprs(GSE1124_1)
Expr_Test3_2 <- exprs(GSE1124_2)

Pheno_Test3_1 <- pData(GSE1124_1)
Pheno_Test3_2 <- pData(GSE1124_2)

FeatData_Test3_1 <- fData(GSE1124_1)
FeatData_Test3_2 <- fData(GSE1124_2)


############################
## Annotation

## 
head(rownames(Expr_Test3_1))
rownames(Expr_Test3_1) <- FeatData_Test3_1$`Gene symbol`
summary(is.na(rownames(Expr_Test3_1)))
#rownames(Expr_Cerebral) <- gsub("-","", rownames(Expr_Cerebral))
#rownames(Expr_Cerebral) <- gsub("_","",rownames(Expr_Cerebral))
sel <- which(apply(Expr_Test3_1, 1, function(x) all(is.finite(x)) ))
Expr_Test3_1 <- Expr_Test3_1[sel, ]
Expr_Test3_1 <- Expr_Test3_1[!is.na(rownames(Expr_Test3_1)),]
dim(Expr_Test3_1)

range(Expr_Test3_1)

Expr_Test3_1 <- log2(Expr_Test3_1)
Expr_Test3_1 <- t(scale(t(Expr_Test3_1), center = TRUE, scale = TRUE))

######3
head(rownames(Expr_Test3_2))
rownames(Expr_Test3_2) <- FeatData_Test3_2$`Gene symbol`
summary(is.na(rownames(Expr_Test3_2)))
#rownames(Expr_Cerebral) <- gsub("-","", rownames(Expr_Cerebral))
#rownames(Expr_Cerebral) <- gsub("_","",rownames(Expr_Cerebral))
sel <- which(apply(Expr_Test3_2, 1, function(x) all(is.finite(x)) ))
Expr_Test3_2 <- Expr_Test3_2[sel, ]
Expr_Test3_2 <- Expr_Test3_2[!is.na(rownames(Expr_Test3_2)),]
dim(Expr_Test3_2)

range(Expr_Test3_2)

Expr_Test3_2 <- log2(Expr_Test3_2)
Expr_Test3_2 <- t(scale(t(Expr_Test3_2), center = TRUE, scale = TRUE))

CommonGns <- intersect(rownames(Expr_Test3_2), rownames(Expr_Test3_1))
CommonGns <- make.names(CommonGns)

rownames(Expr_Test3_1) <- make.names(rownames(Expr_Test3_1))
rownames(Expr_Test3_2) <- make.names(rownames(Expr_Test3_2))

Expr_Test3_1 <- Expr_Test3_1[CommonGns, ]
Expr_Test3_2 <- Expr_Test3_2[CommonGns, ]

Expr_Test3 <- cbind(Expr_Test3_1, Expr_Test3_2)
Expr_Test3 <- normalizeBetweenArrays(Expr_Test3, method = "quantile")
####################################

### Modify the phenotype

# Cerebral vs non-cerebral
# Pheno1
Pheno_Test3_1 <- Pheno_Test3_1[!(Pheno_Test3_1$`diesease status:ch1` == "healthy"), ]
Pheno_Test3_1$DiseaseStatus <- as.factor(Pheno_Test3_1$`diesease status:ch1`)
levels(Pheno_Test3_1$DiseaseStatus) <- c("nonCerebral", "cerebral", "nonCerebral", "nonCerebral") 
table(Pheno_Test3_1$DiseaseStatus)
Pheno_Test3_1$DiseaseStatus <- factor(Pheno_Test3_1$DiseaseStatus, levels = c("nonCerebral", "cerebral"))

Expr_Test3_1 <- Expr_Test3_1[, colnames(Expr_Test3_1) %in% rownames(Pheno_Test3_1)]
all(rownames(Pheno_Test3_1) == colnames(Expr_Test3_1))

ClassCerebralVsNonCerebral1 <- Pheno_Test3_1$DiseaseStatus

########
# Pheno2
Pheno_Test3_2 <- Pheno_Test3_2[!(Pheno_Test3_2$`diesease status:ch1` == "healthy"), ]
Pheno_Test3_2 <- Pheno_Test3_2[!is.na(Pheno_Test3_2$`diesease status:ch1`), ]
Pheno_Test3_2$DiseaseStatus <- as.factor(Pheno_Test3_2$`diesease status:ch1`)
levels(Pheno_Test3_2$DiseaseStatus) <- c("nonCerebral", "cerebral", "nonCerebral", "nonCerebral") 
table(Pheno_Test3_2$DiseaseStatus)
Pheno_Test3_2$DiseaseStatus <- factor(Pheno_Test3_2$DiseaseStatus, levels = c("nonCerebral", "cerebral"))

Expr_Test3_2 <- Expr_Test3_2[, colnames(Expr_Test3_2) %in% rownames(Pheno_Test3_2)]
all(rownames(Pheno_Test3_2) == colnames(Expr_Test3_2))

ClassCerebralVsNonCerebral2 <- Pheno_Test3_2$DiseaseStatus


###################3
# Combine together
ClassCerebralVsNonCerebral3 <- c(ClassCerebralVsNonCerebral1, ClassCerebralVsNonCerebral2)
table(ClassCerebralVsNonCerebral3)
ClassCerebralVsNonCerebral3 <- as.factor(ClassCerebralVsNonCerebral3)
levels(ClassCerebralVsNonCerebral3) <- c("nonCerebral", "cerebral")
Expr_Test3 <- cbind(Expr_Test3_1, Expr_Test3_2)
Expr_Test3 <- normalizeBetweenArrays(Expr_Test3, method = "quantile")
save(Expr_Test3, ClassCerebralVsNonCerebral3, file = "./Objs/CerebralExtraValidation2.rda")


#################################################################
# complicated vs non-complicated
# Pheno1
Pheno_Test3_1 <- Pheno_Test3_1[!(Pheno_Test3_1$`diesease status:ch1` == "healthy"), ]
Pheno_Test3_1$DiseaseStatus <- as.factor(Pheno_Test3_1$`diesease status:ch1`)
levels(Pheno_Test3_1$DiseaseStatus) <- c("unComplicated", "Complicated", "Complicated", "unComplicated") 
table(Pheno_Test3_1$DiseaseStatus)
Pheno_Test3_1$DiseaseStatus <- factor(Pheno_Test3_1$DiseaseStatus, levels = c("unComplicated", "Complicated"))

#Expr_Test3_1 <- Expr_Test3_1[, colnames(Expr_Test3_1) %in% rownames(Pheno_Test3_1)]
#all(rownames(Pheno_Test3_1) == colnames(Expr_Test3_1))

ClassComplicatedVSunComplicated1 <- Pheno_Test3_1$DiseaseStatus

########
# Pheno2
Pheno_Test3_2 <- Pheno_Test3_2[!(Pheno_Test3_2$`diesease status:ch1` == "healthy"), ]
Pheno_Test3_2 <- Pheno_Test3_2[!is.na(Pheno_Test3_2$`diesease status:ch1`), ]
Pheno_Test3_2$DiseaseStatus <- as.factor(Pheno_Test3_2$`diesease status:ch1`)
levels(Pheno_Test3_2$DiseaseStatus) <- c("unComplicated", "Complicated", "Complicated", "unComplicated") 
table(Pheno_Test3_2$DiseaseStatus)
Pheno_Test3_2$DiseaseStatus <- factor(Pheno_Test3_2$DiseaseStatus, levels = c("unComplicated", "Complicated"))

#Expr_Test3_2 <- Expr_Test3_2[, colnames(Expr_Test3_2) %in% rownames(Pheno_Test3_2)]
#all(rownames(Pheno_Test3_2) == colnames(Expr_Test3_2))

ClassComplicatedVSunComplicated2 <- Pheno_Test3_2$DiseaseStatus


###################3
# Combine together
ClassComplicatedVSunComplicated3 <- c(ClassComplicatedVSunComplicated1, ClassComplicatedVSunComplicated2)
table(ClassComplicatedVSunComplicated3)
ClassComplicatedVSunComplicated3 <- as.factor(ClassComplicatedVSunComplicated3)
levels(ClassComplicatedVSunComplicated3) <- c("unComplicated", "Complicated")

#Expr_Test3 <- cbind(Expr_Test3_1, Expr_Test3_2)
#Expr_Test3 <- normalizeBetweenArrays(Expr_Test3, method = "quantile")
save(Expr_Test3, ClassComplicatedVSunComplicated3, file = "./Objs/CompExtraValidation2.rda")

