#############################################################################
#############################################################################
## Mohamed Omar
## 10/4/2019
## Goal: Getting and organizing the data for the next step 
############################################################################
rm(list = ls())

setwd("/Volumes/Macintosh/Research/Projects/Malaria")


library(GEOquery)
library(Biobase)
library(sampling)
library(limma)
library(genefilter)
library(edgeR)
library(caret)

### Getting the data 
# dataset1_2 <- getGEO("GSE1124", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset1 <- dataset1_2$`GSE1124-GPL96_series_matrix.txt.gz`
# dataset2 <- dataset1_2$`GSE1124-GPL97_series_matrix.txt.gz`
# 
# dataset3 <- getGEO("GSE117613", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset3 <- dataset3$GSE117613_series_matrix.txt.gz
# 
# dataset4 <- getGEO("GSE35858", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset4 <- dataset4$GSE35858_series_matrix.txt.gz
# 
# dataset5 <- getGEO("GSE34404", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset5 <- dataset5$GSE34404_series_matrix.txt.gz
# 
# dataset6 <- getGEO("GSE116306", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset6 <- dataset6$GSE116306_series_matrix.txt.gz
# 
# dataset7 <- getGEO("GSE119150", GSEMatrix = TRUE, AnnotGPL = TRUE) # Problem
# dataset7 <- dataset7$GSE119150_series_matrix.txt.gz
# 
# dataset8 <- getGEO("GSE16463", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset8 <- dataset8$GSE16463_series_matrix.txt.gz
# 
# dataset9 <- getGEO("GSE72058", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset9 <- dataset9$GSE72058_series_matrix.txt.gz
# 
# save(dataset1, dataset2, dataset3, dataset4, dataset5, dataset6, dataset7, dataset8, dataset9, file = "./Data/MalariaData.rda")


load("./Data/MalariaData.rda")

##################
## Expression

expr1 <- exprs(dataset1)
expr2 <- exprs(dataset2)
expr3 <- exprs(dataset3)
expr4 <- exprs(dataset4)
expr5 <- exprs(dataset5)
expr6 <- exprs(dataset6)
expr7 <- exprs(dataset7)
expr8 <- exprs(dataset8)
expr9 <- exprs(dataset9)

####################
## Feature data
featData1 <- fData(dataset1)
featData2 <- fData(dataset2)
featData3 <- fData(dataset3)
featData4 <- fData(dataset4)
featData5 <- fData(dataset5)
featData6 <- fData(dataset6)
featData7 <- fData(dataset7)
featData8 <- fData(dataset8)
featData9 <- fData(dataset9)

#############################3
## Phenotype
pheno1 <- pData(dataset1)
pheno2 <- pData(dataset2)
pheno3 <- pData(dataset3)
pheno4 <- pData(dataset4)
pheno5 <- pData(dataset5)
pheno6 <- pData(dataset6)
pheno7 <- pData(dataset7)
pheno8 <- pData(dataset8)
pheno9 <- pData(dataset9)

############################
## Annotation

## Expr1
head(rownames(expr1))
rownames(expr1) <- featData1$`Gene symbol`
summary(is.na(rownames(expr1)))
#rownames(expr1) <- gsub("-","", rownames(expr1))
#rownames(expr1) <- gsub("_","",rownames(expr1))
sel <- which(apply(expr1, 1, function(x) all(is.finite(x)) ))
expr1 <- expr1[sel, ]
expr1 <- expr1[!is.na(rownames(expr1)),]
dim(expr1)

expr1 <- log2(expr1 + 1)
range(expr1)
plot(density(expr1))
boxplot(expr1)
# X1 <- expr1
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt1 <- genefilter(2^X1,ffun)
# expr1 <- expr1[filt1,]
# 
expr1 <- t(scale(t(expr1), center = TRUE, scale = TRUE))

#############################
## Expr2
head(rownames(expr2))
rownames(expr2) <- featData2$`Gene symbol`
expr2 <- expr2[!(rownames(expr2) == ""), ]
#rownames(expr2) <- gsub("-","", rownames(expr2))
# rownames(expr2) <- gsub("_","",rownames(expr2))
sel <- which(apply(expr2, 1, function(x) all(is.finite(x)) ))
expr2 <- expr2[sel, ]
expr2 <- expr2[!is.na(rownames(expr2)),]
dim(expr2)

expr2 <- log2(expr2 + 1)
range(expr2)
plot(density(expr2))
boxplot(expr2)

# X2 <- expr2
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt2 <- genefilter(2^X2,ffun)
# expr2 <- expr2[filt2,]
# dim(expr2)
# 
expr2 <- t(scale(t(expr2), center = TRUE, scale = TRUE))


##############################
## Expr3
head(rownames(expr3))
rownames(expr3) <- featData3$`Gene symbol`
expr3 <- expr3[!(rownames(expr3) == ""), ]
#rownames(expr3) <- gsub("-","", rownames(expr3))
# rownames(expr3) <- gsub("_","",rownames(expr3))
sel <- which(apply(expr3, 1, function(x) all(is.finite(x)) ))
expr3 <- expr3[sel, ]
expr3 <- expr3[!is.na(rownames(expr3)),]
dim(expr3)

range(expr3)
plot(density(expr3))
boxplot(expr3[,1:15])

# X3 <- expr3
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt3 <- genefilter(2^X3,ffun)
# expr3 <- expr3[filt3,]
# dim(expr3)
# 
expr3 <- t(scale(t(expr3), center = TRUE, scale = TRUE))

#############################3
## Expr4
head(rownames(expr4))
rownames(expr4) <- featData4$`Composite Element Database Entry[Gene Symbol]`
expr4 <- expr4[!(rownames(expr4) == ""), ]
#rownames(expr4) <- gsub("-","", rownames(expr4))
# rownames(expr4) <- gsub("_","",rownames(expr4))
sel <- which(apply(expr4, 1, function(x) all(is.finite(x)) ))
expr4 <- expr4[sel, ]
expr4 <- expr4[!is.na(rownames(expr4)),]
dim(expr4)

range(expr4)  ## Z scored
plot(density(expr4))
boxplot(expr4[,1:15])

# X3 <- expr4
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt3 <- genefilter(2^X3,ffun)
# expr4 <- expr4[filt3,]
# dim(expr4)
# 
# expr4 <- t(scale(t(expr4), center = TRUE, scale = TRUE))

#############################3
## Expr5
head(rownames(expr5))
rownames(expr5) <- featData5$`Gene symbol`
expr5 <- expr5[!(rownames(expr5) == ""), ]
#rownames(expr5) <- gsub("-","", rownames(expr5))
# rownames(expr5) <- gsub("_","",rownames(expr5))
sel <- which(apply(expr5, 1, function(x) all(is.finite(x)) ))
expr5 <- expr5[sel, ]
expr5 <- expr5[!is.na(rownames(expr5)),]
dim(expr5)

expr5 <- log2(expr5 + 1)
range(expr5)
plot(density(expr5))
boxplot(expr5[,1:15])

# X3 <- expr5
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt3 <- genefilter(2^X3,ffun)
# expr5 <- expr5[filt3,]
# dim(expr5)
# 
expr5 <- t(scale(t(expr5), center = TRUE, scale = TRUE))

#############################3
## Expr6
head(rownames(expr6))
rownames(expr6) <- featData6$GENE_SYMBOL
expr6 <- expr6[!(rownames(expr6) == ""), ]
#rownames(expr6) <- gsub("-","", rownames(expr6))
# rownames(expr6) <- gsub("_","",rownames(expr6))
sel <- which(apply(expr6, 1, function(x) all(is.finite(x)) ))
expr6 <- expr6[sel, ]
expr6 <- expr6[!is.na(rownames(expr6)),]
dim(expr6)

range(expr6)
plot(density(expr6))
boxplot(expr6[,1:15])

# X3 <- expr6
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt3 <- genefilter(2^X3,ffun)
# expr6 <- expr6[filt3,]
# dim(expr6)
# 
expr6 <- t(scale(t(expr6), center = TRUE, scale = TRUE))

#############################3
## Expr7
head(rownames(expr7))
rownames(expr7) <- featData7$`Gene Symbol`
expr7 <- expr7[!(rownames(expr7) == ""), ]
#rownames(expr7) <- gsub("-","", rownames(expr7))
# rownames(expr7) <- gsub("_","",rownames(expr7))
sel <- which(apply(expr7, 1, function(x) all(is.finite(x)) ))
expr7 <- expr7[sel, ]
expr7 <- expr7[!is.na(rownames(expr7)),]
dim(expr7)

range(expr7)
plot(density(expr7))
boxplot(expr7[,1:10])

# X3 <- expr7
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt3 <- genefilter(2^X3,ffun)
# expr7 <- expr7[filt3,]
# dim(expr7)
# 
expr7 <- t(scale(t(expr7), center = TRUE, scale = TRUE))

#############################3
## Expr8
head(rownames(expr8))
rownames(expr8) <- featData8$`Gene symbol`
expr8 <- expr8[!(rownames(expr8) == ""), ]
#rownames(expr8) <- gsub("-","", rownames(expr8))
# rownames(expr8) <- gsub("_","",rownames(expr8))
sel <- which(apply(expr8, 1, function(x) all(is.finite(x)) ))
expr8 <- expr8[sel, ]
expr8 <- expr8[!is.na(rownames(expr8)),]
dim(expr8)

range(expr8)  ## Z-scored
plot(density(expr8))
boxplot(expr8[,1:10])

# X3 <- expr8
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt3 <- genefilter(2^X3,ffun)
# expr8 <- expr8[filt3,]
# dim(expr8)
# 
# expr8 <- t(scale(t(expr8), center = TRUE, scale = TRUE))

#############################3
## Expr9
head(rownames(expr9))
rownames(expr9) <- featData9$`Gene symbol`
expr9 <- expr9[!(rownames(expr9) == ""), ]
#rownames(expr9) <- gsub("-","", rownames(expr9))
# rownames(expr9) <- gsub("_","",rownames(expr9))
sel <- which(apply(expr9, 1, function(x) all(is.finite(x)) ))
expr9 <- expr9[sel, ]
expr9 <- expr9[!is.na(rownames(expr9)),]
dim(expr9)

expr9 <- log2(expr9 + 1)
range(expr9)  
plot(density(expr9))
boxplot(expr9[,1:10])

# X3 <- expr9
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt3 <- genefilter(2^X3,ffun)
# expr9 <- expr9[filt3,]
# dim(expr9)
# 
expr9 <- t(scale(t(expr9), center = TRUE, scale = TRUE))

################################################################################
###############################################################################
### Modify the phenotype
# Remove controls

# Pheno1
pheno1 <- pheno1[!(pheno1$`diesease status:ch1` == "healthy"), ]
pheno1$DiseaseStatus <- as.factor(pheno1$`diesease status:ch1`)
levels(pheno1$DiseaseStatus) <- c("nonCerebral", "cerebral", "nonCerebral", "nonCerebral") 
table(pheno1$DiseaseStatus)

expr1 <- expr1[, colnames(expr1) %in% rownames(pheno1)]
all(rownames(pheno1) == colnames(expr1))

# Pheno2
pheno2 <- pheno2[!(pheno2$`diesease status:ch1` == "healthy"), ]
pheno2$DiseaseStatus <- as.factor(pheno2$`diesease status:ch1`)
levels(pheno2$DiseaseStatus) <- c("nonCerebral", "cerebral", "nonCerebral", "nonCerebral") 
table(pheno2$DiseaseStatus)

expr2 <- expr2[, colnames(expr2) %in% rownames(pheno2)]
all(rownames(pheno2) == colnames(expr2))

# Pheno3
pheno3 <- pheno3[!(pheno3$`diagnosis:ch1` == "no Plasmodium falciparum infection"), ]
pheno3$DiseaseStatus <- as.factor(pheno3$`diagnosis:ch1`)
levels(pheno3$DiseaseStatus) <- c("cerebral", "nonCerebral") 
table(pheno3$DiseaseStatus)

expr3 <- expr3[, colnames(expr3) %in% rownames(pheno3)]
all(rownames(pheno3) == colnames(expr3))

# Pheno4
pheno4 <- pheno4[!(pheno4$`disease state:ch2` == "Healthy"), ]
pheno4$DiseaseStatus <- as.factor(pheno4$`disease state:ch2`)
levels(pheno4$DiseaseStatus) <- c("nonCerebral", "nonCerebral") 
table(pheno4$DiseaseStatus)

expr4 <- expr4[, colnames(expr4) %in% rownames(pheno4)]
all(rownames(pheno4) == colnames(expr4))

# Pheno5
pheno5 <- pheno5[!(pheno5$source_name_ch1 == "Whole blood, age-matched control"), ]
pheno5$DiseaseStatus <- as.factor(pheno5$source_name_ch1)
levels(pheno5$DiseaseStatus) <- c("nonCerebral") 
table(pheno5$DiseaseStatus)

expr5 <- expr5[, colnames(expr5) %in% rownames(pheno5)]
all(rownames(pheno5) == colnames(expr5))

# Pheno6
pheno6$DiseaseStatus <- as.factor(pheno6$`disease state:ch1`)
levels(pheno6$DiseaseStatus) <- c("cerebral", "nonCerebral", "nonCerebral") 
table(pheno6$DiseaseStatus)

expr6 <- expr6[, colnames(expr6) %in% rownames(pheno6)]
all(rownames(pheno6) == colnames(expr6))

# Pheno7
pheno7 <- pheno7[!(pheno7$`subject status:ch1` == "normal, healthy subject"), ]
pheno7$DiseaseStatus <- as.factor(pheno7$`subject status:ch1`)
levels(pheno7$DiseaseStatus) <- c("nonCerebral") 
table(pheno7$DiseaseStatus)

expr7 <- expr7[, colnames(expr7) %in% rownames(pheno7)]
all(rownames(pheno7) == colnames(expr7))

# Pheno8
pheno8 <- pheno8[pheno8$`disease group:ch1` %in% c("Malaria"), ]
pheno8$DiseaseStatus <- as.factor(pheno8$`disease group:ch1`)
levels(pheno8$DiseaseStatus) <- c("nonCerebral") 
table(pheno8$DiseaseStatus)

expr8 <- expr8[, colnames(expr8) %in% rownames(pheno8)]
all(rownames(pheno8) == colnames(expr8))

# Pheno9
pheno9$DiseaseStatus <- rep("cerebral", nrow(pheno9))
table(pheno9$DiseaseStatus)

########################################################################
########################################################################

allpheno <- list(pheno1, pheno2, pheno3, pheno4, pheno5, pheno6, pheno7, pheno8, pheno9)
names(allpheno) <- c("GSE1124-GPL96", "GSE1124-GPL97", "GSE117613", "GSE35858", "GSE34404", "GSE116306", "GSE119150", "GSE16463", "GSE72058")

allExpr <- list(expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8, expr9)
names(allExpr) <- c("GSE1124-GPL96", "GSE1124-GPL97", "GSE117613", "GSE35858", "GSE34404", "GSE116306", "GSE119150", "GSE16463", "GSE72058")


### Filter phenotype information for the required samples
DiseaseStatus <- mapply(x=allpheno, FUN=function(x) {
  x <- x[,"DiseaseStatus"]
  out <- factor(x, levels=c("nonCerebral", "cerebral"))
  out
})

###################################################################################
##################################################################################

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(allExpr, rownames))


### Filter expression for the required samples
exprsMalaria <- mapply(x=allExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

### Check
all(names(exprsMalaria) == names(DiseaseStatus))

### Check order
all(rownames(allpheno$`GSE1124-GPL96`) == colnames(allExpr$`GSE1124-GPL96`))
all(rownames(allpheno$`GSE1124-GPL97`) == colnames(allExpr$`GSE1124-GPL97`))
all(rownames(allpheno$GSE117613) == colnames(allExpr$GSE117613))
all(rownames(allpheno$GSE35858) == colnames(allExpr$GSE35858))
all(rownames(allpheno$GSE34404) == colnames(allExpr$GSE34404))
all(rownames(allpheno$GSE116306) == colnames(allExpr$GSE116306))
all(rownames(allpheno$GSE119150) == colnames(allExpr$GSE119150))
all(rownames(allpheno$GSE16463) == colnames(allExpr$GSE16463))
all(rownames(allpheno$GSE72058) == colnames(allExpr$GSE72058))


###################################################################
#############
# ## Cross-study validation
# # Leave GSE1124 out
train <- c("GSE117613","GSE35858", "GSE34404", "GSE116306", "GSE119150", "GSE16463", "GSE72058")
test <- c("GSE1124-GPL96", "GSE1124-GPL97")
# 
# ## Training
trainMat <- do.call("cbind", exprsMalaria[train])
trainGroup <- factor(do.call("c", DiseaseStatus[train]))
levels(trainGroup) <- c("nonCerebral", "cerebral")
table(trainGroup)

## Testing
testMat <- do.call("cbind", exprsMalaria[test])
testGroup <- factor(do.call("c", DiseaseStatus[test]))
levels(testGroup) <- c("nonCerebral", "cerebral")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/MalariaDataGood_GSE1124Out_NCvsC.rda")
# 
# #############
# # Leave GSE117613 out
train <- c("GSE1124-GPL96", "GSE1124-GPL97", "GSE35858", "GSE34404", "GSE116306", "GSE119150", "GSE16463", "GSE72058")
test <- c("GSE117613")
# 
# ## Training
trainMat <- do.call("cbind", exprsMalaria[train])
trainGroup <- factor(do.call("c", DiseaseStatus[train]))
levels(trainGroup) <- c("nonCerebral", "cerebral")
table(trainGroup)
length(trainGroup)
# 
# ## Testing
testMat <- exprsMalaria$GSE117613
testGroup <- factor(do.call("c", DiseaseStatus[test]))
levels(testGroup) <- c("nonCerebral", "cerebral")
table(testGroup)
length(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/MalariaDataGood_GSE117613Out_NCvsC.rda")
# 
# #############
# # Leave GSE116306 out
train <- c("GSE1124-GPL96", "GSE1124-GPL97", "GSE117613", "GSE35858", "GSE34404", "GSE119150", "GSE16463", "GSE72058")
test <- c("GSE116306")
# 
# ## Training
trainMat <- do.call("cbind", exprsMalaria[train])
trainGroup <- factor(do.call("c", DiseaseStatus[train]))
levels(trainGroup) <- c("nonCerebral", "cerebral")
table(trainGroup)
length(trainGroup)
# 
# ## Testing
testMat <- exprsMalaria$GSE116306
testGroup <- factor(do.call("c", DiseaseStatus[test]))
levels(testGroup) <- c("nonCerebral", "cerebral")
table(testGroup)
length(testGroup)
 
# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/MalariaDataGood_GSE116306Out_NCvsC.rda")
# 
# #############
# # Leave GSE57813 out
# train <- c("GSE13507","pmid15930337", "GSE32894", "E-MTAB-4321")
# test <- c("GSE57813")
# 
# ## Training
# trainMat <- do.call("cbind", exprsProgression[train])
# trainGroup <- factor(do.call("c", groupProgression_ALL[train]))
# levels(trainGroup) <- c("NoProgression", "Progression")
# table(trainGroup)
# length(trainGroup)
# 
# ## Testing
# testMat <- exprsProgression$GSE57813
# testGroup <- factor(do.call("c", groupProgression_ALL[test]))
# levels(testGroup) <- c("NoProgression", "Progression")
# table(testGroup)
# length(testGroup)
# 
# # Save for cross-study validation
# save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/progressionDataGood_GSE57813Out.rda")
# 
# #############
# # Leave EMTAB-4321 out
# train <- c("GSE13507","pmid15930337", "GSE32894", "GSE57813")
# test <- c("E-MTAB-4321")
# 
# ## Training
# trainMat <- do.call("cbind", exprsProgression[train])
# trainGroup <- factor(do.call("c", groupProgression_ALL[train]))
# levels(trainGroup) <- c("NoProgression", "Progression")
# table(trainGroup)
# length(trainGroup)
# 
# ## Testing
# testMat <- exprsProgression$`E-MTAB-4321`
# testGroup <- factor(do.call("c", groupProgression_ALL[test]))
# levels(testGroup) <- c("NoProgression", "Progression")
# table(testGroup)
# length(testGroup)
# 
# # Save for cross-study validation
# save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/progressionDataGood_EMTABOut.rda")

##################################################################
#####################
## All combined
allMat <- do.call("cbind", exprsMalaria)
allGroup <- unlist(DiseaseStatus)
allStudies <- names(allGroup)

names(allGroup) <- colnames(allMat)
all(colnames(allMat) == names(allGroup))

#############################################################
### WBC count

allpheno$GSE117613$WBC <- as.character(allpheno$GSE117613$`wbc.count:ch1`)
allpheno$GSE117613$WBC <- as.numeric(allpheno$GSE117613$WBC)

allpheno$GSE34404$WBC <- as.character(allpheno$GSE34404$`white blood cells:ch1`)
allpheno$GSE34404$WBC <- as.numeric(allpheno$GSE34404$WBC)

allpheno$GSE116306$WBC <- as.character(allpheno$GSE116306$`leucocytes count (giga/l):ch1`)
allpheno$GSE116306$WBC <- as.numeric(allpheno$GSE116306$WBC)

allpheno$GSE119150$WBC <- as.character(allpheno$GSE119150$`wbc (×10^9/l):ch1`)
allpheno$GSE119150$WBC <- as.numeric(allpheno$GSE119150$WBC)

### Covariates of relevance select complete cases: WBC count
allWBC <- lapply(allpheno, function(x) {
  i <- grep("WBC", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- as.numeric(x[, i  ])
})


allWBC <- unlist(allWBC)


#################################################################################
###############################
### Age
allpheno$GSE117613$AGE <- as.character(allpheno$GSE117613$characteristics_ch1.3)
allpheno$GSE117613$AGE <- gsub("age: ", "", allpheno$GSE117613$AGE)
allpheno$GSE117613$AGE <- as.numeric(allpheno$GSE117613$AGE)


allpheno$GSE34404$AGE <- as.character(allpheno$GSE34404$`age (years):ch1`)
allpheno$GSE34404$AGE <- as.numeric(allpheno$GSE34404$AGE)

allpheno$GSE116306$AGE <- as.character(allpheno$GSE116306$`age:ch1`)
allpheno$GSE116306$AGE <- as.numeric(allpheno$GSE116306$AGE)

allpheno$GSE119150$AGE <- as.character(allpheno$GSE119150$`age (years):ch1`)
allpheno$GSE119150$AGE <- as.numeric(allpheno$GSE119150$AGE)

allpheno$GSE16463$AGE <- as.character(allpheno$GSE16463$`age (years):ch1`)
allpheno$GSE16463$AGE <- as.numeric(allpheno$GSE16463$AGE)

### Covariates of relevance select complete cases: AGE
allAGE <- lapply(allpheno, function(x) {
  i <- grep("^AGE$", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- as.numeric(x[, i  ])
})
allAGE <- unlist(allAGE)

##################################################################################
################################
### Sex
allpheno$GSE117613$GENDER <- as.character(allpheno$GSE117613$`Sex:ch1`)

allpheno$GSE34404$GENDER <- as.character(allpheno$GSE34404$`gender:ch1`)
allpheno$GSE34404$GENDER[allpheno$GSE34404$GENDER == "M"] <- "Male"
allpheno$GSE34404$GENDER[allpheno$GSE34404$GENDER == "F"] <- "Female"

allpheno$GSE116306$GENDER <- as.character(allpheno$GSE116306$`gender:ch1`)
allpheno$GSE116306$GENDER[allpheno$GSE116306$GENDER == "M"] <- "Male"
allpheno$GSE116306$GENDER[allpheno$GSE116306$GENDER == "F"] <- "Female"

allpheno$GSE119150$GENDER <- as.character(allpheno$GSE119150$`gender:ch1`)
allpheno$GSE119150$GENDER[allpheno$GSE119150$GENDER == "male"] <- "Male"

allpheno$GSE16463$GENDER <- as.character(allpheno$GSE16463$`gender:ch1`)

### Covariates of relevance select complete cases: SEX
allGENDER <- lapply(allpheno, function(x) {
  i <- grep("GENDER", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- factor(x[, i  ])
})
allGENDER <- factor(unlist(allGENDER))


#########################################################################

### Assemble in one data.frame and turn numeric
covs <- data.frame(STUDIES=allStudies,
                   WBC=allWBC,
                   GENDER=allGENDER, 
                   AGE=allAGE)

### Prepare vocs for sampling
covs <- sapply(covs , function(x) as.numeric(factor(paste(x))) )


###########################################################################
###SAMPLING

### Balanced stratification
set.seed(333)
trainingOrTesting <- balancedstratification(
  covs[ , , drop=FALSE], strata=1*(allGroup == "cerebral"),
  pik=inclusionprobabilities(1:nrow(covs), nrow(covs) * 0.3),
  comment=TRUE, method=1)

### Show
apply(covs[, -ncol(covs),drop=FALSE], 2, table, allGroup, trainingOrTesting)

### Subset Training
mixTrainMat <- allMat[ , trainingOrTesting == 0]
mixTrainGroup <- allGroup[ trainingOrTesting == 0]
mixTrainStudy <- allStudies[ trainingOrTesting == 0]

### Subset Testing
mixTestMat <- allMat[ , trainingOrTesting == 1]
mixTestGroup <- allGroup[ trainingOrTesting == 1]
mixTestStudy <- allStudies[ trainingOrTesting == 1]

table(mixTrainGroup)
table(mixTestGroup)

###########################################################################
### Save
save(exprsMalaria, 
     mixTrainMat, mixTrainGroup, mixTrainStudy,
     mixTestMat, mixTestGroup, mixTestStudy,
     file="./Objs/MalariaDataGood_NCvsC.rda")

#########################################################################
#########################################################################
########################################################################


