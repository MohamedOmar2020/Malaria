
rm(list = ls())

## Check for duplicated samples
library(doppelgangR)

load("./Data/MalariaData.rda")


testesets <- list(GSE1124_GPL96=dataset1,
                  GSE1124_GPL97=dataset2, 
                  GSE117613 = dataset3,
                  GSE35858 = dataset4,
                  GSE34404 = dataset5,
                  GSE116306 = dataset6,
                  GSE119150 = dataset7,
                  GSE16463 = dataset8,
                  GSE72058 = dataset9)

results1 <- doppelgangR(testesets, BPPARAM = MulticoreParam(workers = 14))

save(results1, file = "./Objs/DuplicateCheck.rda")
load("./Objs/DuplicateCheck.rda")


View(results1@summaryresults)

par(mfrow=c(2,2), las=1)
plot(results1)