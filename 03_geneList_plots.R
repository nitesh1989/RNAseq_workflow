

###clean
rm(list=ls())


###Set Working directory
sample = "BRCA"

# Set path
#my.path = file.path("~/Documents/JHMI-Research/TCGA-Expression-Gene",sample)
my.path = file.path("~/TestRun/TCGA-Expression-Gene",sample)
setwd(my.path)


###libraries
library(geneplotter)

###################################################
### chunk number 2: ###load design and fit objects
load('./objs/lmObjs.ls.rda')
load("expressionMat.rda")

###################################################
###select annotation
annSel <- rownames(expression_mat)

topListOfGenes = toptable(fit=eb.ls,genelist=annSel,adjust.method="fdr",number=100)
save(topListOfGenes,file = "objs/listOfGenes.rda")

write.csv(topListOfGenes,file = "objs/listOfGenes.csv")
