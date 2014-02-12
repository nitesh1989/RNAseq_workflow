# Nitesh Turaga
# TCGA- Expression-Gene UNC Agilent analysis

library(plyr)
library(reshape)

sample = "BRCA"

# Set path
my.path = file.path("~/TestRun/TCGA-Expression-Gene",sample)
setwd(my.path)

dataDir = file.path(sample,"data")

files = list.files(dataDir,full.names=TRUE,recursive=T)
sdad