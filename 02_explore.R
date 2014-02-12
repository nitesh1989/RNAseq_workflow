# Nitesh Turaga
# TCGA- Expression-Gene UNC Agilent analysis

library(plyr)
library(reshape)

sample = "BRCA"

# Set path
my.path = file.path("~/TestRun/TCGA-Expression-Gene",sample)
setwd(my.path)

# Set data directory
dataDir ="data"

# list cancer and normal files
cancer = list.files(file.path(dataDir,"Cancer"),full.names=TRUE,recursive=T)
normal = list.files(file.path(dataDir,"Normal"),full.names=TRUE,recursive=T)

# Read in data

f1 = read.delim(cancer[1],header =T,stringsAsFactors =FALSE,sep="\t")
cancer_dat = cbind(f1,sappl)