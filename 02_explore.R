# Nitesh Turaga
# TCGA- Expression-Gene UNC Agilent analysis

library(plyr)
library(reshape)
library(data.table)


sample = "BRCA"

# Set path
my.path = file.path("~/Documents/JHMI-Research/TCGA-Expression-Gene",sample)
setwd(my.path)

# Set data directory
dataDir ="data"

# list cancer and normal files

cancer = list.files(file.path(dataDir,"Cancer"),full.names=TRUE,recursive=T)
normal = list.files(file.path(dataDir,"Normal"),full.names=TRUE,recursive=T)
f1 = read.delim(cancer[1],header =T,stringsAsFactors =FALSE,sep="\t")

# Read in data and build the expression matrix,
# where rownames are genes and columns are samples,
# each value is log loess normalized.

build.expression.matrix = function(list.of.files) {
    f1 = read.delim(list.of.files[1],header =T,stringsAsFactors =FALSE,sep="\t")
    rownames(f1) = f1[,1]
    f1 = f1[-1,]
    f1 = as.data.frame(f1)
    f1$Hybridization.REF = NULL
    for (i in 2:length(list.of.files)){
        f2 = read.delim(list.of.files[i],header =T,stringsAsFactors =FALSE,sep = "\t")
        f2 = f2[-1,]
        f2$Hybridization.REF = NULL    
        f1 = cbind(f1,f2)
    }
    return(f1)
}

# Build the expression matrix for both types.
normal_dat = build.expression.matrix(normal)
cancer_dat = build.expression.matrix(cancer)

rownames(cancer_dat) = rownames(cancer_dat)
colnames(cancer_dat) = colnames(cancer_dat)
cancer_data = data.matrix(cancer_dat) 

