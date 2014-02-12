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
for (i in 2:length(cancer)){
    f2 = read.delim(cancer[i],header =T,stringsAsFactors =FALSE,sep = "\t")
    f1 = cbind(f1,f2[,2])
}

cancer_dat = f1;
rm(f1)

concat_array = function(files){
    f = lapply(cancer,read.delim,stringsAsFactors=F,header=T,sep="\t")
    x = as.data.frame(f)
    reshape::merge_all(dfs=f,by = "log2.lowess.normalized..cy5.cy3..collapsed.by.gene.symbol")
    y = do.call(what="cbind",f)
#     f1 = read.delim(cancer[1],header =T,stringsAsFactors =FALSE,sep="\t")
#     for (i in 2:length(files)){
#         f2 = read.delim(files[i],header =T,stringsAsFactors =FALSE,sep = "\t")
#         f1 = cbind(f1,f2[,2])
#     }
#     return(f1)    
}

