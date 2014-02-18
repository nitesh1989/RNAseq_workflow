# Nitesh Turaga
# TCGA- Expression-Gene UNC Agilent analysis


# Library laod
library(plyr)
library(reshape)
library(limma)
library(data.table)
library(BiocInstaller)
library(gcrma)
require(gmodels)

sample = "BRCA"

# Set path
#my.path = file.path("~/Documents/JHMI-Research/TCGA-Expression-Gene",sample)
my.path = file.path("~/TestRun/TCGA-Expression-Gene",sample)
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
    # Just to make sure the names are cast on rows and columns
    rownames(f1) = rownames(f1)
    colnames(f1) = colnames(f1)    
    # Return a data matrix
    f1 = data.matrix(f1)
    return(f1)
}

# Build the expression matrix for both types.
normal_dat = build.expression.matrix(normal)
cancer_dat = build.expression.matrix(cancer)

save(normal_dat,cancer_dat,file = "expressionMat.rda",compress =TRUE)
load("expressionMat.rda")

# Phenotype of expression matrix
phenotype = data.frame(filenames = c(colnames(cancer_dat),colnames(normal_dat)),
                       status = c(rep(x="Cancer",length(cancer)),rep(x="Normal",length(normal))))
phenotype$status = as.factor(phenotype$status)

# Design matrix  and make contrasts
forDesign = factor(phenotype$status,levels = unique(phenotype$status))

dMat = model.matrix(~0+forDesign)
colnames(dMat) = gsub('forDesign',"",colnames(dMat))
levels(dMat) = colnames(dMat)
colnames(dMat) = levels(dMat)
cMat = makeContrasts(levels = dMat,CvsN = Cancer-Normal)


# Expression matrix
expression_mat = cbind(cancer_dat,normal_dat)

# Linear model fit
fit = lmFit(expression_mat,design)


###fit model
fit.ls <- lmFit(expression_mat,dMat,method="ls")
fit.ls <- contrasts.fit(fit.ls,cMat)
eb.ls <- eBayes(fit.ls)



# Plot SA of the linear model fit.
#pdf(file = "fit_plot.pdf")
#plotSA(fit,main = "Gene level")
#dev.off()
summary(eb.ls)

save(list=ls(pattern="\\.ls"),file='./objs/lmObjs.ls.rda')



# 
# 
# ##################################################################
# # Previous stuff
# ##################################################################
# #Seperate fits
# fit = ebayes(fit)
# 
# #topTable(fit$p.value[,2], coef="phenotype$statusNormal", adjust="BH")
# 
# #class(fit)
# #names(fit)
# 
# fit$p.corr = p.adjust(fit$p.value[,2],method="bonferroni")
# 
# subset.fit = fit$p.corr[fit$p.corr<=0.0005]
# pdf(file = "heatmap.pdf")
# heatmap(expression_mat[names(subset.fit),])
# dev.off()
# ### If FDR does not give a subset of <5000 genes, use a more stringent correction metho
# summary(subset.fit$p.corr)
# # subset.df = subset.fit
# 
# 

