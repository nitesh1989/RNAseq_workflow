# Nitesh Turaga
# TCGA- Expression-Gene UNC Agilent analysis


# Set path
my.path = "~/TestRun/TCGA-Expression-Gene/"
setwd(my.path)

# Install packages required
packageList = c("EDASeq","edgeR","DESeq","BitSeq","Rsubread","easyRNASeq","goseq","DSS")
biocLite(packageList)

# Load packages
require("DESeq")
require("Biostrings")


