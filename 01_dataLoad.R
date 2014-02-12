# Nitesh Turaga
# TCGA- Expression-Gene UNC Agilent analysis

library(plyr)

sample = "BRCA"

# Set path
my.path = file.path("~/TestRun/TCGA-Expression-Gene",sample)
setwd(my.path)




# Add unzipping step
tar.file = list.files(my.path)
untar(tar.file,exdir=my.path)

#Read in file_sample_map
sdrf = read.table("FILE_SAMPLE_MAP.txt",header = T)


# Add phenotype to file map
sdrf$sample = substr(x=sdrf$barcode.s.,start =14,stop=15)
sdrf$phenotype = ifelse(sdrf$sample == "01","Cancer","Normal")

# files present in level 1 folder

files = list.files(path = my.path,pattern =".data.txt",recursive=TRUE,full.names=TRUE)
files = data.frame(filename = gsub(".+/","",files),path = files)

# Match the sdrf file with the Level 1 files
matched.list = plyr::join(x= files ,y=sdrf,by="filename")


#Make new Cancer and Normal folder
dir.create("./data")
dir.create(path="data/Cancer");dir.create(path = "data/Normal")

#Populate folders
dataDir = "data"

# Make matched list character vectors
matched.list = data.frame(lapply(matched.list,as.character),stringsAsFactors=FALSE)


# Rename from and to
from=matched.list$path
to = file.path(dataDir,matched.list$phenotype,matched.list$filename)


# Move files
file.rename(from= from,to = to)



