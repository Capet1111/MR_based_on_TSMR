library(rio); library(tidyverse)
#Key in your GWAS file name.
GWAS.File = dir()
ann.File = "D:\\10.GWAS\\UKB\\00.×¢ÊÍÎÄ¼þ\\UKB_rsids.tsv"
#import the annotation and GWAS file 
ann = import(ann.File)
n = length(GWAS.File)

for (i in GWAS.File) {
    message(n," GWAS left.")
    message(" Running: ",i)
    outName = str_c("Annoted",i,sep = "_")
    data = import(i)
    #annotate the file
    data = cbind(ann[,1:5],data)
    #export the annotated file
    export(data,outName)
    n = n-1
    rm(data)
}

