library(rio); library(tidyverse)
#Key in your GWAS file name.
GWAS.File = "I10.gwas.imputed_v3.both_sexes.tsv"

#===============================================================================
ann.File = "D:\\10.GWAS\\UKB\\00.×¢ÊÍÎÄ¼þ\\UKB_rsids.tsv"
outName = str_c("Anno",GWAS.File,sep = "_")
#import the annotation and GWAS file 
ann = import(ann.File);data = import(GWAS.File)
#annotate the file
data = cbind(ann[,1:5],data)
#export the annotated file
export(data,outName)
