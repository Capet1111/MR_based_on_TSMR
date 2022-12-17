#filter SNP by p value and eaf

#key in the exposure GWAS file name
expoFile = "Annoted_BMI.21001_irnt.gwas.imputed_v3.both_sexes.tsv"
p.filter = 5e-08

#use library
library(rio); library(tidyverse); library(TwoSampleMR)
#extract SNPs of exposure
expo = import(expoFile)
#filter p value
p.sig <- expo %>% filter(pval <=p.filter) %>% filter(minor_AF >= 0.05)
rs = grep("^rs.",p.sig$rsids,value = T); p.sig = p.sig %>% filter(rsids %in% rs)

#export the IVs
export(p.sig, "exposure_psig.csv")

#programme finished
#next