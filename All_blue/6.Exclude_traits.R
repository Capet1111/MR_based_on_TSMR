#excluding traits
rm(list = ls())
#define excluding traits
traits = c("body mass index", "weight" , "Obesity",
           "Type 2 diabetes",  "Type II diabetes", "Diabetes diagnosed by doctor"
           )

#===============================================================================
library(rio); library(tidyverse); library(TwoSampleMR)
traits = gsub(" ","",toupper(traits))
#import data
dat = import("phenoscan.txt")
dat$trait = gsub(" ","",toupper(dat$trait))
#excluding traits
del_scan = dat %>% filter(trait %in% traits)
#export data
export(del_scan, "pleiotropy_SNP.txt")
#-------------------------------------------------------------------------------
mr_keep = import("mr_keep.csv")
pleiotropy_SNP = import("pleiotropy_SNP.txt")
ple_SNP = unique(pleiotropy_SNP$snp)
#-------------------------------------------------------------------------------
data_for_uniMR = mr_keep %>% filter(!(SNP %in% ple_SNP))
pleiotropy_SNP = merge(pleiotropy_SNP,mr_keep,by.x = "snp",by.y = "SNP",all = F)
export(pleiotropy_SNP,"Del_pleiotropy_SNP.CSV")
export(data_for_uniMR,"data_for_uniMR.csv")
