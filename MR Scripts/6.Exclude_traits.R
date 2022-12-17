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


