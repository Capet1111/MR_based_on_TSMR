library(rio);library(tidyverse);library(TwoSampleMR)
rm(list = ls())
#-------------------------------------------------------------------------------
mr_keep = import("mr_keep.csv")
pleiotropy_SNP = import("pleiotropy_SNP.txt")
ple_SNP = unique(pleiotropy_SNP$snp)
#-------------------------------------------------------------------------------
data_for_uniMR = mr_keep %>% filter(!(SNP %in% ple_SNP))
pleiotropy_SNP = merge(pleiotropy_SNP,mr_keep,by.x = "snp",by.y = "SNP",all = F)
export(pleiotropy_SNP,"Del_pleiotropy_SNP.CSV")
export(data_for_uniMR,"data_for_uniMR.csv")