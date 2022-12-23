##harmonize SNPs
rm(list = ls())
expoSNPFile = grep("^expo_clumped",dir(),value = T)
outSNPFile = "outcome_SNPs.csv"

#use library
library(rio);library(tidyverse);library(TwoSampleMR)
#import data--------------------------------------------------------------------
exposure_dat = import(expoSNPFile); N_SNP = length(unique(exposure_dat$SNP))
outcome_dat = import(outSNPFile)

#===============================================================================
#===============================================================================
#harmonize the two data
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
#===============================================================================
#===============================================================================

#Scan palindromic SNP
palindromic_del <- dat %>% filter(mr_keep == FALSE)
mr_keep <- dat %>% filter(mr_keep == TRUE)
SNPs = data.frame(SNP = unique(mr_keep$SNP))
#export the data for MR performing----------------------------------------------
export(palindromic_del,"Del_palindromic.csv")
export(mr_keep,"mr_keep.csv")
export(SNPs,"mr_keep_SNP.csv")



