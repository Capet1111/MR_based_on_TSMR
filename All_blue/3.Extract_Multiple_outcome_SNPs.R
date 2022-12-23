#extract SNPs of outcome
rm(list = ls())
#multiple exposure and multiple outcome
outcome_File_dir = "D:\\10.GWAS\\Sleep_trait_both_sex"
exposureFile = grep("^expo_clumped",dir(),value = T)
p.filter = 5e-8

#===============================================================================
library(rio); library(tidyverse); library(TwoSampleMR)
GWAS = dir(outcome_File_dir)
expSNP = import(exposureFile)

all_lost = data.frame(); all_sig = data.frame(); all_out = data.frame()
n = length(GWAS)

#extract SNP from outcome data--------------------------------------------------
for (i in GWAS) {
  message(" ")
  message("[ ",n ," ] GWAS left")
  message(i, " is running.")
  endpoint <- import(paste0(outcome_File_dir,"\\",i))#import data
  endpoint$pval = as.numeric(endpoint$pval)
  
  #extract SNP from outcome data
  outSNP = endpoint %>% 
    filter(SNP %in% expSNP$SNP) %>% 
    select(exposure,SNP,beta,se,pval,effect_allele,other_allele)

  lostSNP = expSNP %>% filter(!(SNP %in% outSNP$SNP))
  out_come_sig = outSNP %>% filter(pval <= p.filter)
  outSNP = outSNP %>% filter(pval >= p.filter)
  #merge data
  all_lost = rbind(all_lost, lostSNP)
  all_sig = rbind(all_sig, out_come_sig)
  all_out = rbind(all_out, outSNP)
  n = n-1
}

#export the outcome data--------------------------------------------------------
export(all_lost, "Del_SNP_not_in_outcome.csv")
export(all_sig,"Del_SNP_out_come_sig.CSV")
export(all_out, "outcome_SNPs.csv")

#read exposure SNP into a formal table------------------------------------------
outcome_dat <- read_outcome_data(snps = expSNP$SNP,filename = "outcome_SNPs.csv",
   phenotype_col = "exposure", sep = ",", snp_col = "SNP", beta_col = "beta",
   effect_allele_col = "effect_allele",other_allele_col = "other_allele",
   se_col = "se", pval_col = "pval")
#export the outcome data--------------------------------------------------------
export(outcome_dat,"outcome_SNPs.csv")


#  
###
###next