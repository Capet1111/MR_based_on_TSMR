rm(list = ls())
#clump exposure SNPs
filename = "all.expo.csv"
#set arguments
clump_r2 = 0.001
clump_kb = 10000

#use library
library(rio); library(tidyverse); library(TwoSampleMR)
outTable = data.frame()
##import the exposure SNP-------------------------------------------------------
expSNP <- read_exposure_data( filename = filename,  sep = ",", 
            phenotype_col = "exposure", snp_col = "SNP",  beta_col = "beta",
            se_col = "se", pval_col = "pval", effect_allele_col = "effect_allele",
            other_allele_col = "other_allele",clump = F  )
#if cycle stopped,try again from this site--------------------------------------
{
expo = levels(as.factor(expSNP$exposure)); n = length(expo)
for (i in expo) {
  message(" ");
  message(paste0("------------------- [ ",n," ] left ------------------"))
  message(paste0(i," is running "))
  data = expSNP[expSNP$exposure == i,]
  #clumping the IVs
  clumped <- clump_data(data, clump_r2=clump_r2,clump_kb=clump_kb)
  outTable = rbind(outTable,clumped)
  #write data
  export(outTable,"clumped.csv")
  expSNP = expSNP %>% filter(exposure != i)
  export(expSNP,"left.csv")
  n=n-1;  
  if(n==0){ message(" ");message("------------ run out of exposure ----------")}
}
}
#export clumped exposure SNPs file----------------------------------------------
export(outTable,paste0("expo_clumped(r2=",clump_r2," kb=",clump_kb,").csv"))



