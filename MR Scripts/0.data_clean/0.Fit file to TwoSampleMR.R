#Fit file to Two Sample MR package
expFile = "kidney.CSV"

#use library
library(rio); library(tidyverse); library(TwoSampleMR)
#-------------------------------------------------------------------------------
##import the exposure SNP
expSNP <- read_exposure_data( filename = expFile,  sep = ",",
                              phenotype_col = "exposure",
                              snp_col = "SNP",
                              beta_col = "beta",
                              se_col = "se",
                              pval_col = "p",
                              effect_allele_col = "effect_allele",
                              other_allele_col = "other_allele",
                              eaf_col = "eaf",
                              clump = F  )
#-------------------------------------------------------------------------------
export(expSNP,"expo_SNPs.csv")

