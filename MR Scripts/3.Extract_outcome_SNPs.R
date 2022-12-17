#extract SNPs of outcome
endpointFile = "D:\\10.GWAS\\FinnGene\\finngen_R7_C3_BREAST.tsv"
exposureFile = grep("^expo_clumped",dir(),value = T)
p.filter = 5e-8
#use library
library(rio); library(tidyverse); library(TwoSampleMR)
#import endpoint data-----------------------------------------------------------
expSNP = import(exposureFile)
endpoint <- import(endpointFile)
N_SNP = length(unique(expSNP$SNP))
endpoint$pval = as.numeric(endpoint$pval)

#extract SNP from outcome data--------------------------------------------------
outSNP = endpoint %>% filter(SNP %in% expSNP$SNP)
lostSNP = expSNP %>% filter(!(SNP %in% outSNP$SNP))
out_come_sig = outSNP %>% filter(pval <= p.filter)
outSNP = outSNP %>% filter(pval >= p.filter)
#export the outcome data
export(lostSNP, "Del_SNP_not_in_outcome.csv")
export(out_come_sig,"out_come_sig.CSV")
export(outSNP, "outcome_SNPs.csv")

#read exposure SNP into a formal table------------------------------------------
outcome_dat <- read_outcome_data(snps = expSNP$SNP,filename = "outcome_SNPs.csv",
    phenotype_col = "exposure", sep = ",", snp_col = "SNP", beta_col = "beta",
    effect_allele_col = "effect_allele",other_allele_col = "other_allele",
    se_col = "se", pval_col = "pval")
#export the outcome data
export(outcome_dat,"outcome_SNPs.csv")


#  
###
###next