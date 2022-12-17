#use library
library(rio)
library(tidyverse)
library(TwoSampleMR)

expoSNPFile = c("expo-Getting up in morning.csv",
                "expo-Morningevening person (chronotype).csv",
                "expo-Nap during day.csv",
                "expo-narcolepsy.csv",
                "expo-Sleep duration.csv",
                "expo-Sleeplessness  insomnia.csv",
                "expo-Snoring.csv")



exposure_dat = mv_extract_exposures_local(
                    filenames_exposure = expoSNPFile,
                    sep = ",",
                    phenotype_col = "exposure",
                    snp_col = "SNP",
                    beta_col = "beta",
                    se_col = "se",
                    eaf_col = "eaf",
                    effect_allele_col = "effect_allele",
                    other_allele_col = "other_allele",
                    pval_col = "p",
                    pval_threshold = 5e-08,
                    clump_r2 = 0.001,
                    clump_kb = 10000,
                    harmonise_strictness = 2
)

outSNPFile = "outcome_SNPs.csv"

#read exposure SNP into a formal table
outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,
                                 filename = outSNPFile,
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "beta",
                                 se_col = "se",
                                 effect_allele_col = "effect_allele",
                                 other_allele_col = "other_allele",
                                 pval_col = "p" )


mvdat <- mv_harmonise_data(exposure_dat,  outcome_dat) # 对数据进行harmonisation

res <- mv_multiple(mvdat) #计算结果，大家也可以试试mv_residual()函数
res
