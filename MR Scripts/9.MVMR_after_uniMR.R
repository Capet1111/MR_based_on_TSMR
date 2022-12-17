#set file direction
GWAS_Address = "D:\\10.GWAS\\Sleep_trait"
endpointFile = "D:\\10.GWAS\\FinnGene\\finngen_R7_C3_BREAST.tsv"
mv_SNP_data = "data_for_uniMR.csv"

#Use library--------------------------------------------------------------------
library(rio); library(tidyverse); library(TwoSampleMR)
#import data used for uni_variable MR
dat = import(mv_SNP_data); mv_SNP = unique(dat$SNP); rm(dat)
#get exposure GWAS data name
GWAS = dir(GWAS_Address)

#extract  SNP of all exposure---------------------------------------------------
n = length(GWAS); all.expo = data.frame()
for (i in GWAS) { 
  message(" "); message(paste0(n," GWAS left"))
  data = import(paste0(GWAS_Address,"\\",i))
  data = data %>% filter(SNP %in% mv_SNP) %>% 
    select(exposure,SNP,effect_allele,other_allele,beta,se,pval,eaf)
  all.expo = rbind(all.expo, data)
  n=n-1
  if(n==0){message(" "); message("run out of GWAS")}
}
export(all.expo, "mv_exposure_data.csv")

##export the formal exposure SNP------------------------------------------------
mv_exposure_dat <- read_exposure_data( filename = "mv_exposure_data.csv",
                              sep = ",",  phenotype_col = "exposure",
                              snp_col = "SNP", beta_col = "beta",
                              se_col = "se",  pval_col = "pval",
                              effect_allele_col = "effect_allele",
                              other_allele_col = "other_allele", clump = F  )
export(mv_exposure_dat, "mv_exposure_data.csv")

##extract the outcome SNP-------------------------------------------------------
endpoint = import(endpointFile)
##filter SNPs
outSNP = endpoint %>% filter(SNP %in% mv_exposure_dat$SNP)
lostSNP = mv_exposure_dat %>% filter(!(SNP %in% outSNP$SNP))
export(lostSNP, "SNP_not_in_outcome.csv");export(outSNP,"mv_outcome_data.csv")

##export the formal outcome SNP-------------------------------------------------
outSNPFile = "mv_outcome_data.csv"
mv_outcome_dat <- read_outcome_data(snps = mv_exposure_dat$SNP, filename = outSNPFile,
                                    sep = ",", snp_col = "SNP",beta_col = "beta", 
                                    se_col = "se", phenotype_col = "exposure",
                                    effect_allele_col = "effect_allele",
                                    other_allele_col = "other_allele",
                                    pval_col = "pval"
)
export(mv_outcome_dat,"mv_outcome_data.csv")



#harmonize and MR---------------------------------------------------------------
rm(list = ls())
mv_exposure_dat <- import("mv_exposure_data.csv")
mv_outcome_dat <- import("mv_outcome_data.csv")
#harmonize
mvdat <- mv_harmonise_data(mv_exposure_dat,mv_outcome_dat)

#multiple_MVMR------------------------------------------------------------------
res_multiple <- mv_multiple(mvdat,intercept = F,  instrument_specific = FALSE,
                            pval_threshold = 1); res_multiple 

#lasso_MVMR---------------------------------------------------------------------
res_lasso = mv_subset(mvdat,  features = mv_lasso_feature_selection(mvdat),
     intercept = F, instrument_specific = FALSE, pval_threshold = 1); res_lasso$result

#other methods------------------------------------------------------------------
res_residual = mv_residual(mvdat,intercept = FALSE, instrument_specific = FALSE,
                           pval_threshold = 1); res_residual$result

res_ivw = mv_ivw(mvdat,pval_threshold = 1); res_ivw$result

res_basic = mv_basic(mvdat,pval_threshold = 1); res_basic$result



#export result------------------------------------------------------------------
res_multiple_table <- res_multiple$result %>% arrange( pval ) %>% 
  mutate( OR = exp(b), uci = exp(b+1.96*se), lci = exp(b-1.96*se), p = pval )
export(res_multiple_table, "MVMR_multiple_result.CSV")

res_lasso_table <- res_lasso$result %>% arrange( pval ) %>%
  mutate( OR = exp(b), uci = exp(b+1.96*se), lci = exp(b-1.96*se), p = pval )
export(res_lasso_table, "MVMR_lasso_result.CSV")

view(res_multiple_table)
view(res_lasso_table)
