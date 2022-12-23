rm(list = ls())
#set file direction
exposure_GWAS_Address = "D:\\10.GWAS\\Sleep_trait_both_sex"
outcome_GWAS_Address = "D:\\10.GWAS\\FinnGene\\R7_K11"
mv_SNP_data = "data_for_uniMR.csv"
#get GWAS data name
expo_GWAS = dir(exposure_GWAS_Address)
ouct_GWAS = dir(outcome_GWAS_Address)

#Use library--------------------------------------------------------------------
library(rio); library(tidyverse); library(TwoSampleMR)
#import data used for uni_variable MR
dat = import(mv_SNP_data); mv_SNP = unique(dat$SNP); rm(dat)

#extract  SNP of all exposure---------------------------------------------------
n = length(expo_GWAS); all.expo = data.frame()
for (i in expo_GWAS) { 
  message(" "); message("[[ ",n," ]] GWAS left")
  message("Extracting [[", i, " ]]")
  data = import(paste0(exposure_GWAS_Address,"\\",i))
  data = data %>% 
    filter(SNP %in% mv_SNP) %>% 
    select(exposure,SNP,effect_allele,other_allele,beta,se,pval,eaf)
  all.expo = rbind(all.expo, data)
  n=n-1
  if(n==0){message(" "); message("===== run out of GWAS =====")}
}
export(all.expo, "mv_exposure_data.csv")

##export the formal exposure SNP------------------------------------------------
mv_exposure_dat <- read_exposure_data( filename = "mv_exposure_data.csv",
   sep = ",",  phenotype_col = "exposure", snp_col = "SNP", beta_col = "beta",
   se_col = "se",  pval_col = "pval",   effect_allele_col = "effect_allele",
   other_allele_col = "other_allele", clump = F  )
export(mv_exposure_dat, "mv_exposure_data.csv")


##extract the outcome SNP-------------------------------------------------------
mv_exposure_dat = import("mv_exposure_data.csv")
all_MVMR_multiple = data.frame()
all_MVMR_lasso = data.frame()

#cycle for outcome
for (j in ouct_GWAS) {
  message(" ");  message("[[ ", j, " ]] is runing")
  endpoint = import(paste0(outcome_GWAS_Address,"\\",j))
  
  j = unique(endpoint$exposure)
  ##filter SNPs
  outSNP = endpoint %>% filter(SNP %in% mv_exposure_dat$SNP)
  lostSNP = mv_exposure_dat %>% filter(!(SNP %in% outSNP$SNP))
  export(outSNP, sprintf("%s/mv_outcome_data.csv",j))
  
  ##export the formal outcome SNP
  outSNPFile = paste0(j,"\\","mv_outcome_data.csv")
  mv_outcome_dat <- read_outcome_data(snps = mv_exposure_dat$SNP, filename = outSNPFile,
      sep = ",", snp_col = "SNP",  beta_col = "beta",     se_col = "se", 
      phenotype_col = "exposure",   effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",    pval_col = "pval"
  )
  export(mv_outcome_dat, sprintf("%s/mv_outcome_data.csv",j))
  
  #harmonize
  mvdat <- mv_harmonise_data(mv_exposure_dat,mv_outcome_dat)
  
  #multiple_MVMR
  res_multiple <- mv_multiple(mvdat,intercept = F,  instrument_specific = FALSE,
                              pval_threshold = 1); res_multiple 
  all_MVMR_multiple = rbind(all_MVMR_multiple,res_multiple$result)
  
  #export result
  res_multiple_table <- res_multiple$result %>%
    mutate( OR = exp(b), 
            OR.95L = exp(b-1.96*se), 
            OR.95H = exp(b+1.96*se),
            method = rep("MVMR",times = length(exposure))) %>% 
    select(outcome, exposure, method, OR, OR.95L, OR.95H, pval)
  
  export(res_multiple_table, sprintf("%s/MVMR_multiple_result.CSV",j))
  
  export(res_multiple_table, sprintf("%s/mr_OR.csv",j), append = T)
  
  temp = import(sprintf("%s/mr_OR.csv",j))
  temp = temp %>% arrange(exposure)
  export(temp, sprintf("%s/mr_OR.csv",j), append = F)
  
  #lasso_MVMR
  feature = mv_lasso_feature_selection(mvdat)
  if(nrow(feature) != 0 ){
    res_lasso = mv_subset(mvdat,  features = feature,
                          intercept = F, instrument_specific = FALSE, pval_threshold = 1); res_lasso$result
    all_MVMR_lasso = rbind(all_MVMR_lasso,res_multiple$result)
    
    res_lasso_table <- res_lasso$result %>% arrange( pval ) %>%
      mutate( OR = exp(b), 
              OR.95L = exp(b-1.96*se), 
              OR.95H = exp(b+1.96*se),
              method = rep("MVMR",times = length(exposure))) %>% 
      select(outcome, exposure, method, OR, OR.95L, OR.95H, pval)
    export(res_lasso_table, sprintf("%s/MVMR_lasso_result.CSV",j))
  }
}
export(all_MVMR_multiple,"all_MVMR_multiple.csv")
export(all_MVMR_lasso,"all_MVMR_lasso.csv")





