#use library
rm(list = ls());library(rio); library(tidyverse)
#import data--------------------------------------------------------------------
Filename = "D:\\10.GWAS\\FinnGene\\NEW\\finngen_R7_T2D_WIDE.tsv"
data = import(Filename);head(data)
#mutate exposure name
exposure_name = "finngen_R7_T2D_WIDE"

#rename variable name-----------------------------------------------------------
{
  SNP = c("SNP","rsids","SNP","MarkerName")
  effect_allele = c("effect_allele","alt","EA")
  other_allele = c("other_allele","ref","OA")
  beta = c("beta","Beta","b","Effect")
  se = c("sebeta","se","SE","StdErr")
  eaf = c("eaf","af_alt","minor_AF","MinFreq")
  pval = c("pval","p","P-value")
}
#data clean---------------------------------------------------------------------
data = data %>% 
  mutate(exposure=rep(exposure_name,nrow(data))) %>% 
  select(exposure, SNP=any_of(SNP), effect_allele=any_of(effect_allele),
         other_allele=any_of(other_allele), beta=any_of(beta), se=any_of(se), 
          pval=any_of(pval), everything())
if(is.numeric(data$beta) == F){data$beta = as.numeric(data$beta)}
if(is.numeric(data$se) == F){data$se = as.numeric(data$se)}
if(is.numeric(data$pval) == F){data$pval = as.numeric(data$pval)}

if( is.na( intersect(colnames(data),eaf) ) == F ){
data = data %>% select(exposure, SNP, effect_allele, other_allele, beta, se, pval, 
                       eaf=any_of(eaf),everything());head(data)
if(is.numeric(data$eaf)==F){data$eaf = as.numeric(data$eaf)}
}
head(data)

#export data--------------------------------------------------------------------
export(data,Filename); message("program done")



