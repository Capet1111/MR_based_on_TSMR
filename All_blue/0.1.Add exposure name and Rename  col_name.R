#use library
rm(list = ls());library(rio); library(tidyverse)
#mutate exposure name
Filename = "PGC_data/TS_Oct2018.tsv"
exposure_name = "PGC_TS"
#import data--------------------------------------------------------------------
data = import(Filename);head(data)

#rename variable name-----------------------------------------------------------
SNP = c("SNP","rsids","SNP","MarkerName","ID")
effect_allele = c("effect_allele","alt","EA","ALT","A1")
other_allele = c("other_allele","ref","OA","REF","A2")
beta = c("beta","BETA","Beta","b","Effect","LogOR")
se = c("sebeta","se","SE","StdErr","StdErrLogOR")
eaf = c("eaf","EAF","af_alt","minor_AF","MinFreq","Freq")
pval = c("pval","PVAL","p","P","P-value")
OR = c("OR")
chr = c("CHR")
pos = c("POS")


#data clean---------------------------------------------------------------------
#convert OR into beta
if(!is_empty(intersect(ls(data),OR))){
  data = data %>% 
    select(OR = any_of(OR),everything()) %>% 
    mutate(beta = log(OR))
}
#mutate a "chr:pos" SNP
if((!is_empty(intersect(ls(data),chr)))&(!is_empty(intersect(ls(data),pos)))){
  data = data %>% 
    select(chr = any_of(chr),pos = any_of(pos),everything()) %>% 
    mutate('chr:pos' = str_c("chr",chr,":",pos))
}
#rename colnames
data = data %>% 
  mutate(exposure=rep(exposure_name,nrow(data))) %>% 
  select(exposure, SNP=any_of(SNP), effect_allele=any_of(effect_allele),
         other_allele=any_of(other_allele), beta=any_of(beta), se=any_of(se), 
          pval=any_of(pval), everything()
         )
#make sure numeric status
if(!is.numeric(data$beta)){data$beta = as.numeric(data$beta)}
if(!is.numeric(data$se)){data$se = as.numeric(data$se)}
if(!is.numeric(data$pval)){data$pval = as.numeric(data$pval)}
#rename EAF if there is.
if(!is_empty(intersect(colnames(data),eaf))){
  data = data %>% 
    select(exposure, SNP, effect_allele, other_allele, 
           beta, se, pval, eaf=any_of(eaf),everything())
  if(!is.numeric(data$eaf)){
    data$eaf = as.numeric(data$eaf)
    }
}
head(data)

#export data--------------------------------------------------------------------
export(data,Filename); message("program done")



