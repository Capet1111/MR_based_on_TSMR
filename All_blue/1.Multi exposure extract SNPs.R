rm(list = ls())
#set arguments------------------------------------------------------------------
p.filter = 5e-08
eaf.filter = 0.01
GWAS_Address = "D:\\10.GWAS\\PGC_data"

#use library
library(rio); library(tidyverse); library(TwoSampleMR)
#Extract SNPs-------------------------------------------------------------------
GWAS =dir(GWAS_Address); n = length(GWAS) 
all.expo = data.frame()
for (i in GWAS) {
  message(" "); message("[ ", n, " ] GWAS left"); message("[ ",i, " ] is running")
  #extract SNPs of exposure
  expo = import(paste0(GWAS_Address,"\\",i))      
  message("extracting data ")
  if("eaf" %in% colnames(expo))
  {
    p.sig <- expo %>% filter(pval < p.filter) %>% 
      filter(eaf >= eaf.filter) %>% 
      select(exposure,SNP,effect_allele,other_allele,beta,se,pval)
  }else{
    p.sig <- expo %>% filter(pval < p.filter) %>% 
      select(exposure,SNP,effect_allele,other_allele,beta,se,pval)
  }
  #Loosen P filter if no SNP was available.
  n_10 = 1
  while( nrow(p.sig) < 3){
    if(nrow(p.sig) < 3){
      message("Less than 3 SNP was available under p_filter = ", p.filter*n_10)    
      }
    n_10 = n_10*10
    if("eaf" %in% colnames(expo))
    {
      p.sig <- expo %>% 
        filter(pval < p.filter*n_10) %>% 
        filter(eaf >= eaf.filter) %>% 
        select(exposure,SNP,effect_allele,other_allele,beta,se,pval)
    }else{
      p.sig <- expo %>% 
        filter(pval < p.filter*n_10) %>% 
        select(exposure,SNP,effect_allele,other_allele,beta,se,pval)
    }
  }
  #merge SNP into a table
  all.expo = rbind(all.expo,p.sig)
  #get a name for single exposure file
  outName = sprintf("expo-%s_%s.csv",unique(expo$exposure), p.filter*n_10)
  #export the IVs
  export(p.sig, file = outName)
  message(nrow(p.sig), " SNPs selected under p_filter = ", p.filter*n_10)
  message(" ");   n=n-1
  if(n == 0){message(""); message("Run out of GWAS data. ")}
}
#export data--------------------------------------------------------------------
export(all.expo, "all.expo.csv")
#
#next



