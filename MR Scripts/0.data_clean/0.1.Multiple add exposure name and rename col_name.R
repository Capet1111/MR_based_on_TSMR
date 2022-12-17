
library(rio); library(tidyverse)

anno = import("E:\\00.GDMs\\Serum ~ Gallstone\\GWAS\\metaboliteMap.txt")

GWAS =dir(); n = length(GWAS) 


for (i in GWAS) {
  data = import(i)
  
  id = gsub(".xenobiotics.txt","",i)
  exposure_name = anno[anno$metabolonID == id,2]
  
  data = data %>% 
    mutate(exposure=rep(exposure_name,nrow(data))) %>% 
    select(exposure, 
           SNP = MarkerName, 
           effect_allele = Allele1,
           other_allele = Allele2, 
           beta = Effect, 
           se = StdErr, 
           eaf = Freq1,
           pval = 'P-value', 
           everything()
           )
  Filename = str_c("Formated_",i)
  export(data,Filename)
}


