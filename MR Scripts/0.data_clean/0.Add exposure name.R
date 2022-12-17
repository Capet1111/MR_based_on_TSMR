#use library
library(rio); library(tidyverse)

Filename = "D:\\00sleep\\GWAS\\FinnGen_gallstone.tsv"
expName = "gallstone"

d1 = import(Filename)
d1 = d1 %>% mutate(exposure = rep(expName,nrow(d1))) %>% 
            select(exposure, everything())
head(d1)

export(d1,Filename)

