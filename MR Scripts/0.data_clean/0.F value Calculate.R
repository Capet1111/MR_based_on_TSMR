#calculate F value
FileName = "tea45 - ¸±±¾.CSV"
k = 45
expo_type = "continuous"
expo_type = "binary"

#use library
library(TwoSampleMR); library(rio); library(tidyverse)
#import data
IVs = import(FileName)
#calculate R2-------------------------------------------------------------------
if(expo_type == "continuous"){
  IVs$r2 = get_r_from_pn(p = IVs$p,  n = IVs$N)
  IVs = IVs %>% mutate(Fval = ((N-k-1)*r2)/((1-r2)*k))
}
if(expo_type == "binary"){
  IVs$r2 = get_r_from_lor(lor = "beta", af = "eaf",
                          ncase = "ncase",  ncontrol = "ncontrol",
                          model = "logit",  correction = FALSE)
  IVs = IVs %>% mutate(Fval = (((ncase+ncontrol)-k-1)*r2)/((1-r2)*k))
}
#export data--------------------------------------------------------------------
export(IVs, "IVs.csv")

#programme finished
#next


