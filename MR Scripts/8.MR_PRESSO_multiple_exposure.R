
MR_data = "data_for_uniMR.csv"

#use library
library(rio); library(tidyverse); library(MRPRESSO)
rm(list = ls())
#import data 
data = import(MR_data)
expo = levels(as.factor(data$exposure)); n = length(expo)

#cycle every exposure
for (i in expo) {
  dat = data[data$exposure == i,]
  n=n-1;  message("------------------------------------------------------")
  message(paste0(i," is running, [ ",n," ] exposure is waiting."))
  PRESSO = mr_presso( data = dat, 
                      BetaOutcome ="beta.outcome", 
                      BetaExposure = "beta.exposure", 
                      SdOutcome ="se.outcome", 
                      SdExposure = "se.exposure", 
                      OUTLIERtest = T,
                      DISTORTIONtest = T, 
                      NbDistribution = 5000,  
                      SignifThreshold = 0.05  )
  #calculate OR and 95ci
  PRESSO_OR = mr_res_PRESSO[[1]][["Main MR results"]]
  PRESSO_OR = PRESSO_OR %>% 
    select(Exposure, `MR Analysis`, beta = `Causal Estimate`, se = Sd, `T-stat`, p = `P-value`) %>% 
    mutate(OR = exp(beta),uci = exp(beta + 1.96 * se),lci = exp(beta - 1.96 * se), pval = `P-value`)
  mr_res_PRESSO[[1]][["Main MR results"]] <-  PRESSO_OR
  #export result
  export(PRESSO_OR, file = paste0("PRESSO_",i,"_OR.CSV"))
  capture.output(PRESSO,file = paste0("PRESSO_",i,"_result.txt"))
  if(n==0){message(" ");message("Run out of exposure")}
}


