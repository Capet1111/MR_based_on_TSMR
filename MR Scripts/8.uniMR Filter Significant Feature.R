##perform MR 
rm(list = ls())
MR_data = "data_for_uniMR.csv"
#use library
library(rio); library(tidyverse); library(TwoSampleMR)
#import data
dat <- import(MR_data)  

#performing MR 
mr_res <- mr(dat, method = c("mr_ivw","mr_weighted_median","mr_egger_regression"))

#Sensitivity analysis-----------------------------------------------------------
mr_plt = mr_pleiotropy_test(dat);  mr_het = mr_heterogeneity(dat)
mr_LOO = mr_leaveoneout(dat);      mr_Odd = generate_odds_ratios(mr_res)
mr_OR  = data.frame(Method = mr_Odd$method, OR=mr_Odd$or,
                    OR.95L=mr_Odd$or_lci95, OR.95H = mr_Odd$or_uci95, pvalue=mr_Odd$pval)

#p.adjust
mr_res_ivw = mr_res %>% filter(method=="Inverse variance weighted") %>%
  mutate(FDR = p.adjust(pval,method = "BH")) %>% arrange(pval)

#find significant exposure
sig_ivw = mr_res_ivw %>% filter(pval <= 0.05) %>% arrange(exposure)
expo.sig = sig_ivw$exposure

#significant sensitivity test
sig_res = mr_res %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
sig_plt = mr_plt %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
sig_het = mr_het %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
sig.dat = dat %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
sig_LOO = mr_leaveoneout(sig.dat)

#export significant
export(sig_res, "sig_uniMR_res.csv"); export(sig_plt, "sig_plt.csv")
export(sig_het, "sig_het.csv"); export(sig_LOO, "sig_LOO.csv")
export(sig.dat, "sig.dat.csv"); export(sig_ivw, "sig_ivw.csv")
export(mr_res,"mr_res.csv")
export(sig.dat,"sig_data_for_uniMR.csv")



#PRESSO-------------------------------------------------------------------------
expo = levels(as.factor(sig.dat$exposure))
for (i in expo) {
  message(i,"  is running")
  data = sig.dat %>% filter(exposure == i)
  mr_res_PRESSO = run_mr_presso(data, 
                                NbDistribution = 5000, 
                                SignifThreshold = 0.05)
  PRESSO_OR = mr_res_PRESSO[[1]][["Main MR results"]]
  PRESSO_OR = PRESSO_OR %>% 
    select(Exposure,`MR Analysis`,beta=`Causal Estimate`,se=Sd,`T-stat`,p=`P-value`) %>% 
    mutate(OR = exp(beta),uci = exp(beta + 1.96 * se),lci = exp(beta - 1.96 * se) )
  mr_res_PRESSO[[1]][["Main MR results"]] <-  PRESSO_OR
  i = gsub("\\*","",i);   i = gsub("-","_",i)
  export(PRESSO_OR, sprintf("PRESSO_OR_%s.CSV",i))
  capture.output(mr_res_PRESSO, file = sprintf("mr_res_PRESSO_%s.txt",i))
}




