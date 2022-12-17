##perform MR 
MR_data = "harmonised.data.csv"
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
mr_res_ivw = mr_res %>% filter(method == "Inverse variance weighted") %>% arrange(pval)
ivw$FDR = p.adjust(ivw$pval,method = "BH")

#find significant exposure
sig_ivw = mr_res_ivw %>% filter(pval <= 0.05) %>% arrange(exposure)
#
exposure.sig = sig_ivw$exposure

#significant sensitivity test
sig = mr_res %>% filter(exposure %in% exposure.sig) %>% arrange(exposure)
sig_plt = mr_plt %>% filter(exposure %in% exposure.sig) %>% arrange(exposure)
sig_het = mr_het %>% filter(exposure %in% exposure.sig) %>% arrange(exposure)
sig.dat = dat %>% filter(exposure %in% exposure.sig) %>% arrange(exposure)
sig_LOO = mr_leaveoneout(sig.dat)

#export significant
export(sig, "sig_MR_res.csv"); export(sig_plt, "sig_plt.csv")
export(sig_het, "sig_het.csv"); export(sig_LOO, "sig_LOO.csv")
export(sig.dat, "sig.dat.csv"); export(sig_ivw, "sig_ivw.csv")
export(mr_res,"mr_res.csv")

##
##program
##next
