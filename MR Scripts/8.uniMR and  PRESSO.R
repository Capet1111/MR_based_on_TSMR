#perform MR
formed_MR_data = "data_for_uniMR.csv"
#use library
rm(list = ls())
library(rio);library(tidyverse);library(TwoSampleMR)
#Mendelian Randomization--------------------------------------------------------
dat = import(formed_MR_data)
mr_res = mr(dat,
            method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))

#FDR p adjust-------------------------------------------------------------------
mr_res_ivw = mr_res %>% filter(method=="Inverse variance weighted") %>%
             mutate(FDR = p.adjust(pval,method = "BH")) %>% arrange(pval)

sig_ivw = mr_res_ivw %>% filter(pval <= 0.05) %>% arrange(exposure)

#Sensitivity analysis-----------------------------------------------------------
mr_plt = mr_pleiotropy_test(dat);  mr_het = mr_heterogeneity(dat)
mr_LOO = mr_leaveoneout(dat);      mr_Odd = generate_odds_ratios(mr_res)
mr_OR  = data.frame(Method = mr_Odd$method, OR=mr_Odd$or,
    OR.95L=mr_Odd$or_lci95, OR.95H = mr_Odd$or_uci95, pvalue=mr_Odd$pval)

#export data--------------------------------------------------------------------
export(mr_res_ivw, "mr_res_ivw.csv")
export(mr_plt, "mr_plt.csv"); export(mr_het, "mr_het.csv")
export(mr_LOO, "mr_LOO.csv"); export(mr_res, "mr_res.csv")
export(mr_Odd, "mr_Odd.csv"); export(mr_OR,"mr_OR.csv")

#draw picture-------------------------------------------------------------------
mr_scatter_plot(mr_res, dat)
mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
mr_leaveoneout_plot(mr_LOO)

#PRESSO-------------------------------------------------------------------------
dat1 = dat %>% filter(exposure == "Nap during day")
mr_res_PRESSO = run_mr_presso(dat1, NbDistribution = 5000, SignifThreshold = 0.05);mr_res_PRESSO

PRESSO_OR = mr_res_PRESSO[[1]][["Main MR results"]]
PRESSO_OR = PRESSO_OR %>% 
  select(Exposure,`MR Analysis`,beta=`Causal Estimate`,se=Sd,`T-stat`,p=`P-value`) %>% 
  mutate(OR = exp(beta),uci = exp(beta + 1.96 * se),lci = exp(beta - 1.96 * se) )
mr_res_PRESSO[[1]][["Main MR results"]] <-  PRESSO_OR
export(PRESSO_OR, "PRESSO_OR.CSV")
capture.output(mr_res_PRESSO, file = "mr_res_PRESSO.txt")










#Adjust SNP
adj_dat = dat[-c(22),]


adj_mr_res = mr(adj_dat,method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression"))
adj_mr_plt = mr_pleiotropy_test(adj_dat)
adj_mr_het = mr_heterogeneity(adj_dat)
adj_mr_LOO = mr_leaveoneout(adj_dat)
adj_mr_Odd = generate_odds_ratios(adj_mr_res)

adj_PRESSO = run_mr_presso(adj_dat, NbDistribution = 5000, SignifThreshold = 0.05)
capture.output(adj_PRESSO, "adj_PRESSO.txt")

export(adj_mr_plt, "adj_mr_plt.csv"); export(adj_mr_het, "adj_mr_het.csv")
export(adj_mr_LOO, "adj_mr_LOO.csv"); export(adj_mr_res,"adj_mr_res.csv")
export(adj_mr_Odd,"adj_mr_Odd.csv")

mr_scatter_plot(mr_results = mr(adj_dat,method_list = c("mr_ivw","mr_weighted_median","mr_egger_regression")),adj_dat)
mr_funnel_plot(singlesnp_results = mr_singlesnp(adj_dat))
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(adj_dat))


adj_OR = data.frame()
adj_OR = cbind(Method = adj_mr_Odd$method,OR=adj_mr_Odd$or,OR.95L=adj_mr_Odd$or_lci95,OR.95H=adj_mr_Odd$or_uci95,pvalue=adj_mr_Odd$pval)
adj_OR = as.data.frame(adj_OR)
export(adj_OR,"adj_OR.csv")






