rm(list = ls())
#perform MR
# formed_MR_data = "data_for_uniMR.csv"
formed_MR_data = "mr_keep.csv"
#MR method
method_list = c("mr_ivw",  "mr_weighted_median",  "mr_egger_regression" )
PRESSO_for_all = T


#Mendelian Randomization--------------------------------------------------------
#use library
library(rio);library(tidyverse);library(TwoSampleMR)
all_ivw_sig = data.frame(); all_sig_res = data.frame(); all_sig_plt = data.frame()
all_sig_het = data.frame(); all_sig_LOO = data.frame(); all_sig_dat = data.frame()

##Mendelian Randomization
dat = import(formed_MR_data)
for (i in levels(as.factor(dat$outcome))) {
  if(dir.exists(i) == F){dir.create(i)}
  data = dat %>% filter(outcome == i)
  #Mendelian Randomization
  mr_res = mr(dat = data, method_list = method_list)
  #FDR p adjust
  mr_res_ivw = mr_res %>% 
    filter(method=="Inverse variance weighted") %>%
    mutate(FDR = p.adjust(pval,method = "BH")) %>% arrange(pval)
  sig_ivw = mr_res_ivw %>% filter(pval <= 0.05) %>% arrange(exposure)
  
  #Sensitivity analysis
  mr_plt = mr_pleiotropy_test(data);  mr_het = mr_heterogeneity(data)
  mr_LOO = mr_leaveoneout(data);  
  mr_Odd = generate_odds_ratios(mr_res)
  mr_OR  = data.frame(
    outcome = mr_Odd$outcome, 
    exposure = mr_Odd$exposure, 
    method = mr_Odd$method, 
    OR = mr_Odd$or,
    OR.95L = mr_Odd$or_lci95, 
    OR.95H = mr_Odd$or_uci95, 
    pval = mr_Odd$pval
  )
  #export data
  export(mr_res_ivw, sprintf("%s/mr_res_ivw.csv",i))
  export(mr_plt, sprintf("%s/mr_plt.csv",i)) 
  export(mr_het, sprintf("%s/mr_het.csv",i))
  export(mr_LOO, sprintf("%s/mr_LOO.csv",i))
  export(mr_res, sprintf("%s/mr_res.csv",i))
  export(mr_Odd, sprintf("%s/mr_Odd.csv",i))
  export(mr_OR,  sprintf("%s/mr_OR.csv", i))
  
  #find significant exposure
  sig_ivw = mr_res_ivw %>% filter(pval <= 0.05) %>% arrange(exposure)
  expo.sig = sig_ivw$exposure
  #significant sensitivity test
  sig_res = mr_res %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
  sig_plt = mr_plt %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
  sig_het = mr_het %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
  sig_LOO = mr_LOO %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
  sig_dat = data %>% filter(exposure %in% expo.sig) %>% arrange(exposure)
  
  all_ivw_sig = rbind(all_ivw_sig,sig_ivw)
  all_sig_res = rbind(all_sig_res,sig_res)
  all_sig_plt = rbind(all_sig_plt,sig_plt)
  all_sig_het = rbind(all_sig_het,sig_het)
  all_sig_LOO = rbind(all_sig_LOO,sig_LOO)
  all_sig_dat = rbind(all_sig_dat,sig_dat)
  
  #export significant
  export(sig_res, sprintf("%s/sig_uniMR_res.csv",i))
  export(sig_plt, sprintf("%s/sig_plt.csv",i))
  export(sig_het, sprintf("%s/sig_het.csv",i))
  export(sig_LOO, sprintf("%s/sig_LOO.csv",i))
  export(sig_dat, sprintf("%s/sig_dat.csv",i))
  export(sig_ivw, sprintf("%s/sig_ivw.csv",i))
  export(sig_dat, sprintf("%s/sig_data_for_MVMR.csv",i))
  
  #rm
  rm(list = c("sig_ivw", "sig_res", "sig_plt", "sig_het", "sig_LOO", "sig_dat"))
  rm(list = c("mr_het", "mr_plt", "mr_het", "mr_LOO", "mr_res", "mr_Odd", "mr_OR", "mr_res_ivw"))
}
view(all_ivw_sig); view(all_sig_res); view(all_sig_plt); view(all_sig_het)

export(all_ivw_sig,"all_ivw_sig.csv")
export(all_sig_res,"all_sig_res.csv")
export(all_sig_plt,"all_sig_plt.csv")
export(all_sig_het,"all_sig_het.csv")
export(all_sig_LOO,"all_sig_LOO.csv")
export(all_sig_dat,"all_sig_dat.csv")

#draw picture-------------------------------------------------------------------
mr_scatter_plot(all_sig_res, all_sig_dat)
mr_funnel_plot(singlesnp_results = mr_singlesnp(all_sig_dat))
mr_leaveoneout_plot(all_sig_LOO)




#PRESSO-------------------------------------------------------------------------
if(PRESSO_for_all == T){
  dat_1 = import(formed_MR_data)
}else{
  dat_1 = import("all_sig_dat.csv")
}
#cycle for outcome
for (j in levels(as.factor(dat_1$outcome))) {
  data_1 = dat_1 %>% filter(outcome == j)
  message(" "); message(j)
  #cycle for exposure
  all_res = data.frame()
  for (i in levels(as.factor(data_1$exposure))) {
    data_2 = data_1 %>% filter(exposure == i)
    mr_res_PRESSO = run_mr_presso(data_2, 
                                  NbDistribution = 5000, 
                                  SignifThreshold = 0.05)
    PRESSO_OR = mr_res_PRESSO[[1]][["Main MR results"]]
    PRESSO_OR = PRESSO_OR %>% 
      select(exposure = Exposure,
             method = `MR Analysis`,
             beta=`Causal Estimate`,
             se=Sd,
             T_stat =`T-stat`,
             pval=`P-value`) %>% 
      mutate(exposure = c(i,i),
             outcome = c(j,j),
             method = c("MR PRESSO Raw","MR PRESSO Corrected"),
             OR = exp(beta),
             OR.95L = exp(beta-1.96*se),
             OR.95H = exp(beta+1.96*se)) %>% 
      select(outcome, exposure, method, OR, OR.95L, OR.95H, pval)
    message("writed into mr_OR.csv")
    export(PRESSO_OR[1,], sprintf("%s/mr_OR.csv",j), append = TRUE)
    message("writed full result")
    capture.output(mr_res_PRESSO, file = sprintf("%s/PRESSO_res_%s_with_%s.txt",j,j,i))
  }
  temp = import(sprintf("%s/mr_OR.csv",j))
  temp = temp %>% arrange(exposure)
  export(temp, sprintf("%s/mr_OR.csv",j), append = F)
}




