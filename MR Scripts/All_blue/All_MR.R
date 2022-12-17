#use library
rm(list = ls())
library(rio); library(tidyverse); library(TwoSampleMR)
#set arguments------------------------------------------------------------------
p.filter = 5e-08
eaf.filter = 0.01
GWAS_Address = "D:\\10.GWAS\\Sleep_trait_both_sex"

#Extract SNPs-------------------------------------------------------------------
GWAS =dir(GWAS_Address); n = length(GWAS) 
all.expo = data.frame()
for (i in GWAS) {
  n=n-1; message(" "); message(paste0(n," GWAS left"))
  #extract SNPs of exposure
  message("importing data ")
  expo = import(paste0(GWAS_Address,"\\",i))      
  message("extracting data ")
  p.sig <- expo %>% filter(pval < p.filter) %>% 
    filter(eaf >= eaf.filter) %>% 
    select(exposure,SNP,effect_allele,other_allele,beta,se,pval)
  #merge SNP into a table
  all.expo = rbind(all.expo,p.sig)
  #get a name for single exposure file
  expoName = unique(p.sig$exposure)
  outName = str_c("expo-",expoName,".csv")
  #export the IVs
  message(paste0(nrow(p.sig), " SNPs selected")); message("writing data ")
  export(p.sig, file=outName)
  message("next "); message(" ")
}
#export data--------------------------------------------------------------------
export(all.expo, "all.expo.csv")
#
#next


#clump exposure SNPs
rm(list = ls())
filename = "all.expo.csv"
#set arguments
clump_r2 = 0.001
clump_kb = 10000

#use library
library(rio); library(tidyverse); library(TwoSampleMR)
outTable = data.frame()
##import the exposure SNP-------------------------------------------------------
expSNP <- read_exposure_data( filename = filename,  sep = ",", 
                              phenotype_col = "exposure", snp_col = "SNP",  beta_col = "beta",
                              se_col = "se", pval_col = "pval", effect_allele_col = "effect_allele",
                              other_allele_col = "other_allele",clump = F  )
#if cycle stopped,try again from this site--------------------------------------
{
  expo = levels(as.factor(expSNP$exposure)); n = length(expo)
  for (i in expo) {
    message(" ");
    message(paste0("------------------- [ ",n," ] left ------------------"))
    message(paste0(i," is running "))
    data = expSNP[expSNP$exposure == i,]
    #clumping the IVs
    clumped <- clump_data(data, clump_r2=clump_r2,clump_kb=clump_kb)
    outTable = rbind(outTable,clumped)
    #write data
    export(outTable,"clumped.csv")
    expSNP = expSNP %>% filter(exposure != i)
    export(expSNP,"left.csv")
    n=n-1;  
    if(n==0){ message(" ");message("------------ run out of exposure ----------")}
  }
}
#export clumped exposure SNPs file----------------------------------------------
export(outTable,paste0("expo_clumped(r2=",clump_r2," kb=",clump_kb,").csv"))





#extract SNPs of outcome
#multiple exposure and multiple outcome
rm(list = ls())
outcome_File_dir = "D:\\10.GWAS\\FinnGene\\R7_K11"
exposureFile = grep("^expo_clumped",dir(),value = T)
p.filter = 5e-8

#===============================================================================
library(rio); library(tidyverse); library(TwoSampleMR)
GWAS = dir(outcome_File_dir)
expSNP = import(exposureFile)

all_lost = data.frame(); all_sig = data.frame(); all_out = data.frame()
n = length(GWAS)

#extract SNP from outcome data--------------------------------------------------
for (i in GWAS) {
  message(" ")
  message("[ ",n ," ] GWAS left")
  message(i, " is running.")
  endpoint <- import(paste0(outcome_File_dir,"\\",i))#import data
  endpoint$pval = as.numeric(endpoint$pval)
  
  #extract SNP from outcome data
  outSNP = endpoint %>% filter(SNP %in% expSNP$SNP)
  lostSNP = expSNP %>% filter(!(SNP %in% outSNP$SNP))
  out_come_sig = outSNP %>% filter(pval <= p.filter)
  outSNP = outSNP %>% filter(pval >= p.filter)
  #merge data
  all_lost = rbind(all_lost, lostSNP)
  all_sig = rbind(all_sig, out_come_sig)
  all_out = rbind(all_out, outSNP)
  n = n-1
}

#export the outcome data--------------------------------------------------------
export(all_lost, "Del_SNP_not_in_outcome.csv")
export(all_sig,"Del_SNP_out_come_sig.CSV")
export(all_out, "outcome_SNPs.csv")

#read exposure SNP into a formal table------------------------------------------
outcome_dat <- read_outcome_data(snps = expSNP$SNP,filename = "outcome_SNPs.csv",
                                 phenotype_col = "exposure", sep = ",", snp_col = "SNP", beta_col = "beta",
                                 effect_allele_col = "effect_allele",other_allele_col = "other_allele",
                                 se_col = "se", pval_col = "pval")
#export the outcome data--------------------------------------------------------
export(outcome_dat,"outcome_SNPs.csv")


#  
###
###next

##harmonize SNPs
expoSNPFile = grep("^expo_clumped",dir(),value = T)
outSNPFile = "outcome_SNPs.csv"

#use library
library(rio);library(tidyverse);library(TwoSampleMR)
#-------------------------------------------------------------------------------
#import data
exposure_dat = import(expoSNPFile); N_SNP = length(unique(exposure_dat$SNP))
outcome_dat = import(outSNPFile)

#===============================================================================
#===============================================================================
#harmonize the two data
dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
#===============================================================================
#===============================================================================

#Scan palindromic SNP
palindromic_del <- dat %>% filter(mr_keep == FALSE)
mr_keep <- dat %>% filter(mr_keep == TRUE)
SNPs = data.frame(SNP = unique(mr_keep$SNP))
#-------------------------------------------------------------------------------
#export the data for MR performing.
export(palindromic_del,"Del_palindromic.csv")
export(mr_keep,"mr_keep.csv")
export(SNPs,"mr_keep_SNP.csv")


##phenoscann
rm(list = ls())
#use library
library(rio); library(tidyverse); library(TwoSampleMR); library(phenoscanner)

dat = import("mr_keep_SNP.csv")
SNPs = dat$SNP
n = ceiling(length(SNPs)/10)
all_scan = data.frame()
n_temp = c(1:10)
#phenoscanner snp query---------------------------------------------------------
for (i in 1:n) {
  message(" "); message("NO: ", i,"/", n, " batch is running. ")
  temp_SNP = na.omit(SNPs[n_temp])
  scan <- phenoscanner(snpquery = temp_SNP, genequery = NULL, regionquery = NULL,
                       catalogue = "GWAS", pvalue = 5e-08, proxies = "None", 
                       r2 = 0.8, build = 37)
  all_scan = rbind(all_scan, scan$results)
  n_temp = n_temp + 10
  if(i == n ){message("run out of SNP")}
}
#export data
export(all_scan, "phenoscan.txt")




#excluding traits
#define excluding traits
traits = c("body mass index", "weight" , "Obesity",
           "Type 2 diabetes",  "Type II diabetes", "Diabetes diagnosed by doctor"
)

#===============================================================================
library(rio); library(tidyverse); library(TwoSampleMR)
traits = gsub(" ","",toupper(traits))
#import data
dat = import("phenoscan.txt")
dat$trait = gsub(" ","",toupper(dat$trait))

#excluding traits
del_scan = dat %>% filter(trait %in% traits)

#export data
export(del_scan, "pleiotropy_SNP.txt")






library(rio);library(tidyverse);library(TwoSampleMR)
#-------------------------------------------------------------------------------
mr_keep = import("mr_keep.csv")
pleiotropy_SNP = import("pleiotropy_SNP.txt")
ple_SNP = unique(pleiotropy_SNP$snp)
#-------------------------------------------------------------------------------
data_for_uniMR = mr_keep %>% filter(!(SNP %in% ple_SNP))
pleiotropy_SNP = merge(pleiotropy_SNP,mr_keep,by.x = "snp",by.y = "SNP",all = F)
export(pleiotropy_SNP,"Del_pleiotropy_SNP.CSV")
export(data_for_uniMR,"data_for_uniMR.csv")





rm(list = ls())
#perform MR
formed_MR_data = "data_for_uniMR.csv"
#MR method
method_list = c("mr_ivw",
                "mr_weighted_median",
                "mr_egger_regression"
)
#use library
library(rio);library(tidyverse);library(TwoSampleMR)
#Mendelian Randomization--------------------------------------------------------
all_ivw_sig = data.frame(); all_sig_res = data.frame(); all_sig_plt = data.frame()
all_sig_het = data.frame(); all_sig_LOO = data.frame(); all_sig_dat = data.frame()
dat = import(formed_MR_data)
endpoint = levels(as.factor(dat$outcome))

##Mendelian Randomization
for (i in endpoint) {
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
  mr_OR  = data.frame(Method = mr_Odd$method, OR=mr_Odd$or,
                      OR.95L=mr_Odd$or_lci95, OR.95H = mr_Odd$or_uci95, pvalue=mr_Odd$pval)
  #export data
  export(mr_res_ivw, paste0(i,"\\","mr_res_ivw.csv"))
  export(mr_plt, paste0(i,"\\","mr_plt.csv")) 
  export(mr_het, paste0(i,"\\","mr_het.csv"))
  export(mr_LOO, paste0(i,"\\","mr_LOO.csv"))
  export(mr_res, paste0(i,"\\","mr_res.csv"))
  export(mr_Odd, paste0(i,"\\","mr_Odd.csv"))
  export(mr_OR,  paste0(i,"\\","mr_OR.csv"))
  
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
  export(sig_res, paste0(i,"\\","sig_uniMR_res.csv") )
  export(sig_plt, paste0(i,"\\","sig_plt.csv") )
  export(sig_het, paste0(i,"\\","sig_het.csv") )
  export(sig_LOO, paste0(i,"\\","sig_LOO.csv") )
  export(sig_dat, paste0(i,"\\","sig_dat.csv") )
  export(sig_ivw, paste0(i,"\\","sig_ivw.csv") )
  export(sig_dat, paste0(i,"\\","sig_data_for_MVMR.csv"))
  
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
expo = levels(as.factor(dat$exposure))

for (i in expo) {
  dat1 = dat %>% filter(exposure == i)
  outc = levels(as.factor(dat1$outcome))
  
  message(" "); message(i)
  
  for (j in outc) {
    data = dat1 %>% filter(outcome == j)
    mr_res_PRESSO = run_mr_presso(data, 
                                  NbDistribution = 5000, 
                                  SignifThreshold = 0.05)
    PRESSO_OR = mr_res_PRESSO[[1]][["Main MR results"]]
    PRESSO_OR = PRESSO_OR %>% 
      select(Exposure,`MR Analysis`,beta=`Causal Estimate`,se=Sd,`T-stat`,pval=`P-value`) %>% 
      mutate(OR = exp(beta),uci = exp(beta + 1.96*se),lci = exp(beta - 1.96*se))
    mr_res_PRESSO[[1]][["Main MR results"]] <-  PRESSO_OR
    
    i = gsub("\\*","",i);   i = gsub("-","_",i)
    j = gsub("\\*","",j);   j = gsub("-","_",j)
    
    export(PRESSO_OR, sprintf("PRESSO_OR_%s_ON_%s.CSV",i,j))
    capture.output(mr_res_PRESSO, file = sprintf("PRESSO_res_%s_ON_%s.txt",i,j))
  }
}




rm(list = ls())
#set file direction
exposure_GWAS_Address = "D:\\10.GWAS\\Sleep_trait_both_sex"
outcome_GWAS_Address = "D:\\10.GWAS\\FinnGene\\R7_K11"
mv_SNP_data = "data_for_uniMR.csv"
#get GWAS data name
expo_GWAS = dir(exposure_GWAS_Address)
ouct_GWAS = dir(outcome_GWAS_Address)

#Use library--------------------------------------------------------------------
library(rio); library(tidyverse); library(TwoSampleMR)
#import data used for uni_variable MR
dat = import(mv_SNP_data); mv_SNP = unique(dat$SNP); rm(dat)

#extract  SNP of all exposure---------------------------------------------------
n = length(expo_GWAS); all.expo = data.frame()
for (i in expo_GWAS) { 
  message(" "); message(paste0(n," GWAS left"))
  data = import(paste0(exposure_GWAS_Address,"\\",i))
  data = data %>% filter(SNP %in% mv_SNP) %>% 
    select(exposure,SNP,effect_allele,other_allele,beta,se,pval,eaf)
  all.expo = rbind(all.expo, data)
  n=n-1
  if(n==0){message(" "); message("run out of GWAS")}
}
export(all.expo, "mv_exposure_data.csv")

##export the formal exposure SNP------------------------------------------------
mv_exposure_dat <- read_exposure_data( filename = "mv_exposure_data.csv",
                                       sep = ",",  phenotype_col = "exposure",   snp_col = "SNP", beta_col = "beta",
                                       se_col = "se",  pval_col = "pval",   effect_allele_col = "effect_allele",
                                       other_allele_col = "other_allele", clump = F  )
export(mv_exposure_dat, "mv_exposure_data.csv")




##extract the outcome SNP-------------------------------------------------------
mv_exposure_dat = import("mv_exposure_data.csv")
all_MVMR_multiple = data.frame()
all_MVMR_lasso = data.frame()


for (j in ouct_GWAS) {
  message(" ");  message(j," is runing")
  endpoint = import(paste0(outcome_GWAS_Address,"\\",j))
  
  j = unique(endpoint$exposure)
  
  ##filter SNPs
  outSNP = endpoint %>% filter(SNP %in% mv_exposure_dat$SNP)
  lostSNP = mv_exposure_dat %>% filter(!(SNP %in% outSNP$SNP))
  export(lostSNP, paste0(j,"\\","SNP_not_in_outcome.csv"))
  export(outSNP,paste0(j,"\\","mv_outcome_data.csv"))
  
  ##export the formal outcome SNP
  outSNPFile = paste0(j,"\\","mv_outcome_data.csv")
  mv_outcome_dat <- read_outcome_data(snps = mv_exposure_dat$SNP, filename = outSNPFile,
                                      sep = ",", snp_col = "SNP",  beta_col = "beta",     se_col = "se", 
                                      phenotype_col = "exposure",   effect_allele_col = "effect_allele",
                                      other_allele_col = "other_allele",    pval_col = "pval"
  )
  export(mv_outcome_dat,paste0(j,"\\","mv_outcome_data.csv"))
  
  #harmonize
  mvdat <- mv_harmonise_data(mv_exposure_dat,mv_outcome_dat)
  
  #multiple_MVMR
  res_multiple <- mv_multiple(mvdat,intercept = F,  instrument_specific = FALSE,
                              pval_threshold = 1); res_multiple 
  all_MVMR_multiple = rbind(all_MVMR_multiple,res_multiple$result)
  
  #export result
  res_multiple_table <- res_multiple$result %>% arrange( pval ) %>% 
    mutate( OR = exp(b), uci = exp(b+1.96*se), lci = exp(b-1.96*se), p = pval )
  export(res_multiple_table, paste0(j,"\\", "MVMR_multiple_result.CSV"))
  
  #lasso_MVMR
  feature = mv_lasso_feature_selection(mvdat)
  if(nrow(feature) != 0 ){
    res_lasso = mv_subset(mvdat,  features = feature,
                          intercept = F, instrument_specific = FALSE, pval_threshold = 1); res_lasso$result
    all_MVMR_lasso = rbind(all_MVMR_lasso,res_multiple$result)
    
    res_lasso_table <- res_lasso$result %>% arrange( pval ) %>%
      mutate( OR = exp(b), uci = exp(b+1.96*se), lci = exp(b-1.96*se), p = pval )
    export(res_lasso_table, paste0(j,"\\", "MVMR_lasso_result.CSV") )
  }
}

export(all_MVMR_multiple,"all_MVMR_multiple.csv")
export(all_MVMR_lasso,"all_MVMR_lasso.csv")
















