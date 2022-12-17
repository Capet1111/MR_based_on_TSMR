
Extract_SNPs <- function( p.filter, eaf.filter, GWAS_Address ){
  #use library
  library(rio); library(tidyverse); library(TwoSampleMR)
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
  export(all.expo, "all.expo.csv")
}




clump_exposure_SNP <- function(r2 = 0.001, kb = 10000){
  library(rio); library(tidyverse); library(TwoSampleMR)
  outTable = data.frame()
  expSNP <- read_exposure_data( filename = "all.expo.csv",  sep = ",", 
                                phenotype_col = "exposure", snp_col = "SNP",  beta_col = "beta",
                                se_col = "se", pval_col = "pval", effect_allele_col = "effect_allele",
                                other_allele_col = "other_allele",clump = F  )
  expo = levels(as.factor(expSNP$exposure)); n = length(expo)
  for (i in expo) {
    message(" ");
    message(paste0("------------------- [ ",n," ] left ------------------"))
    message(paste0(i," is running "))
    data = expSNP[expSNP$exposure == i,]
    #clumping the IVs
    clumped <- clump_data(data, clump_r2 = r2,clump_kb = kb)
    outTable = rbind(outTable,clumped)
    #write data
    export(outTable,"clumped.csv")
    expSNP = expSNP %>% filter(exposure != i)
    export(expSNP,"left.csv")
    n=n-1;  
    if(n==0){ message(" ");message("------------ run out of exposure ----------")}
  }
  export(outTable,paste0("expo_clumped(r2=",clump_r2," kb=",clump_kb,").csv"))
}

Extract_outcome_SNP <- function(exposureFile, endpointFile){
  #use library
  library(rio); library(tidyverse); library(TwoSampleMR)
  expSNP <-  import(exposureFile)
  endpoint <- import(endpointFile)
  N_SNP <-  length(unique(expSNP$SNP))
  
  outSNP = endpoint %>% filter(SNP %in% expSNP$SNP)
  lostSNP = expSNP %>% filter(!(SNP %in% outSNP$SNP))
  out_come_sig = outSNP %>% filter(pval <= 5e-8)
  outSNP = outSNP %>% filter(pval >= 5e-8)
  #export the outcome data
  export(lostSNP, "Del_SNP_not_in_outcome.csv")
  export(out_come_sig,"out_come_sig.CSV")
  export(outSNP, "outcome_SNPs.csv")
  
  outcome_dat <- read_outcome_data(
    snps = expSNP$SNP,filename = "outcome_SNPs.csv",se_col = "se", pval_col = "pval",
    phenotype_col = "exposure", sep = ",", snp_col = "SNP", beta_col = "beta",
    effect_allele_col = "effect_allele",other_allele_col = "other_allele"
  )
  #export the outcome data
  export(outcome_dat,"outcome_SNPs.csv")
}



Harmonise <- function(expoSNPFile, outSNPFile){
  #use library
  library(rio);library(tidyverse);library(TwoSampleMR)

  #import data
  exposure_dat = import(expoSNPFile); N_SNP = length(unique(expSNP$SNP))
  outcome_dat = import(outSNPFile)

  #harmonize the two data
  dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)

  #Scan palindromic SNP
  palindromic_del <- dat %>% filter(mr_keep == FALSE)
  mr_keep <- dat %>% filter(mr_keep == TRUE)
  SNPs = data.frame(SNP = unique(mr_keep$SNP))

  #export the data for MR performing.
  export(palindromic_del,"Del_palindromic.csv")
  export(mr_keep,"mr_keep.csv"); export(SNPs,"mr_keep_SNP.csv")
  
  #create file and directory for PhonoScan
  if(!dir.exists("pleiotropy")){  dir.create("pleiotropy")}
  file.create("pleiotropy_SNP.txt")
  
  n = ceiling( length(mr_keep$SNP)/100 )
  for (i in 1:n) {  do.call(file.create, str_c("pleiotropy\\SNP_",i,".txt"))}
  file.copy("mr_keep_SNP.csv","pleiotropy\\mr_keep_SNP.csv")
}



MR_expo_on_outc <- function(MR_data, method_list){
  #use library
  library(rio);library(tidyverse);library(TwoSampleMR)
  dat = import(MR_data)
  mr_res = mr(dat, method_list = method_list)
  
  mr_res_ivw = mr_res %>% filter(method=="Inverse variance weighted") %>%
    mutate(FDR = p.adjust(pval,method = "BH")) %>% arrange(pval)
  
  mr_plt = mr_pleiotropy_test(dat);  mr_het = mr_heterogeneity(dat)
  mr_LOO = mr_leaveoneout(dat);      mr_Odd = generate_odds_ratios(mr_res)
  mr_OR  = data.frame(Method = mr_Odd$method, OR=mr_Odd$or,
                      OR.95L=mr_Odd$or_lci95, OR.95H = mr_Odd$or_uci95, pvalue=mr_Odd$pval)
  
  export(mr_res_ivw, "mr_res_ivw.csv")
  export(mr_plt, "mr_plt.csv"); export(mr_het, "mr_het.csv")
  export(mr_LOO, "mr_LOO.csv"); export(mr_res, "mr_res.csv")
  export(mr_Odd, "mr_Odd.csv"); export(mr_OR,"mr_OR.csv") 
}

#MR_PRSSO
PRESSO <- function(MR_data){
  #use library
  library(rio); library(tidyverse); library(MRPRESSO)
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
}





