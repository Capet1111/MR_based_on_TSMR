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


