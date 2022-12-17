library(TwoSampleMR)
id_exposure <- c("ukb-b-3957",
                 "ukb-b-5776", 
                 "ukb-b-4424", 
                 "ukb-b-4956", 
                 "ukb-b-4616",
                 "ukb-b-17400",
                 "ukb-b-2772") 
id_outcome <- "finn-b-K11_CHOLELITH"
expo = extract_instruments("ukb-b-13806",
                    p1 = 5e-08,
                    clump = T,
                    p2 = 5e-08,
                    r2 = 0.001,
                    kb = 10000,
                    access_token = ieugwasr::check_access_token(),
                    force_server = FALSE)


outc = extract_outcome_data(expo$SNP,
                     id_outcome,
                     proxies = T,
                     rsq = 0.8,
                     align_alleles = 1,
                     palindromes = 1,
                     maf_threshold = 0.3,
                     access_token = ieugwasr::check_access_token(),
                     splitsize = 10000,
                     proxy_splitsize = 500)
dat = harmonise_data(exposure_dat = expo, outcome_dat = outc, action = 2)


mr_res = mr(dat)

ivw = mr_res %>% filter(method == "Inverse variance weighted") %>% 
                 mutate(FDR = p.adjust(pval,method = "BH"))

mr_res














mr_plt = mr_pleiotropy_test(dat)
mr_het = mr_heterogeneity(dat)
PRESSO = run_mr_presso(dat, NbDistribution = 10000, SignifThreshold = 0.05)
PRESSO

adj_dat = dat[-c(40),]


adj_mr_res = mr(adj_dat)
adj_mr_plt = mr_pleiotropy_test(adj_dat)
adj_mr_het = mr_heterogeneity(adj_dat)
adj_PRESSO = run_mr_presso(adj_dat, NbDistribution = 10000, SignifThreshold = 0.05)
adj_PRESSO


