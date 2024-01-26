#00.5_gwas_twosample.R

##This script creates covariate and phenotype input files for Regenie;
##GWAS based on these files are run on the cluster

#load required packages
library(tidyverse)
library(data.table)
library(purrr)

#read in files
load(file="data/datlist_tr.RData")
linkage_IDs <- fread("N:/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt")

#take mean of imputed values for missing data on symptom phenotypes
#due to missing values in mother-reported outcomes (adhd, odd)
domains <- c("smfq_dep_c_14c","scared_anx_c_14c","rsdbd_cd_c_14c","rsdbd_odd_c_14m","rsdbd_adhd_c_14m")

means <- map(datlist, ~ .x %>%
             group_by(ind_id) %>%
             summarise(across(all_of(domains), ~ mean(., na.rm = TRUE), .names = "mean_{.col}")))

combined_means <- bind_rows(means)

final_means <- combined_means %>%
  group_by(ind_id) %>%
  summarise(across(starts_with("mean"), mean, na.rm = TRUE))

#create covariate file for GWAS
covar_gwas <- merge(final_means, linkage_IDs, by.x="ind_id", by.y="ID_2306", all.x=T) %>%
  select(FID,IID,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,
         PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,
         batch_harvest12a,batch_harvest12b,batch_harvest24,batch_rotterdam1,
         batch_rotterdam2,batch_adhd1,batch_adhd2,batch_norment_may2016,
         batch_norment_feb2018,batch_norment_feb2020_v1,batch_norment_feb2020_v3,
         batch_norment_aug2020_996,batch_norment_aug2020_1029,
         batch_norment_nov2020_1066,batch_norment_nov2020_1077,
         batch_norment_nov2020_1108,batch_norment_nov2020_1109,
         batch_norment_nov2020_1135,batch_norment_nov2020_1146,
         batch_norment_mar2021_1273,batch_norment_mar2021_1409,
         batch_norment_mar2021_1413,batch_norment_mar2021_1531) %>%
  na.omit()

#create phenotype file for dep GWAS
pheno_gwas_dep <- merge(final_means, linkage_IDs, by.x="ind_id", by.y="ID_2306", all.x=T) %>%
  select(FID,IID,mean_smfq_dep_c_14c) %>%
  na.omit() %>%
  mutate(mean_smfq_dep_c_14c = scale(mean_smfq_dep_c_14c)) %>%
  mutate(mean_smfq_dep_c_14c = as.numeric(mean_smfq_dep_c_14c))

#create phenotype file for anx GWAS
pheno_gwas_anx <- merge(final_means, linkage_IDs, by.x="ind_id", by.y="ID_2306", all.x=T) %>%
  select(FID,IID,mean_scared_anx_c_14c) %>%
  na.omit() %>%
  mutate(mean_scared_anx_c_14c = scale(mean_scared_anx_c_14c)) %>%
  mutate(mean_scared_anx_c_14c = as.numeric(mean_scared_anx_c_14c))

#create phenotype file for cd GWAS
pheno_gwas_cd <- merge(final_means, linkage_IDs, by.x="ind_id", by.y="ID_2306", all.x=T) %>%
  select(FID,IID,mean_rsdbd_cd_c_14c) %>%
  na.omit() %>%
  mutate(mean_rsdbd_cd_c_14c = scale(mean_rsdbd_cd_c_14c)) %>%
  mutate(mean_rsdbd_cd_c_14c = as.numeric(mean_rsdbd_cd_c_14c))

#create phenotype file for odd GWAS
pheno_gwas_odd <- merge(final_means, linkage_IDs, by.x="ind_id", by.y="ID_2306", all.x=T) %>%
  select(FID,IID,mean_rsdbd_odd_c_14m) %>%
  na.omit() %>%
  mutate(mean_rsdbd_odd_c_14m = scale(mean_rsdbd_odd_c_14m)) %>%
  mutate(mean_rsdbd_odd_c_14m = as.numeric(mean_rsdbd_odd_c_14m))

#create phenotype file for adhd GWAS
pheno_gwas_adhd <- merge(final_means, linkage_IDs, by.x="ind_id", by.y="ID_2306", all.x=T) %>%
  select(FID,IID,mean_rsdbd_adhd_c_14m) %>%
  na.omit() %>%
  mutate(mean_rsdbd_adhd_c_14m = scale(mean_rsdbd_adhd_c_14m)) %>%
  mutate(mean_rsdbd_adhd_c_14m = as.numeric(mean_rsdbd_adhd_c_14m))

#write out covariate and phenotype files
write_delim(covar_gwas,file="./scripts/GWAS/input_files/covar_gwas.txt",col_names=TRUE)
write_delim(pheno_gwas_dep,file="./scripts/GWAS/input_files/pheno_gwas_dep_tr.txt",col_names=TRUE)
write_delim(pheno_gwas_anx,file="./scripts/GWAS/input_files/pheno_gwas_anx_tr.txt",col_names=TRUE)
write_delim(pheno_gwas_cd,file="./scripts/GWAS/input_files/pheno_gwas_cd_tr.txt",col_names=TRUE)
write_delim(pheno_gwas_odd,file="./scripts/GWAS/input_files/pheno_gwas_odd_tr.txt",col_names=TRUE)
write_delim(pheno_gwas_adhd,file="./scripts/GWAS/input_files/pheno_gwas_adhd_tr.txt",col_names=TRUE)

