#00_prepare_data.R

#Data preparation for the project

#The purpose of this script is to source and prepare variables for analysis. 

#Most of the preparation of raw variables in the MoBa dataset and linked
#registry sources is done using the `phenotools` R package, code and 
#documentation for which is available at https://github.com/psychgen/phenotools.

#install.packages("//ess01/p471/data/durable/common/software/phenotools_0.3.2.zip", 
#         repos=NULL,
#         type = "binary")

#MoBa/MBRN variables required:

##Phenotypic data:

#  - self-reported age at menarche and pubertal stage 14 yr;
#  - depressive sx (SMFQ) 8 & 14 yr; 
#  - anxiety sx (SCARED) 8 & 14 yr; 
#  - conduct disorder sx (RS-DBD) 8 & 14 yr;
#  - oppositional defiant disorder (RS-DBD) sx 8 & 14 yr;
#  - ADHD traits (RS-DBD) 8 & 14 yr;
#  - covariates: BMI (derived) 8 & 14 yr; maternal/paternal age; child age at
#  - questionnaire completion 8 & 14 yr; parental education; parental income; 
#  - parental cohabitation; parity; financial problems; maternal prenatal 
#  - depression (SCL); maternal postnatal depression (EPD)

##Genetic instruments:

#  - age at menarche 10.1038/ng.3841
#  - childhood body size 10.1093/hmg/ddv472
#  - adult BMI 10.1093/hmg/ddy271

##Registry variables required:

#Primary care diagnoses (KUHR; ICPC-2):
#  - depressive disorders: P76 
#  - anxiety disorders: P74, P79, P82
#  - adhd: P81
#  - disruptive behaviour disorders: P23 
  
#Secondary care diagnoses (NPR; ICD-10):
#  - depressive disorders: F32-F33, F34.1 
#  - anxiety disorders: F40-F44, F93.0-F93.2
#  - adhd: F90
#  - disruptive behaviour disorders: F91, F92   
  
#Load required packages
library(phenotools)
library(tidyverse)
library(genotools)
library(MendelianRandomization)

##Curate dataset with variables from MoBa phenotypic data
covars <- c("KJONN","PARITET_5",
            "AGE_RETURN_MTHS_Q8AAR",
            "AGE_YRS_UB",
            "AA1315", "AA1316",
            "AA1124", "AA1126",
            "EE488", "GG383", "EE584",
            "FARS_ALDER", "MORS_ALDER", 
            "scl_short_m_q1",
            "scl_full_m_q3",
            "epds_short_m_6m",
            "bmi_derived_c_8yr",
            "bmi_derived_c_14c")

aux_vars <- c("VEKT","SVLEN_DG")

q8yr_vars <- c("smfq_dep_c_8yr",
               "scared_anx_c_8yr",
               "rsdbd_adhd_c_8yr",
               "rsdbd_cd_c_8yr",
               "rsdbd_odd_c_8yr")

q14yr_vars <- c("UB231",
                "pds_menarche_c_14c",
                "pds_breast_c_14c",
                "pds_growth_c_14c",
                "pds_hair_c_14c",
                "pds_skin_c_14c",
                "smfq_dep_c_14c",
                "scared_anx_c_14c",
                "rsdbd_adhd_c_14m",
                "rsdbd_cd_c_14c",
                "rsdbd_odd_c_14m")

if(!file.exists("data/pheno_data.RData")){
  pheno_data <- curate_dataset(variables_required=list(moba=c(covars,
                                                              aux_vars,
                                                              q8yr_vars,
                                                              q14yr_vars),
                                                       npr=c("dep=F32,F33,F341",
                                                             "anx=F40,F41,F42,F43,F44,F930,F931,F932",
                                                             "adhd=F90",
                                                             "beh=F91,F92"),
                                                       kuhr=c("dep=P76",
                                                              "anx=P74, P79, P82",
                                                              "adhd=P81",
                                                              "beh=P23")),
                               out_format = "merged_df",
                               recursive=T,
                               dx_owners="child")
  pheno_data <- pheno_data %>% 
    rename("sex"=KJONN_raw) %>%
    filter(sex ==2)
  
  save(pheno_data, file="data/pheno_data.RData")
}else{
  load(file="data/pheno_data.RData")
}

#source script to process the diagnostic data
if(!file.exists("data/dx_data.RData")){
  source("./scripts/00.1_process_dx_data.R")
}else{
  load(file="./data/dx_data.RData")
}

##Join questionnaire data with processed diagnostic data

#select/rename phenotypic variables and create 
#individual identifier ("ind_id") for merging 
pheno_data <- pheno_data %>%
  select(preg_id,BARN_NR,m_id,sex,birth_yr,parity=PARITET_5_raw,
         age_8yr=AGE_RETURN_MTHS_Q8AAR_raw,age_14yr=AGE_YRS_UB_raw,
         m_age_at_birth=MORS_ALDER_raw,p_age_at_birth=FARS_ALDER_raw,
         parent_cohab_18m=EE488_raw,parent_cohab_3yr=GG383_raw,
         smfq_dep_c_8yr,scared_anx_c_8yr,rsdbd_adhd_c_8yr,rsdbd_cd_c_8yr,
         rsdbd_odd_c_8yr,smfq_dep_c_14c,scared_anx_c_14c,rsdbd_cd_c_14c,
         rsdbd_odd_c_14m,rsdbd_adhd_c_14m,financ_probs_18m=EE584_raw,
         AA1315_raw,AA1316_raw,AA1124_raw,AA1126_raw,scl_short_m_q1,
         scl_full_m_q3,epds_short_m_6m,bmi_derived_c_8yr,bmi_derived_c_14c,
         menarche_c_14c=pds_menarche_c_14c,breast_c_14c=pds_breast_c_14c,
         growth_c_14c=pds_growth_c_14c,skin_c_14c=pds_skin_c_14c,
         hair_c_14c=pds_hair_c_14c,aam_c_14c=UB231_raw) %>%
  mutate(preg_id = as.numeric(preg_id)) %>%
  mutate(ind_id = paste0(preg_id,"_",BARN_NR))

#summarise yearwise child diagnostic data into two periods: ch = 0-8; ad = 10-17
summary_dx <- dx_data %>%
  right_join(pheno_data %>% select(ind_id)) %>%  #first restrict to those in the analytic dataset to preserve RAM
  pivot_longer(matches("_dx_")) %>% 
  separate(name, into=c("pheno","null","age"), sep="_") %>% 
  mutate(age=as.numeric(age),
         period = case_when(age<9 ~ "ch",
                            age>9 ~ "ad",
                            TRUE  ~ NA_character_)) %>% 
  drop_na(period) %>% 
  group_by(ind_id, pheno, period) %>% 
  summarise(n_dx=sum(value,na.rm=T),
            n_na=sum(is.na(value)), 
            prop_na=sum(is.na(value))/n())

#Below we implement the rule that IN CHILDHOOD, any non-missing follow-up is sufficient to get a non-missing summary, because with left-censoring in the context of an ever/never summary, any info 
#is informative as to the correct value for the summary (i.e., diagnosis or no diagnosis at age 8, even if all prior years are missing, is representative of the value individuals should get in
#the summary in all but rare cases where diagnoses are given early in childhood but not repeated); in contrast, IN ADOLESCENCE, individuals must have follow up for at least 5 out of 8 years to
#receive a non-missing summary value unless the available years include a diagnosis (because the data are right censored, meaning that missing years can only change the ever/never summary value
#an individual should receive if they include new diagnoses - the likelihood of which diminishes the more available years an individual has with no diagnosis)

dx_data_reduced <- summary_dx %>% 
  mutate(dx=ifelse(n_dx==0,0,1),
         dx=ifelse(prop_na==1,NA,dx),
         dx=ifelse(prop_na>0.6 & n_dx==0 & period=="ad" , NA, dx)) %>% 
  select(-n_dx,-prop_na, -n_na) %>% 
  unite("name",pheno, period, sep="_dx_") %>% 
  pivot_wider(id_cols=ind_id, names_from = name, values_from = dx)

#combine with pheno_data
pheno_data <- pheno_data %>%
  left_join(dx_data_reduced)

rm(dx_data)

##Manipulate variables

#change response ordering and derive education variable (combine '2' & '3')
pheno_data <- pheno_data %>%
  mutate(education_m_q1 = case_when(AA1124_raw==1~1,
                                    AA1124_raw==2~2, 
                                    AA1124_raw==3~2,
                                    AA1124_raw==4~3,
                                    AA1124_raw==5~4,
                                    AA1124_raw==6~5,
                                    AA1124_raw==0~NA_real_)) %>%
  mutate(education_f_q1 = case_when(AA1126_raw==1~1,
                                    AA1126_raw==2~2, 
                                    AA1126_raw==3~2,
                                    AA1126_raw==4~3,
                                    AA1126_raw==5~4,
                                    AA1126_raw==6~5,
                                    AA1126_raw==0~NA_real_)) %>%
  rowwise() %>%
  mutate(parent_education_q1 = mean(c(education_m_q1,education_f_q1),na.rm=TRUE)) %>%
  ungroup() %>% 
  mutate_all(~ifelse(is.nan(.), NA, .))

#recode and compute derived income variable (set '0' & '8' to NA)
pheno_data <- pheno_data %>%
  mutate(income_m_q1 = case_when(AA1315_raw==1~1,
                                 AA1315_raw==2~2, 
                                 AA1315_raw==3~3,
                                 AA1315_raw==4~4,
                                 AA1315_raw==5~5,
                                 AA1315_raw==6~6,
                                 AA1315_raw==7~7,
                                 AA1315_raw==0~NA_real_)) %>%
  mutate(income_f_q1 = case_when(AA1316_raw==1~1,
                                 AA1316_raw==2~2, 
                                 AA1316_raw==3~3,
                                 AA1316_raw==4~4,
                                 AA1316_raw==5~5,
                                 AA1316_raw==6~6,
                                 AA1316_raw==7~7, 
                                 AA1316_raw==8~NA_real_,
                                 AA1316_raw==0~NA_real_)) %>%
  rowwise() %>%
  mutate(parent_income_q1 = mean(c(income_m_q1,income_f_q1),na.rm=TRUE)) %>%
  ungroup() %>% 
  mutate_all(~ifelse(is.nan(.), NA, .))

#recode parental cohabitation and financial prob variables
pheno_data <- pheno_data %>%
  mutate(parent_cohab_18m = case_when(parent_cohab_18m == 0 ~ NA_real_,
                                      parent_cohab_18m == 1 ~ 0,
                                      parent_cohab_18m == 2 ~ 1),
         parent_cohab_3yr = case_when(parent_cohab_3yr == 0 ~ NA_real_,
                                      parent_cohab_3yr == 1 ~ 1,
                                      parent_cohab_3yr == 2 ~ 0),
         financ_probs_18m = case_when(financ_probs_18m == 0 ~ NA_real_,
                                      financ_probs_18m == 1 ~ 0,
                                      financ_probs_18m == 2 ~ 1,
                                      financ_probs_18m == 3 ~ 2,
                                      financ_probs_18m == 4 ~ 3))

#filter to 14-year responders (based on the SMFQ)
adol_data <- pheno_data %>%
  filter(!is.na(smfq_dep_c_14c)) %>%
  select(-sex,-matches(c("_raw","education_m_q1","education_f_q1",
                         "income_m_q1","income_f_q1"))) 

#remove outliers from items not processed through phenotools
vars <- c("m_age_at_birth", "p_age_at_birth")

adol_data <- adol_data %>%
  mutate(across(all_of(vars), ~replace(.,abs(.-mean(.,na.rm=T)) >3*sd(.,na.rm=T),NA)))

#save out
save(adol_data, file="./data/adol_data.RData")

##Fetch and process genetic instrument

#source script to create genetic instrument in PRSice2 (note that some of 
#the script needs to be run on the cluster initially)
source("./scripts/00.2_construct_grs.R")

##Join genetic and phenotypic data for complete analytic dataset

#load genetic data
load(file="data/pgs_procd.RData")
load(file="./data/adol_data.RData")

#merge with genetic data
fulldata <- adol_data %>%
  left_join(pgs_procd) 

#save out analysis dataset (before imputation)
save(fulldata, file="data/analysis_dataset.RData")

##Conduct multiple imputation

#source MI script
source("./scripts/00.3_multiple_imputation.R")

#process MI data (transform and standardise)
source("./scripts/00.4_process_imputed_data.R")

##Prepare GWAS summary statistics for two-sample MR

#create covariate and phenotype input files for Regenie
source("./scripts/00.5_gwas_twosample.R")

#the Regenie GWAS is run on the cluster (shell scripts in `GWAS/input_files/`)

#process and plot GWAS results, saving out the sumstats
source("./scripts/00.6_process_gwas_data.R")

##Validate genetic instrument for age at menarche

#load processed MI data
load(file="./data/standardised_MI_tr.RData")

#check prediction of aam across multiply imputed datasets
model_aam <- with(standardised_MI, lm(aam_c_14c~PGS_AaM))

summary(pool(model_aam))

pool.r.squared(model_aam)#R2 = 0.0692

#calculate mean F-statistic across imputed datasets
#as an approximation, and R-squared for validation

#load list of imputed datasets
load(file="data/datlist_tr.RData")

#make vector to store F-stats and R-squared
f_stats <- numeric(length(datlist))
r_squared <- numeric(length(datlist))

#loop through each imputed dataset in the list
for(i in seq_along(datlist)) {
  #fit the linear model on each imputed dataset
  model <- lm(aam_c_14c~PGS_AaM, data = datlist[[i]])
  
  #extract the F-statistic
  f_stats[i] <- summary(model)$fstatistic[1]
  r_squared[i] <- summary(model)$r.squared
}

#calculate the mean of the F-statistics
mean(f_stats)#F = 996.4
mean(r_squared)#R2 = 0.0692 (nearly identical)

#source script for prediction of covariates using the
#genetic instrument for age at menarche
source("./scripts/00.7_sensitivity_analyses_covs.R")

#load results from running the sensitivity analyses
load(file="output/sensitivity_cov_results.RData")

#print results
print(all_results)
