#04.1_aim4_exploratory_analyses.R

##This script conducts exploratory analyses primarily related to Aim 4.
##This includes adding adolescent depression dx and pre-pubertal dx for 
##each domain as covariates, and running one-sample MR analysis with dx
##outcomes divided into pre-adolescence/early/mid-late adolescence.

#Read in the data and load packages:
load(file="data/datlist_tr.RData")
library(tidyverse)
library(ivreg)
library(sandwich)
library(effectsize)

#Function to perform two-stage regression and adjust SEs for adolescent dx
two_stage_with_adjSEs_ad <- function(data) {
  
  # First stage: linear regression
  lm_model <- lm(aam_c_14c ~ PGS_AaM, data = data)
  predictions <- predict(lm_model, newdata = data)
  data$pred_aam <- predictions
  
  # Second stage: logistic regression with the predictions from the first stage
  logit_model <- glm(as.formula(paste(domain[[2]], "~ pred_aam + dep_dx_ad + ", domain[[1]])), family = binomial(link = "logit"), data = data)
  
  # Get robust standard errors accounting for uncertainty in first stage
  robust_se <- sqrt(diag(vcovHC(logit_model, type = "HC1"))['pred_aam'])
  
  # Extract coefficient for the predictor from the second stage
  estimate_iv <- coef(logit_model)['pred_aam']
  
  # Calculate z-value, two-tailed p-value, and confidence interval for the IV estimate
  z_iv <- estimate_iv / robust_se
  pval_iv <- 2 * pnorm(abs(z_iv), lower.tail = FALSE)
  ci_low_iv <- estimate_iv - qnorm(0.975) * robust_se
  ci_high_iv <- estimate_iv + qnorm(0.975) * robust_se
  
  # Combine IV results into a data frame
  iv_results <- data.frame(
    Estimate_IV = estimate_iv,
    Robust_SE_IV = robust_se,
    Z_IV = z_iv,
    P_Value_IV = pval_iv,
    CI_Low_IV = ci_low_iv,
    CI_High_IV = ci_high_iv
  )
  
  return(iv_results)
}

#list the diagnosis variables
domains <-list(anx=c("anx_dx_ch","anx_dx_ad"),
               adhd = c("adhd_dx_ch","adhd_dx_ad"),
               beh=c("beh_dx_ch","beh_dx_ad"))

#run adolescent dx models for all domains across the 50 imputed datasets
domains_res_ad <- list()
n = 0

for(domain in domains) {
  n = n + 1
  
  # apply this function to each imputed dataset
  results_list <- lapply(datlist, function(data) two_stage_with_adjSEs_ad(data))
  
  m <- length(datlist)
  
  # pool the results from each imputed dataset using Rubin's rules 
  pooled_estimate <- mean(sapply(results_list, function(x) x$Estimate_IV))
  within_var <- mean(sapply(results_list, function(x) x$Robust_SE_IV)^2)*m/(m - 1)
  between_var <- var(sapply(results_list, function(x) x$Estimate_IV))*(m - 1)/m
  total_var <- within_var + (1 + 1/m) * between_var
  pooled_se <- sqrt(total_var)
  
  #calculate degrees of freedom using Barnard and Rubin's adjustment
  df <- (m - 1) * (1 + (within_var / total_var))^2
  
  #compute test statistics for the pooled estimate
  z_val <- pooled_estimate / pooled_se
  
  #calculate one-tailed and two-tailed p-values
  two_tailed_p <- 2*(pt(-abs(z_val), df))
  
  #constructing pooled results
  pooled_results <- data.frame(
    Estimate = pooled_estimate,
    se = pooled_se,
    z = pooled_estimate / pooled_se,
    p = two_tailed_p,
    CI_low = pooled_estimate - qnorm(0.975) * pooled_se,
    CI_high = pooled_estimate + qnorm(0.975) * pooled_se,
    ROPE_low = d_to_oddsratio(-0.20),
    ROPE_high = d_to_oddsratio(0.20)
  )
  
  domains_res_ad[[names(domains)[[n]]]] <- pooled_results
  
}

#reformat adolescent dx results
anx_dx_ad <- domains_res_ad$anx %>%
  mutate(domain = "ANX")
dbd_dx_ad <- domains_res_ad$beh %>%
  mutate(domain = "DBD")
adhd_dx_ad <- domains_res_ad$adhd %>%
  mutate(domain = "ADHD")

#exponentiate the diagnostic outcome results ad
ests_4b_ad_exp <- anx_dx_ad %>%
  full_join(dbd_dx_ad) %>%
  full_join(adhd_dx_ad) %>%
  mutate(Estimate = exp(Estimate),
         CI_low = exp(CI_low),
         CI_high = exp(CI_high),
         outcome = "Diagnoses age 10-17")

#print results
print(ests_4b_ad_exp)

#exploration of the causal effect on diagnostic outcomes in adolescence

#load data
load(file="data/pheno_data.RData")
load(file="./data/dx_data.RData")

#create individual identifier ("ind_id") for merging 
pheno_data <- pheno_data %>%
  select(preg_id,BARN_NR,m_id,birth_yr,aam_c_14c=UB231_raw,
         menarche_c_14c=pds_menarche_c_14c,smfq_dep_c_14c) %>%
  mutate(preg_id = as.numeric(preg_id)) %>%
  mutate(ind_id = paste0(preg_id,"_",BARN_NR))

#summarise yearwise adolescent diagnostic data into three periods: pre=9-11, early=12-14, and mid-late=15-17
summary_dx <- dx_data %>%
  right_join(pheno_data %>% select(ind_id)) %>%  #first restrict to those in the analytic dataset to preserve RAM
  pivot_longer(matches("_dx_")) %>% 
  separate(name, into=c("pheno","null","age"), sep="_") %>% 
  mutate(age=as.numeric(age),
         period = case_when(age<9 ~ "ch",
                            age>8 & age<12 ~ "pre",
                            age>11 & age<15 ~ "early",
                            age>14 ~ "late",
                            TRUE  ~ NA_character_)) %>% 
  drop_na(period) %>% 
  group_by(ind_id, pheno, period) %>% 
  summarise(n_dx=sum(value,na.rm=T),
            n_na=sum(is.na(value)), 
            prop_na=sum(is.na(value))/n())

#transform dataset to wide
dx_data_reduced <- summary_dx %>% 
  mutate(dx=ifelse(n_dx==0,0,1),
         dx=ifelse(prop_na==1,NA,dx),
         dx=ifelse(prop_na>0.6 & n_dx==0 & period=="late" , NA, dx)) %>% 
  select(-n_dx,-prop_na, -n_na) %>% 
  unite("name",pheno, period, sep="_dx_") %>% 
  pivot_wider(id_cols=ind_id, names_from = name, values_from = dx)

#combine with pheno_data
pheno_data <- pheno_data %>%
  left_join(dx_data_reduced)

rm(dx_data)

#filter to 14-year responders (based on the SMFQ)
adol_data_sensitivity <- pheno_data %>%
  filter(!is.na(smfq_dep_c_14c)) %>%
  select(-smfq_dep_c_14c,-matches("_raw"))

#load genetic data
load(file="data/pgs_procd.RData")

#merge with genetic data
fulldata_sensitivity <- adol_data_sensitivity %>%
  left_join(pgs_procd) 

#set age at menarche >14 to 15 for the purpose of this sensitivity analysis,
#to avoid having to do imputation, standardise the genetic instrument for age
#at menarche and age at menarche, and set NAs in diagnostic variables to 0
fulldata_sensitivity <- fulldata_sensitivity %>%
  mutate(aam_c_14c = as.numeric(aam_c_14c),
         PGS_AaM = scale(PGS_AaM)) %>%
  mutate(aam_c_14c = case_when(menarche_c_14c == 0 ~ 9,
                               menarche_c_14c == 1 ~ aam_c_14c,
                               TRUE ~ NA_real_)) %>%
  mutate(aam_c_14c = scale(aam_c_14c)) %>%
  mutate(across(contains("_dx_"), ~ifelse(is.na(.), 0, .)))

fulldata_counts_adhd <- fulldata_sensitivity %>%
  mutate(
    # New ADHD diagnoses in the pre-adolescent stage, removing childhood diagnoses
    new_adhd_dx_pre = ifelse(adhd_dx_pre == 1 & adhd_dx_ch == 0, 1, 0),
    
    # New ADHD diagnoses in early adolescence, removing childhood and pre-adolescent diagnoses
    new_adhd_dx_early = ifelse(adhd_dx_early == 1 & adhd_dx_ch == 0 & adhd_dx_pre == 0, 1, 0),
    
    # New ADHD diagnoses in mid-late adolescence, removing all previous stages
    new_adhd_dx_late = ifelse(adhd_dx_late == 1 & adhd_dx_ch == 0 & adhd_dx_pre == 0 & adhd_dx_early == 0, 1, 0)
  )

# To get the counts of new diagnoses for each stage ADHD
new_diagnoses_counts_adhd <- fulldata_counts_adhd %>%
  summarise(
    new_adhd_dx_pre_count = sum(new_adhd_dx_pre, na.rm = TRUE),
    new_adhd_dx_early_count = sum(new_adhd_dx_early, na.rm = TRUE),
    new_adhd_dx_late_count = sum(new_adhd_dx_late, na.rm = TRUE)
  )

# Print the result
print(new_diagnoses_counts_adhd)


#Run two-stage IV regression analysis for ADHD diagnoses in pre/early/mid-late adolescence

# First stage: linear regression
lm_model <- lm(aam_c_14c ~ PGS_AaM, data = fulldata_sensitivity)
predictions <- predict(lm_model, newdata = fulldata_sensitivity)
fulldata_sensitivity$pred_aam <- predictions

# Second stage: logistic regression with the predictions from the first stage
logit_model1 <- glm(as.formula("adhd_dx_pre ~ pred_aam + adhd_dx_ch"), family = binomial(link = "logit"), data = fulldata_sensitivity)
logit_model2 <- glm(as.formula("adhd_dx_early ~ pred_aam + adhd_dx_ch + adhd_dx_pre"), family = binomial(link = "logit"), data = fulldata_sensitivity)
logit_model3 <- glm(as.formula("adhd_dx_late ~ pred_aam + adhd_dx_ch + adhd_dx_pre + adhd_dx_early"), family = binomial(link = "logit"), data = fulldata_sensitivity)

# Extract coefficient for the predictor from the second stage
estimate_iv1 <- coef(logit_model1)['pred_aam']
estimate_iv2 <- coef(logit_model2)['pred_aam']
estimate_iv3 <- coef(logit_model3)['pred_aam']

# Get robust standard errors accounting for uncertainty in first stage
robust_se_iv1 <- sqrt(diag(vcovHC(logit_model1, type = "HC1"))['pred_aam'])
robust_se_iv2 <- sqrt(diag(vcovHC(logit_model2, type = "HC1"))['pred_aam'])
robust_se_iv3 <- sqrt(diag(vcovHC(logit_model3, type = "HC1"))['pred_aam'])

# Calculate z-value, two-tailed p-value, and confidence interval for the IV estimate
z_iv1 <- estimate_iv1 / robust_se_iv1
pval_iv1 <- 2 * pnorm(abs(z_iv1), lower.tail = FALSE)
ci_low_iv1 <- estimate_iv1 - qnorm(0.975) * robust_se_iv1
ci_high_iv1 <- estimate_iv1 + qnorm(0.975) * robust_se_iv1

z_iv2 <- estimate_iv2 / robust_se_iv2
pval_iv2 <- 2 * pnorm(abs(z_iv2), lower.tail = FALSE)
ci_low_iv2 <- estimate_iv2 - qnorm(0.975) * robust_se_iv2
ci_high_iv2 <- estimate_iv2 + qnorm(0.975) * robust_se_iv2

z_iv3 <- estimate_iv3 / robust_se_iv3
pval_iv3 <- 2 * pnorm(abs(z_iv3), lower.tail = FALSE)
ci_low_iv3 <- estimate_iv3 - qnorm(0.975) * robust_se_iv3
ci_high_iv3 <- estimate_iv3 + qnorm(0.975) * robust_se_iv3

# Combine IV results into data frames
iv_results_pre <- data.frame(
  Estimate_IV = exp(estimate_iv1),
  Robust_SE_IV = robust_se_iv1,
  Z_IV = z_iv1,
  P_Value_IV = pval_iv1,
  CI_Low_IV = exp(ci_low_iv1),
  CI_High_IV = exp(ci_high_iv1),
  Outcome = "Diagnoses age 9-11")

iv_results_early <- data.frame(
  Estimate_IV = exp(estimate_iv2),
  Robust_SE_IV = robust_se_iv2,
  Z_IV = z_iv2,
  P_Value_IV = pval_iv2,
  CI_Low_IV = exp(ci_low_iv2),
  CI_High_IV = exp(ci_high_iv2),
  Outcome = "Diagnoses age 12-14")

iv_results_mid_late <- data.frame(
  Estimate_IV = exp(estimate_iv3),
  Robust_SE_IV = robust_se_iv3,
  Z_IV = z_iv3,
  P_Value_IV = pval_iv3,
  CI_Low_IV = exp(ci_low_iv3),
  CI_High_IV = exp(ci_high_iv3),
  Outcome = "Diagnoses age 15-17")
         
#print results
print(iv_results_pre)
print(iv_results_early)
print(iv_results_mid_late)

#check new diagnoses at each stage DEP
fulldata_counts_dep <- fulldata_sensitivity %>%
  mutate(
    new_dep_dx_pre = ifelse(dep_dx_pre == 1 & dep_dx_ch == 0, 1, 0),
    new_dep_dx_early = ifelse(dep_dx_early == 1 & dep_dx_ch == 0 & dep_dx_pre == 0, 1, 0),
    new_dep_dx_late = ifelse(dep_dx_late == 1 & dep_dx_ch == 0 & dep_dx_pre == 0 & dep_dx_early == 0, 1, 0)
  )

# To get the counts of new diagnoses for each stage
new_diagnoses_counts_dep <- fulldata_counts_dep %>%
  summarise(
    new_dep_dx_pre_count = sum(new_dep_dx_pre, na.rm = TRUE),
    new_dep_dx_early_count = sum(new_dep_dx_early, na.rm = TRUE),
    new_dep_dx_late_count = sum(new_dep_dx_late, na.rm = TRUE)
  )

# Print the result
print(new_diagnoses_counts_dep)

#Run two-stage IV regression analysis for DEP diagnoses in pre/early/mid-late adolescence

# First stage: linear regression
lm_model <- lm(aam_c_14c ~ PGS_AaM, data = fulldata_sensitivity)
predictions <- predict(lm_model, newdata = fulldata_sensitivity)
fulldata_sensitivity$pred_aam <- predictions

# Second stage: logistic regression with the predictions from the first stage
logit_model1 <- glm(as.formula("dep_dx_pre ~ pred_aam + dep_dx_ch"), family = binomial(link = "logit"), data = fulldata_sensitivity)
logit_model2 <- glm(as.formula("dep_dx_early ~ pred_aam + dep_dx_ch + dep_dx_pre"), family = binomial(link = "logit"), data = fulldata_sensitivity)
logit_model3 <- glm(as.formula("dep_dx_late ~ pred_aam + dep_dx_ch + dep_dx_pre + dep_dx_early"), family = binomial(link = "logit"), data = fulldata_sensitivity)

# Extract coefficient for the predictor from the second stage
estimate_iv1 <- coef(logit_model1)['pred_aam']
estimate_iv2 <- coef(logit_model2)['pred_aam']
estimate_iv3 <- coef(logit_model3)['pred_aam']

# Get robust standard errors accounting for uncertainty in first stage
robust_se_iv1 <- sqrt(diag(vcovHC(logit_model1, type = "HC1"))['pred_aam'])
robust_se_iv2 <- sqrt(diag(vcovHC(logit_model2, type = "HC1"))['pred_aam'])
robust_se_iv3 <- sqrt(diag(vcovHC(logit_model3, type = "HC1"))['pred_aam'])

# Calculate z-value, two-tailed p-value, and confidence interval for the IV estimate
z_iv1 <- estimate_iv1 / robust_se_iv1
pval_iv1 <- 2 * pnorm(abs(z_iv1), lower.tail = FALSE)
ci_low_iv1 <- estimate_iv1 - qnorm(0.975) * robust_se_iv1
ci_high_iv1 <- estimate_iv1 + qnorm(0.975) * robust_se_iv1

z_iv2 <- estimate_iv2 / robust_se_iv2
pval_iv2 <- 2 * pnorm(abs(z_iv2), lower.tail = FALSE)
ci_low_iv2 <- estimate_iv2 - qnorm(0.975) * robust_se_iv2
ci_high_iv2 <- estimate_iv2 + qnorm(0.975) * robust_se_iv2

z_iv3 <- estimate_iv3 / robust_se_iv3
pval_iv3 <- 2 * pnorm(abs(z_iv3), lower.tail = FALSE)
ci_low_iv3 <- estimate_iv3 - qnorm(0.975) * robust_se_iv3
ci_high_iv3 <- estimate_iv3 + qnorm(0.975) * robust_se_iv3

# Combine IV results into data frames
iv_results_pre <- data.frame(
  Estimate_IV = exp(estimate_iv1),
  Robust_SE_IV = robust_se_iv1,
  Z_IV = z_iv1,
  P_Value_IV = pval_iv1,
  CI_Low_IV = exp(ci_low_iv1),
  CI_High_IV = exp(ci_high_iv1),
  Outcome = "Diagnoses age 9-11")

iv_results_early <- data.frame(
  Estimate_IV = exp(estimate_iv2),
  Robust_SE_IV = robust_se_iv2,
  Z_IV = z_iv2,
  P_Value_IV = pval_iv2,
  CI_Low_IV = exp(ci_low_iv2),
  CI_High_IV = exp(ci_high_iv2),
  Outcome = "Diagnoses age 12-14")

iv_results_mid_late <- data.frame(
  Estimate_IV = exp(estimate_iv3),
  Robust_SE_IV = robust_se_iv3,
  Z_IV = z_iv3,
  P_Value_IV = pval_iv3,
  CI_Low_IV = exp(ci_low_iv3),
  CI_High_IV = exp(ci_high_iv3),
  Outcome = "Diagnoses age 15-17")

#print results
print(iv_results_pre)
print(iv_results_early)
print(iv_results_mid_late)

#check new diagnoses at each stage ANX
fulldata_counts_anx <- fulldata_sensitivity %>%
  mutate(
    new_anx_dx_pre = ifelse(anx_dx_pre == 1 & anx_dx_ch == 0, 1, 0),
    new_anx_dx_early = ifelse(anx_dx_early == 1 & anx_dx_ch == 0 & anx_dx_pre == 0, 1, 0),
    new_anx_dx_late = ifelse(anx_dx_late == 1 & anx_dx_ch == 0 & anx_dx_pre == 0 & anx_dx_early == 0, 1, 0)
  )

# To get the counts of new diagnoses for each stage
new_diagnoses_counts_anx <- fulldata_counts_anx %>%
  summarise(
    new_anx_dx_pre_count = sum(new_anx_dx_pre, na.rm = TRUE),
    new_anx_dx_early_count = sum(new_anx_dx_early, na.rm = TRUE),
    new_anx_dx_late_count = sum(new_anx_dx_late, na.rm = TRUE)
  )

# Print the result
print(new_diagnoses_counts_anx)

#Run two-stage IV regression analysis for anx diagnoses in pre/early/mid-late adolescence

# First stage: linear regression
lm_model <- lm(aam_c_14c ~ PGS_AaM, data = fulldata_sensitivity)
predictions <- predict(lm_model, newdata = fulldata_sensitivity)
fulldata_sensitivity$pred_aam <- predictions

# Second stage: logistic regression with the predictions from the first stage
logit_model1 <- glm(as.formula("anx_dx_pre ~ pred_aam + anx_dx_ch"), family = binomial(link = "logit"), data = fulldata_sensitivity)
logit_model2 <- glm(as.formula("anx_dx_early ~ pred_aam + anx_dx_ch + anx_dx_pre"), family = binomial(link = "logit"), data = fulldata_sensitivity)
logit_model3 <- glm(as.formula("anx_dx_late ~ pred_aam + anx_dx_ch + anx_dx_pre + anx_dx_early"), family = binomial(link = "logit"), data = fulldata_sensitivity)

# Extract coefficient for the predictor from the second stage
estimate_iv1 <- coef(logit_model1)['pred_aam']
estimate_iv2 <- coef(logit_model2)['pred_aam']
estimate_iv3 <- coef(logit_model3)['pred_aam']

# Get robust standard errors accounting for uncertainty in first stage
robust_se_iv1 <- sqrt(diag(vcovHC(logit_model1, type = "HC1"))['pred_aam'])
robust_se_iv2 <- sqrt(diag(vcovHC(logit_model2, type = "HC1"))['pred_aam'])
robust_se_iv3 <- sqrt(diag(vcovHC(logit_model3, type = "HC1"))['pred_aam'])

# Calculate z-value, two-tailed p-value, and confidence interval for the IV estimate
z_iv1 <- estimate_iv1 / robust_se_iv1
pval_iv1 <- 2 * pnorm(abs(z_iv1), lower.tail = FALSE)
ci_low_iv1 <- estimate_iv1 - qnorm(0.975) * robust_se_iv1
ci_high_iv1 <- estimate_iv1 + qnorm(0.975) * robust_se_iv1

z_iv2 <- estimate_iv2 / robust_se_iv2
pval_iv2 <- 2 * pnorm(abs(z_iv2), lower.tail = FALSE)
ci_low_iv2 <- estimate_iv2 - qnorm(0.975) * robust_se_iv2
ci_high_iv2 <- estimate_iv2 + qnorm(0.975) * robust_se_iv2

z_iv3 <- estimate_iv3 / robust_se_iv3
pval_iv3 <- 2 * pnorm(abs(z_iv3), lower.tail = FALSE)
ci_low_iv3 <- estimate_iv3 - qnorm(0.975) * robust_se_iv3
ci_high_iv3 <- estimate_iv3 + qnorm(0.975) * robust_se_iv3

# Combine IV results into data frames
iv_results_pre <- data.frame(
  Estimate_IV = exp(estimate_iv1),
  Robust_SE_IV = robust_se_iv1,
  Z_IV = z_iv1,
  P_Value_IV = pval_iv1,
  CI_Low_IV = exp(ci_low_iv1),
  CI_High_IV = exp(ci_high_iv1),
  Outcome = "Diagnoses age 9-11")

iv_results_early <- data.frame(
  Estimate_IV = exp(estimate_iv2),
  Robust_SE_IV = robust_se_iv2,
  Z_IV = z_iv2,
  P_Value_IV = pval_iv2,
  CI_Low_IV = exp(ci_low_iv2),
  CI_High_IV = exp(ci_high_iv2),
  Outcome = "Diagnoses age 12-14")

iv_results_mid_late <- data.frame(
  Estimate_IV = exp(estimate_iv3),
  Robust_SE_IV = robust_se_iv3,
  Z_IV = z_iv3,
  P_Value_IV = pval_iv3,
  CI_Low_IV = exp(ci_low_iv3),
  CI_High_IV = exp(ci_high_iv3),
  Outcome = "Diagnoses age 15-17")

#print results
print(iv_results_pre)
print(iv_results_early)
print(iv_results_mid_late)

