#02_aim2_analyses.R

#Aim 2: Investigate to what extent effects of age at menarche extend 
#to other domains of mental health

#The purpose of this script is to run primary analyses corresponding to Aim 2. 
#Specifically, we:

# - Run a linear regression model for effect of aam on all sx variables
# - Run a logistic regression for effect of aam on all dx variables 

#Read in the data and load packages:
load(file="data/datlist_tr.RData")
library(dplyr)
library(tibble)
library(patchwork)
library(effectsize)
library(mice)

##Hypotheses 2.1-4a

#list the symptom variables
domains <-list(anx=c("scared_anx_c_8yr","scared_anx_c_14c"),
               adhd = c("rsdbd_adhd_c_8yr","rsdbd_adhd_c_14m"),
               cd=c("rsdbd_cd_c_8yr","rsdbd_cd_c_14c"),
               odd=c("rsdbd_odd_c_8yr","rsdbd_odd_c_14m"))

#create function for linear regression with clustering on maternal ID
cluster_lm <- function(data, formula) {
  data$wgt__ <- 1 #this is set to 1 to bypass a bug in the lm.cluster function
  miceadds::lm.cluster(data=data, formula=formula, cluster="m_id")
}

#run models for all domains and save results in a list
domains_res <- list()
n = 0

for(domain in domains) {
  n = n + 1
  
  #specify completely unadjusted model
  formula_str <- paste(domain[[2]], "~ aam_c_14c")
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  lm_results <- lapply(datlist, FUN=function(data) cluster_lm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(lm_results, FUN=function(rr) coef(rr))
  vars <- lapply(lm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2a_unadj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify covariate adjusted model
  formula_str <- paste(domain[[2]], "~ aam_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c")
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  lm_results <- lapply(datlist, FUN=function(data) cluster_lm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(lm_results, FUN=function(rr) coef(rr))
  vars <- lapply(lm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2a_adj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify depression adjusted and covariate unadjusted model
  formula_str <- paste(domain[[2]], "~ aam_c_14c + smfq_dep_c_14c")
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  lm_results <- lapply(datlist, FUN=function(data) cluster_lm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(lm_results, FUN=function(rr) coef(rr))
  vars <- lapply(lm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2a_dep_unadj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify childhood symptoms adjusted model, not adjusted for covs or depression
  formula_str <- paste(domain[[2]], "~ aam_c_14c +", domain[[1]])
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  lm_results <- lapply(datlist, FUN=function(data) cluster_lm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(lm_results, FUN=function(rr) coef(rr))
  vars <- lapply(lm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2a_8yr_unadj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify childhood symptoms and cov adjusted model, not adjusted for depression
  formula_str <- paste(domain[[2]], "~ aam_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m + 
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth + financ_probs_18m + 
                 scl_short_m_q1 + scl_full_m_q3 + epds_short_m_6m + parity + 
                 bmi_derived_c_8yr + bmi_derived_c_14c +", domain[[1]])
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  lm_results <- lapply(datlist, FUN=function(data) cluster_lm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(lm_results, FUN=function(rr) coef(rr))
  vars <- lapply(lm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2a_8yr_adj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify depression and covariate adjusted model
  formula_str <- paste(domain[[2]], "~ aam_c_14c + smfq_dep_c_14c + age_8yr + 
                 age_14yr + parent_income_q1 + parent_education_q1 + 
                 parent_cohab_18m + parent_cohab_3yr + p_age_at_birth + 
                 m_age_at_birth + financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c")
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  lm_results <- lapply(datlist, FUN=function(data) cluster_lm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(lm_results, FUN=function(rr) coef(rr))
  vars <- lapply(lm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2a_dep_adj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify depression and childhood symptoms adjusted model without covariates
  formula_str <- paste(domain[[2]], "~ aam_c_14c + smfq_dep_c_14c +", domain[[1]])
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  lm_results <- lapply(datlist, FUN=function(data) cluster_lm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(lm_results, FUN=function(rr) coef(rr))
  vars <- lapply(lm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2a_dep_8yr_unadj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify depression and childhood symptoms adjusted model with covariates
  formula_str <- paste(domain[[2]], "~ aam_c_14c + smfq_dep_c_14c + age_8yr + 
                 age_14yr + parent_income_q1 + parent_education_q1 + parent_cohab_18m + 
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth + financ_probs_18m + 
                 scl_short_m_q1 + scl_full_m_q3 + epds_short_m_6m + parity + 
                 bmi_derived_c_8yr + bmi_derived_c_14c +", domain[[1]])
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  lm_results <- lapply(datlist, FUN=function(data) cluster_lm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(lm_results, FUN=function(rr) coef(rr))
  vars <- lapply(lm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2a_dep_8yr_adj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  models <- list(model2a_unadj, model2a_adj,
                 model2a_dep_unadj, model2a_dep_adj, 
                 model2a_8yr_unadj, model2a_8yr_adj,
                 model2a_dep_8yr_unadj, model2a_dep_8yr_adj)
  
  domains_res[[names(domains)[[n]]]] <- models
  
}

est <- map(domains_res, function(x){
  map(x, function(y){
    y[2, ] %>%
      select(results,se,lci=`(lower`,uci=`upper)`,p)
  })
})

est_comb <- map(est, function(z){
  reduce(z, bind_rows) %>%
    mutate(model = c("model2a_unadj", "model2a_adj", "model2a_dep_unadj", "model2a_dep_adj", "model2a_8yr_unadj", "model2a_8yr_adj", "model2a_dep_8yr_unadj", "model2a_dep_8yr_adj")) %>%
    mutate(child_adj = c("Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nadjusted", "Childhood\nadjusted", "Childhood\nadjusted", "Childhood\nadjusted")) %>%
    mutate(cov_adj = c("Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted")) %>%
    mutate(dep_adj = c("Depression\nunadjusted", "Depression\nunadjusted", "Depression\nadjusted", "Depression\nadjusted", "Depression\nunadjusted", "Depression\nunadjusted", "Depression\nadjusted", "Depression\nadjusted"))
})

anx <- est_comb$anx %>%
  mutate(domain = "ANX")
adhd <- est_comb$adhd %>%
  mutate(domain = "ADHD")
cd <- est_comb$cd %>%
  mutate(domain = "CD")
odd <- est_comb$odd %>%
  mutate(domain = "ODD")

ests2a <- anx %>%
  full_join(adhd) %>% 
  full_join(cd) %>% 
  full_join(odd) %>%
  rename(Estimate = results,
         CI95_low = lci,
         CI95_high = uci) %>%
  mutate(ROPE_low =- d_to_r(0.22),
         ROPE_high = d_to_r(0.22),
         outcome = "14-year symptoms")


##Hypotheses 2.1-3b

#list the diagnosis variables
domains <-list(anx=c("anx_dx_ch","anx_dx_ad"),
               adhd = c("adhd_dx_ch","adhd_dx_ad"),
               beh=c("beh_dx_ch","beh_dx_ad"))

#create function for logistic regression with clustering on maternal ID
cluster_glm <- function(data, formula) {
  data$wgt__ <- 1 #this is set to 1 to bypass a bug in the glm.cluster function
  miceadds::glm.cluster(data=data, formula=formula, family=binomial(link="logit"), cluster="m_id")
}

#run models for all domains and save results in a list
domains_res_2b <- list()
n = 0

for(domain in domains) {
  n = n + 1
  
  #specify completely unadjusted model
  formula_str <- paste(domain[[2]], "~ aam_c_14c")
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  glm_results <- lapply(datlist, FUN=function(data) cluster_glm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(glm_results, FUN=function(rr) coef(rr))
  vars <- lapply(glm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2b_unadj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify covariate adjusted model
  formula_str <- paste(domain[[2]], "~ aam_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c")
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  glm_results <- lapply(datlist, FUN=function(data) cluster_glm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(glm_results, FUN=function(rr) coef(rr))
  vars <- lapply(glm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2b_adj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify depression adjusted and covariate unadjusted model
  formula_str <- paste(domain[[2]], "~ aam_c_14c + dep_dx_ad")
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  glm_results <- lapply(datlist, FUN=function(data) cluster_glm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(glm_results, FUN=function(rr) coef(rr))
  vars <- lapply(glm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2b_dep_unadj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify childhood dx adjusted model, not adjusted for covs or depression
  formula_str <- paste(domain[[2]], "~ aam_c_14c +", domain[[1]])
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  glm_results <- lapply(datlist, FUN=function(data) cluster_glm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(glm_results, FUN=function(rr) coef(rr))
  vars <- lapply(glm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2b_8yr_unadj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify childhood dx and cov adjusted model, not adjusted for depression
  formula_str <- paste(domain[[2]], "~ aam_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m + 
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth + financ_probs_18m + 
                 scl_short_m_q1 + scl_full_m_q3 + epds_short_m_6m + parity + 
                 bmi_derived_c_8yr + bmi_derived_c_14c +", domain[[1]])
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  glm_results <- lapply(datlist, FUN=function(data) cluster_glm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(glm_results, FUN=function(rr) coef(rr))
  vars <- lapply(glm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2b_8yr_adj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify depression and covariate adjusted model
  formula_str <- paste(domain[[2]], "~ aam_c_14c + dep_dx_ad + age_8yr + 
                 age_14yr + parent_income_q1 + parent_education_q1 + 
                 parent_cohab_18m + parent_cohab_3yr + p_age_at_birth + 
                 m_age_at_birth + financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c")
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  glm_results <- lapply(datlist, FUN=function(data) cluster_glm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(glm_results, FUN=function(rr) coef(rr))
  vars <- lapply(glm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2b_dep_adj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify depression and childhood dx adjusted model without covariates
  formula_str <- paste(domain[[2]], "~ aam_c_14c + dep_dx_ad +", domain[[1]])
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  glm_results <- lapply(datlist, FUN=function(data) cluster_glm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(glm_results, FUN=function(rr) coef(rr))
  vars <- lapply(glm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2b_dep_8yr_unadj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  #specify depression and childhood dx adjusted model with covariates
  formula_str <- paste(domain[[2]], "~ aam_c_14c + dep_dx_ad + age_8yr + 
                 age_14yr + parent_income_q1 + parent_education_q1 + parent_cohab_18m + 
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth + financ_probs_18m + 
                 scl_short_m_q1 + scl_full_m_q3 + epds_short_m_6m + parity + 
                 bmi_derived_c_8yr + bmi_derived_c_14c +", domain[[1]])
  model_formula <- as.formula(formula_str)
  
  #run the regression on each imputed dataset
  glm_results <- lapply(datlist, FUN=function(data) cluster_glm(data, model_formula))
  
  #extract betas and vars
  betas <- lapply(glm_results, FUN=function(rr) coef(rr))
  vars <- lapply(glm_results, FUN=function(rr) vcov(rr))
  
  #pool the results
  model2b_dep_8yr_adj <- summary(miceadds::pool_mi(qhat=betas, u=vars))
  
  models <- list(model2b_unadj, model2b_adj,
                 model2b_dep_unadj, model2b_dep_adj, 
                 model2b_8yr_unadj, model2b_8yr_adj,
                 model2b_dep_8yr_unadj, model2b_dep_8yr_adj)
  
  domains_res_2b[[names(domains)[[n]]]] <- models
  
}

est_2b <- map(domains_res_2b, function(x){
  map(x, function(y){
    y[2, ] %>%
      select(results,se,lci=`(lower`,uci=`upper)`,p)
  })
})

est_comb_2b <- map(est_2b, function(z){
  reduce(z, bind_rows) %>%
    mutate(model = c("model2b_unadj", "model2b_adj", "model2b_dep_unadj", "model2b_dep_adj", "model2b_8yr_unadj", "model2b_8yr_adj", "model2b_dep_8yr_unadj", "model2b_dep_8yr_adj")) %>%
    mutate(child_adj = c("Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nadjusted", "Childhood\nadjusted", "Childhood\nadjusted", "Childhood\nadjusted")) %>%
    mutate(cov_adj = c("Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted")) %>%
    mutate(dep_adj = c("Depression\nunadjusted", "Depression\nunadjusted", "Depression\nadjusted", "Depression\nadjusted", "Depression\nunadjusted", "Depression\nunadjusted", "Depression\nadjusted", "Depression\nadjusted"))
})

anx_dx <- est_comb_2b$anx %>%
  mutate(domain = "ANX")
adhd_dx <- est_comb_2b$adhd %>%
  mutate(domain = "ADHD")
dbd_dx <- est_comb_2b$beh %>%
  mutate(domain = "DBD")

ests2b <- anx_dx %>%
  full_join(adhd_dx) %>% 
  full_join(dbd_dx) %>% 
  rename(Estimate = results,
         CI_low = lci,
         CI_high = uci) %>%
  mutate(Estimate = exp(Estimate),
         CI_low = exp(CI_low),
         CI_high = exp(CI_high),
         ROPE_low = d_to_oddsratio(-0.22),
         ROPE_high = d_to_oddsratio(0.22),
         outcome = "Diagnoses age 10-17")

#merge symptom and diagnostic outcome results
ests2 <- ests2a %>%
  full_join(ests2b)

#save out
save(ests2, file="./data/aim2_results.RData")
