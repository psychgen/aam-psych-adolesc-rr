#01.1_aim1_exploratory_analyses.R

##This script runs exploratory observational analyses, with a) breast stage
##as an additional exposure, and b) models with a categorised version of the 
##age at menarche exposure and a dichotomised depressive symptom outcome.

#Read in the data and load packages:
load(file="data/datlist_tr.RData")
library(tidyverse)
library(effectsize)
library(patchwork)
library(miceadds)

# Hypothesis 1a (breast stage instead of age at menarche as exposure)

#initialise a list to store the models
lm_results <- list()

#create list of model formulas
models <- list(model1a_unadj_br = smfq_dep_c_14c ~ breast_c_14c,
               model1a_adj_br = smfq_dep_c_14c ~ breast_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c,
               model1a_8yr_unadj_br = smfq_dep_c_14c ~ breast_c_14c + smfq_dep_c_8yr,
               model1a_8yr_adj_br = smfq_dep_c_14c ~ breast_c_14c + smfq_dep_c_8yr + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m + 
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth + financ_probs_18m + 
                 scl_short_m_q1 + scl_full_m_q3 + epds_short_m_6m + parity + 
                 bmi_derived_c_8yr + bmi_derived_c_14c)

#create function for linear regression with clustering on maternal ID
for(model in names(models)) {
  #initialise a sub-list for the current model
  lm_results[[model]] <- list()
  
  #loop through each dataset in the list of imputed datasets
  for(i in seq_along(datlist)) {
    lm_results[[model]][[i]] <- miceadds::lm.cluster(
      data = datlist[[i]],
      formula = models[[model]],
      cluster = datlist[[i]]$m_id
    )
  }
}

#extract betas and vars
betas_br <- lapply(lm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances_br <- lapply(lm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1a_unadj_br <- summary(pool_mi(qhat=betas_br$model1a_unadj, u=variances_br$model1a_unadj))
model1a_adj_br <- summary(pool_mi(qhat=betas_br$model1a_adj, u=variances_br$model1a_adj))
model1a_8yr_unadj_br <- summary(pool_mi(qhat=betas_br$model1a_8yr_unadj, u=variances_br$model1a_8yr_unadj))
model1a_8yr_adj_br <- summary(pool_mi(qhat=betas_br$model1a_8yr_adj, u=variances_br$model1a_8yr_adj))

all_res_br <- list(model1a_unadj_br, model1a_adj_br, model1a_8yr_unadj_br, model1a_8yr_adj_br)

est_1a_br <- map(all_res_br, function(x){
  x[2, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1a_br <- bind_rows(est_1a_br) %>%
  mutate(model = c("model1a_unadj_br", "model1a_adj_br", "model1a_8yr_unadj_br", "model1a_8yr_adj_br")) %>%
  mutate(child_adj = c("Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nadjusted", "Childhood\nadjusted")) %>%
  mutate(cov_adj = c("Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted"))

row.names(est_comb_1a_br) <- NULL

#compute additional statistics and rename 
#to enable merging with MR results
ests1a_br <- est_comb_1a_br %>%
  mutate(Estimate = results,
         z = results / se,
         two_tailed_p = p,
         one_tailed_p = p/2,
         CI_low95 = lci,
         CI_high95 = uci,
         outcome = "14-year symptoms") %>%
  select(Estimate,se,z,two_tailed_p,one_tailed_p,CI_low95,
         CI_high95,outcome,cov_adj,child_adj)

#Hypothesis 1b (breast stage instead of age at menarche as exposure)

#initialise a list to store the models
glm_results <- list()

#create list of model formulas
models <- list(model1b_unadj_br = dep_dx_ad ~ breast_c_14c,
               model1b_adj_br = dep_dx_ad ~ breast_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c)

for(model in names(models)) {
  #initialise a sub-list for the current model
  glm_results[[model]] <- list()
  
  #loop through each dataset in the list of imputed datasets
  for(i in seq_along(datlist)) {
    glm_results[[model]][[i]] <- miceadds::glm.cluster(
      data = datlist[[i]],
      formula = models[[model]],
      family=binomial(link="logit"),
      cluster = datlist[[i]]$m_id
    )
  }
}

#extract betas and vars
betas_1b_br <- lapply(glm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances_1b_br <- lapply(glm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1b_unadj_br <- summary(pool_mi(qhat=betas_1b_br$model1b_unadj, u=variances_1b_br$model1b_unadj))
model1b_adj_br <- summary(pool_mi(qhat=betas_1b_br$model1b_adj, u=variances_1b_br$model1b_adj))

all_res_1b_br <- list(model1b_unadj_br, model1b_adj_br)

est_1b_br <- map(all_res_1b_br, function(x){
  x[2, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1b_br <- bind_rows(est_1b_br) %>%
  mutate(model = c("model1b_unadj_br", "model1b_adj_br")) %>%
  mutate(cov_adj = c("Covariate\nunadjusted", "Covariate\nadjusted"))

row.names(est_comb_1b_br) <- NULL

#compute additional test statistics
ests1b_br <- est_comb_1b_br %>%
  mutate(Estimate = exp(results),
         z = results / se,
         two_tailed_p = p,
         one_tailed_p = p/2,
         CI_low95 = exp(lci),
         CI_high95 = exp(uci),
         outcome = "Diagnoses age 10-17",
         child_adj = "Childhood\nunadjusted") %>%
  select(Estimate,se,z,two_tailed_p,one_tailed_p,CI_low95,
         CI_high95,outcome,cov_adj,child_adj)

#merge
ests1_br <- ests1a_br %>%
  full_join(ests1b_br)

# Hypothesis 1a (age at menarche and breast stage both included as exposures)

#initialise a list to store the models
lm_results <- list()

#create list of model formulas
models <- list(model1a_unadj_br2 = smfq_dep_c_14c ~ aam_c_14c + breast_c_14c,
               model1a_adj_br2 = smfq_dep_c_14c ~ aam_c_14c + breast_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c,
               model1a_8yr_unadj_br2 = smfq_dep_c_14c ~ aam_c_14c + breast_c_14c + smfq_dep_c_8yr,
               model1a_8yr_adj_br2 = smfq_dep_c_14c ~ aam_c_14c + breast_c_14c + smfq_dep_c_8yr + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m + 
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth + financ_probs_18m + 
                 scl_short_m_q1 + scl_full_m_q3 + epds_short_m_6m + parity + 
                 bmi_derived_c_8yr + bmi_derived_c_14c)

#create function for linear regression with clustering on maternal ID
for(model in names(models)) {
  #initialise a sub-list for the current model
  lm_results[[model]] <- list()
  
  #loop through each dataset in the list of imputed datasets
  for(i in seq_along(datlist)) {
    lm_results[[model]][[i]] <- miceadds::lm.cluster(
      data = datlist[[i]],
      formula = models[[model]],
      cluster = datlist[[i]]$m_id
    )
  }
}

#extract betas and vars
betas_br2 <- lapply(lm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances_br2 <- lapply(lm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1a_unadj_br2 <- summary(pool_mi(qhat=betas_br2$model1a_unadj, u=variances_br2$model1a_unadj))
model1a_adj_br2 <- summary(pool_mi(qhat=betas_br2$model1a_adj, u=variances_br2$model1a_adj))
model1a_8yr_unadj_br2 <- summary(pool_mi(qhat=betas_br2$model1a_8yr_unadj, u=variances_br2$model1a_8yr_unadj))
model1a_8yr_adj_br2 <- summary(pool_mi(qhat=betas_br2$model1a_8yr_adj, u=variances_br2$model1a_8yr_adj))

all_res_br2 <- list(model1a_unadj_br2, model1a_adj_br2, model1a_8yr_unadj_br2, model1a_8yr_adj_br2)

est_1a_br2 <- map(all_res_br2, function(x){
  x[2:3, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1a_br2 <- bind_rows(est_1a_br2) %>%
  mutate(predictor = c("Age at menarche","Breast stage","Age at menarche","Breast stage","Age at menarche","Breast stage","Age at menarche","Breast stage")) %>%
  mutate(model = c(rep("model1a_unadj_br2",2), rep("model1a_adj_br2",2), rep("model1a_8yr_unadj_br2",2), rep("model1a_8yr_adj_br2",2))) %>%
  mutate(child_adj = c(rep("Childhood\nunadjusted",2), rep("Childhood\nunadjusted",2), rep("Childhood\nadjusted",2), rep("Childhood\nadjusted",2))) %>%
  mutate(cov_adj = c(rep("Covariate\nunadjusted",2), rep("Covariate\nadjusted",2), rep("Covariate\nunadjusted",2), rep("Covariate\nadjusted",2)))

row.names(est_comb_1a_br2) <- NULL

#compute additional statistics and rename 
#to enable merging with MR results
ests1a_br2 <- est_comb_1a_br2 %>%
  mutate(Estimate = results,
         z = results / se,
         two_tailed_p = p,
         one_tailed_p = p/2,
         CI_low95 = lci,
         CI_high95 = uci,
         outcome = "14-year symptoms") %>%
  select(predictor,Estimate,se,z,two_tailed_p,one_tailed_p,
         CI_low95,CI_high95,outcome,cov_adj,child_adj)

#Hypothesis 1b (age at menarche and breast stage both included as exposures)

#initialise a list to store the models
glm_results <- list()

#create list of model formulas
models <- list(model1b_unadj_br2 = dep_dx_ad ~ aam_c_14c + breast_c_14c,
               model1b_adj_br2 = dep_dx_ad ~ aam_c_14c + breast_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c)

for(model in names(models)) {
  #initialise a sub-list for the current model
  glm_results[[model]] <- list()
  
  #loop through each dataset in the list of imputed datasets
  for(i in seq_along(datlist)) {
    glm_results[[model]][[i]] <- miceadds::glm.cluster(
      data = datlist[[i]],
      formula = models[[model]],
      family=binomial(link="logit"),
      cluster = datlist[[i]]$m_id
    )
  }
}

#extract betas and vars
betas_1b_br2 <- lapply(glm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances_1b_br2 <- lapply(glm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1b_unadj_br2 <- summary(pool_mi(qhat=betas_1b_br2$model1b_unadj, u=variances_1b_br2$model1b_unadj))
model1b_adj_br2 <- summary(pool_mi(qhat=betas_1b_br2$model1b_adj, u=variances_1b_br2$model1b_adj))

all_res_1b_br2 <- list(model1b_unadj_br2, model1b_adj_br2)

est_1b_br2 <- map(all_res_1b_br2, function(x){
  x[2:3, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1b_br2 <- bind_rows(est_1b_br2) %>%
  mutate(predictor = c("Age at menarche","Breast stage","Age at menarche","Breast stage")) %>%
  mutate(model = c(rep("model1b_unadj_br2",2), rep("model1b_adj_br2",2))) %>%
  mutate(cov_adj = c(rep("Covariate\nunadjusted",2), rep("Covariate\nadjusted",2)))

row.names(est_comb_1b_br2) <- NULL

#compute additional test statistics 
ests1b_br2 <- est_comb_1b_br2 %>%
  mutate(Estimate = exp(results),
         z = results / se,
         two_tailed_p = p,
         one_tailed_p = p/2,
         CI_low95 = exp(lci),
         CI_high95 = exp(uci),
         outcome = "Diagnoses age 10-17",
         child_adj = "Childhood\nunadjusted") %>%
  select(predictor,Estimate,se,z,two_tailed_p,one_tailed_p,
         CI_low95,CI_high95,outcome,cov_adj,child_adj)

#merge results for symptoms and diagnoses
ests1_br2 <- ests1a_br2 %>%
  full_join(ests1b_br2)

#use 2 decimal points for reporting in text
ests1_br_2dp <- ests1_br
ests1_br_2dp[] <- lapply(ests1_br_2dp, function(x) if(is.numeric(x)) round(x, 2) else x)
ests1_br_2dp

ests1_br2_2dp <- ests1_br2
ests1_br2_2dp[] <- lapply(ests1_br2_2dp, function(x) if(is.numeric(x)) round(x, 2) else x)
ests1_br2_2dp

#save out results for breast stage and aam+breast stage as exposures
save(ests1_br, file="./data/aim1_breast_stage_results.RData")
save(ests1_br2, file="./data/aim1_aam_and_breast_stage_results.RData")

##Hypothesis 1a (categorised aam as exposure and dichotomised smfq as outcome)

#load data
load(file="data/datlist_categorised.RData")

#set the factor levels in each dataset and create contrast variables
ordered_levels <- c("early", "average", "late")

#use lapply to set the factor levels in each dataset and create the dichotomous variables
datlist_cat <- lapply(datlist_cat, function(df) {
  
  #ensure 'aam_cat' is an ordered factor with the specified levels
  df$aam_cat <- factor(df$aam_cat, levels = ordered_levels, ordered = TRUE)
  
  #set the contrasts for aam_cat, using 'average' as the reference level
  contrasts(df$aam_cat) <- contr.treatment(3, base = 2)
  
  #make dichotomised smfq outcomes numeric
  df$smfq_14yr_dic <- as.numeric(df$smfq_14yr_dic)-1
  df$smfq_8yr_dic <- as.numeric(df$smfq_8yr_dic)-1
  
  return(df)
  
})

#initialise a list to store the models
glm_results <- list()

#create list of model formulas
models <- list(model1a_unadj_cat = smfq_14yr_dic ~ aam_cat,
               model1a_adj_cat = smfq_14yr_dic ~ aam_cat + 
                 age_8yr + age_14yr + parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c,
               model1a_8yr_unadj_cat = smfq_14yr_dic ~ aam_cat + smfq_8yr_dic,
               model1a_8yr_adj_cat = smfq_14yr_dic ~ aam_cat + smfq_8yr_dic + 
                 age_8yr + age_14yr + parent_income_q1 + parent_education_q1 + parent_cohab_18m + 
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth + financ_probs_18m + 
                 scl_short_m_q1 + scl_full_m_q3 + epds_short_m_6m + parity + 
                 bmi_derived_c_8yr + bmi_derived_c_14c)

#create function for logistic regression with clustering on maternal ID
for(model in names(models)) {
  #initialise a sub-list for the current model
  glm_results[[model]] <- list()
  
  #loop through each dataset in the list of imputed datasets
  for(i in seq_along(datlist_cat)) {
    glm_results[[model]][[i]] <- miceadds::glm.cluster(
      data = datlist_cat[[i]],
      formula = models[[model]],
      family=binomial(link="logit"),
      cluster = datlist_cat[[i]]$m_id
    )
  }
}

#extract betas and vars
betas_cat <- lapply(glm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances_cat <- lapply(glm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1a_unadj_cat <- summary(pool_mi(qhat=betas_cat$model1a_unadj_cat, u=variances_cat$model1a_unadj_cat))
model1a_adj_cat <- summary(pool_mi(qhat=betas_cat$model1a_adj_cat, u=variances_cat$model1a_adj_cat))
model1a_8yr_unadj_cat <- summary(pool_mi(qhat=betas_cat$model1a_8yr_unadj_cat, u=variances_cat$model1a_8yr_unadj_cat))
model1a_8yr_adj_cat <- summary(pool_mi(qhat=betas_cat$model1a_8yr_adj_cat, u=variances_cat$model1a_8yr_adj_cat))

all_res_cat <- list(model1a_unadj_cat, model1a_adj_cat, model1a_8yr_unadj_cat, model1a_8yr_adj_cat)

est_1a_cat <- map(all_res_cat, function(x){
  x[2:3, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1a_cat <- bind_rows(est_1a_cat) %>%
  mutate(predictor = c("aam_early","aam_late","aam_early","aam_late","aam_early","aam_late","aam_early","aam_late")) %>%
  mutate(model = c(rep("model1a_unadj_cat",2), rep("model1a_adj_cat",2), rep("model1a_8yr_unadj_cat",2), rep("model1a_8yr_adj_cat",2))) %>%
  mutate(child_adj = c(rep("Childhood\nunadjusted",2), rep("Childhood\nunadjusted",2), rep("Childhood\nadjusted",2), rep("Childhood\nadjusted",2))) %>%
  mutate(cov_adj = c(rep("Covariate\nunadjusted",2), rep("Covariate\nadjusted",2), rep("Covariate\nunadjusted",2), rep("Covariate\nadjusted",2)))

row.names(est_comb_1a_cat) <- NULL

#compute additional statistics and rename 
#to enable merging with MR results
ests1a_cat <- est_comb_1a_cat %>%
  mutate(Estimate = exp(results),
         z = results / se,
         two_tailed_p = p,
         one_tailed_p = p/2,
         CI_low95 = exp(lci),
         CI_high95 = exp(uci),
         outcome = "14-year symptoms\ndichotomised") %>%
  select(predictor,Estimate,se,z,two_tailed_p,one_tailed_p,
         CI_low95,CI_high95,outcome,cov_adj,child_adj)

#use 2 decimal points for reporting in text
ests1a_cat_2dp <- ests1a_cat
ests1a_cat_2dp[] <- lapply(ests1a_cat_2dp, function(x) if(is.numeric(x)) round(x, 2) else x)
ests1a_cat_2dp

#Hypothesis 1b (categorised aam as exposure and depression dx as outcome)

#initialise a list to store the models
glm_results <- list()

#create list of model formulas
models <- list(model1b_unadj_cat = dep_dx_ad ~ aam_cat,
               model1b_adj_cat = dep_dx_ad ~ aam_cat + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c)

for(model in names(models)) {
  #initialise a sub-list for the current model
  glm_results[[model]] <- list()
  
  #loop through each dataset in the list of imputed datasets
  for(i in seq_along(datlist_cat)) {
    glm_results[[model]][[i]] <- miceadds::glm.cluster(
      data = datlist_cat[[i]],
      formula = models[[model]],
      family=binomial(link="logit"),
      cluster = datlist_cat[[i]]$m_id
    )
  }
}

#extract betas and vars
betas_1b_cat <- lapply(glm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances_1b_cat <- lapply(glm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1b_unadj_cat <- summary(pool_mi(qhat=betas_1b_cat$model1b_unadj_cat, u=variances_1b_cat$model1b_unadj_cat))
model1b_adj_cat <- summary(pool_mi(qhat=betas_1b_cat$model1b_adj_cat, u=variances_1b_cat$model1b_adj_cat))

all_res_1b_cat <- list(model1b_unadj_cat, model1b_adj_cat)

est_1b_cat <- map(all_res_1b_cat, function(x){
  x[2:3, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1b_cat <- bind_rows(est_1b_cat) %>%
  mutate(predictor = c("aam_early", "aam_late", "aam_early", "aam_late")) %>%
  mutate(model = c(rep("model1b_unadj_cat",2), rep("model1b_adj_cat",2))) %>%
  mutate(cov_adj = c(rep("Covariate\nunadjusted",2), rep("Covariate\nadjusted",2)))

row.names(est_comb_1b_cat) <- NULL

#compute additional test statistics
ests1b_cat <- est_comb_1b_cat %>%
  mutate(Estimate = exp(results),
         z = results / se,
         two_tailed_p = p,
         one_tailed_p = p/2,
         CI_low95 = exp(lci),
         CI_high95 = exp(uci),
         outcome = "Diagnoses age 10-17",
         child_adj = "Childhood\nunadjusted") %>%
  select(predictor,Estimate,se,z,two_tailed_p,one_tailed_p,
         CI_low95,CI_high95,outcome,cov_adj,child_adj)

#use 2 decimal points for reporting in text
ests1b_cat_2dp <- ests1b_cat
ests1b_cat_2dp[] <- lapply(ests1b_cat_2dp, function(x) if(is.numeric(x)) round(x, 2) else x)
ests1b_cat_2dp
