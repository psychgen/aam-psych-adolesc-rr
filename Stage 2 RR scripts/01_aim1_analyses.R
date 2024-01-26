#01_aim1_analyses.R

#Aim 1: Investigate to what extent earlier age at menarche is associated 
#with elevated rates of adolescent depressive symptoms/depression 

#The purpose of this script is to run the main analyses corresponding to Aim 1. 
#Specifically, we:

# - Run a linear regression model for aam on dep_sx
# - Run a logistic regression for aam on dep_dx

#Read in the data and load packages:
load(file="data/datlist_tr.RData")
library(tidyverse)
library(effectsize)
library(patchwork)
library(miceadds)

##Hypothesis 1a

#initialise a list to store the models
lm_results <- list()

#create list of model formulas
models <- list(model1a_unadj = smfq_dep_c_14c ~ aam_c_14c,
               model1a_adj = smfq_dep_c_14c ~ aam_c_14c + age_8yr + age_14yr + 
                 parent_income_q1 + parent_education_q1 + parent_cohab_18m +
                 parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                 financ_probs_18m + scl_short_m_q1 + scl_full_m_q3 +
                 epds_short_m_6m + parity + bmi_derived_c_8yr + bmi_derived_c_14c,
               model1a_8yr_unadj = smfq_dep_c_14c ~ aam_c_14c + smfq_dep_c_8yr,
               model1a_8yr_adj = smfq_dep_c_14c ~ aam_c_14c + smfq_dep_c_8yr + age_8yr + age_14yr + 
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
betas <- lapply(lm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances <- lapply(lm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1a_unadj <- summary(pool_mi(qhat=betas$model1a_unadj, u=variances$model1a_unadj))
model1a_adj <- summary(pool_mi(qhat=betas$model1a_adj, u=variances$model1a_adj))
model1a_8yr_unadj <- summary(pool_mi(qhat=betas$model1a_8yr_unadj, u=variances$model1a_8yr_unadj))
model1a_8yr_adj <- summary(pool_mi(qhat=betas$model1a_8yr_adj, u=variances$model1a_8yr_adj))

all_res <- list(model1a_unadj, model1a_adj, model1a_8yr_unadj, model1a_8yr_adj)

est_1a <- map(all_res, function(x){
  x[2, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1a <- bind_rows(est_1a) %>%
  mutate(model = c("model1a_unadj", "model1a_adj", "model1a_8yr_unadj", "model1a_8yr_adj")) %>%
  mutate(child_adj = c("Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nadjusted", "Childhood\nadjusted")) %>%
  mutate(cov_adj = c("Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted"))

row.names(est_comb_1a) <- NULL

#compute additional statistics and rename 
#to enable merging with MR results
ests1a <- est_comb_1a %>%
  mutate(Estimate = results,
         z = results / se,
         p_onetailed = ifelse(results<0,p/2,1),
         CI95_low = lci,
         CI95_high = uci,
         ROPE_low = - d_to_r(0.23),
         outcome = "14-year symptoms") %>%
  select(Estimate,se,z,p,p_onetailed,CI95_low,CI95_high,
         ROPE_low,outcome,cov_adj,child_adj)

##test the impact of accounting for clustering on the standard errors

#initialise a list to store the models
lm_results_noclust <- list()

#create function for linear regression without clustering on maternal ID
for(model in names(models)) {
  #initialise a sub-list for the current model
  lm_results_noclust[[model]] <- list()
  
  #loop through each dataset in the list of imputed datasets
  for(i in seq_along(datlist)) {
    lm_results_noclust[[model]][[i]] <- miceadds::lm.cluster(
      data = datlist[[i]],
      formula = models[[model]],
      cluster = datlist[[i]]$ind_id #set individual ID as clustering variable
    )
  }
}

#extract betas and vars
betas_noclust <- lapply(lm_results_noclust, function(models) {
  lapply(models, function(model) coef(model))
})

variances_noclust <- lapply(lm_results_noclust, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1a_unadj_noclust <- summary(pool_mi(qhat=betas_noclust$model1a_unadj, u=variances_noclust$model1a_unadj))
model1a_adj_noclust <- summary(pool_mi(qhat=betas_noclust$model1a_adj, u=variances_noclust$model1a_adj))
model1a_8yr_unadj_noclust <- summary(pool_mi(qhat=betas_noclust$model1a_8yr_unadj, u=variances_noclust$model1a_8yr_unadj))
model1a_8yr_adj_noclust <- summary(pool_mi(qhat=betas_noclust$model1a_8yr_adj, u=variances_noclust$model1a_8yr_adj))

all_res_noclust <- list(model1a_unadj_noclust, model1a_adj_noclust, 
                        model1a_8yr_unadj_noclust, model1a_8yr_adj_noclust)

est_1a_noclust <- map(all_res_noclust, function(x){
  x[2, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1a_noclust <- bind_rows(est_1a_noclust) %>%
  mutate(model = c("model1a_unadj", "model1a_adj", "model1a_8yr_unadj", "model1a_8yr_adj")) %>%
  mutate(child_adj = c("Childhood\nunadjusted", "Childhood\nunadjusted", "Childhood\nadjusted", "Childhood\nadjusted")) %>%
  mutate(cov_adj = c("Covariate\nunadjusted", "Covariate\nadjusted", "Covariate\nunadjusted", "Covariate\nadjusted"))

row.names(est_comb_1a_noclust) <- NULL

#print results without accounting for clustering
print(est_comb_1a_noclust) #se = 0.00833 (adjusted estimate)

#print results with accounting for clustering
print(est_comb_1a) #se = 0.00835 (adjusted estimate)

#differences are negligible and do not influence the interpretation of results

#Hypothesis 1b

#initialise a list to store the models
glm_results <- list()

#create list of model formulas
models <- list(model1b_unadj = dep_dx_ad ~ aam_c_14c,
               model1b_adj = dep_dx_ad ~ aam_c_14c + age_8yr + age_14yr + 
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
betas_1b <- lapply(glm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances_1b <- lapply(glm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1b_unadj <- summary(pool_mi(qhat=betas_1b$model1b_unadj, u=variances_1b$model1b_unadj))
model1b_adj <- summary(pool_mi(qhat=betas_1b$model1b_adj, u=variances_1b$model1b_adj))

all_res_1b <- list(model1b_unadj, model1b_adj)

est_1b <- map(all_res_1b, function(x){
  x[2, ] %>%
    select(results,se,lci=`(lower`,uci=`upper)`,p)
})

est_comb_1b <- bind_rows(est_1b) %>%
  mutate(model = c("model1b_unadj", "model1b_adj")) %>%
  mutate(cov_adj = c("Covariate\nunadjusted", "Covariate\nadjusted"))

row.names(est_comb_1b) <- NULL

#compute additional statistics and rename 
#to enable merging with MR results
ests1b <- est_comb_1b %>%
  mutate(Estimate = exp(results),
         z = results / se,
         p_onetailed = ifelse(results<0,p/2,1),
         CI95_low = exp(lci),
         CI95_high = exp(uci),
         ROPE_low = d_to_oddsratio(-0.23),
         outcome = "Diagnoses age 10-17",
         child_adj = "Childhood\nunadjusted") %>%
  select(Estimate,se,z,p,p_onetailed,CI95_low,CI95_high,
         ROPE_low,outcome,cov_adj,child_adj)

#merge results for symptoms and diagnoses
ests1 <- ests1a %>%
  full_join(ests1b)

#save out
save(ests1, file="./data/aim1_results_tr.RData")

