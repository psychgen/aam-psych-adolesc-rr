# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Aim 4: Determine whether any associations between early age at menarche and adolescent symptoms/diagnoses in other domains are likely to be causal (`05_aim4_analyses.R`)
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("05_aim4_analyses_codemap.png", error = FALSE)
#'
#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to run the primary analyses corresponding to Aim 4. 
#'Specifically, we will:
#'
#'  * Run a 2SLS one-sample MR  with NHS+equivalence tests for causal effect of aam on sx 
#'  * Run a 2SLS one-sample MR  with NHST+equivalence tests for causal effect of aam on earlier sx (negative control)
#'  * Run a 2SLS (logistic 2nd stage) one-sample MR with NHS+equivalence tests for causal effect of aam on on all dx variables 
#'  
#'
#'&nbsp;
#'&nbsp;
#'
#'Read in the data and load packages:
# ---- 
load(file="data/analysis_dataset.RData")
library(tidyverse)
library(parameters)
library(effectsize)
library(AER)
library(modelr)
library(broom)
library(sandwich)
#'
#'## Hypotheses 4.1-4a
#'
# ----

# List the symptom variables
domains <-list(anx=c("scared_anx_c_8yr","scared_anx_c_14s"),
               adhd = c("rsdbd_adhd_c_8yr","rsdbd_adhd_c_14m"),
               cd=c("rsdbd_cd_c_8yr","rsdbd_cd_c_14s"),
               odd=c("rsdbd_odd_c_8yr","rsdbd_odd_c_14m"))

# Scale for standardised effects
fulldata <- fulldata %>% 
  mutate(across( c("smfq_dep_c_14s","aam_c_14s","smfq_dep_c_8yr",unlist(domains) ), scale ))

# Run models for all domains and save results in a list
domains_res <- list()
n=0
for(domain in domains){
  n=n+1
  
  model4a_unadj <- ivreg(get(domain[[2]]) ~ aam_c_14s  | PGS_AaM, data = fulldata) %>% 
    broom::tidy(conf.int=T, conf.level=0.9) %>% 
    mutate(ROPE_high= d_to_r(0.2),
           ROPE_low =- d_to_r(0.2)) 
  
  model4a_adj <- ivreg(get(domain[[2]]) ~ aam_c_14s  +
                         sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                         parent_educ_derived_q1 + parent_cohab_18m +
                         parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                         financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                         epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                         bmi_derived_c_14s  | PGS_AaM  +
                         sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                         parent_educ_derived_q1 + parent_cohab_18m +
                         parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                         financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                         epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                         bmi_derived_c_14s, data = fulldata) %>% 
    broom::tidy(conf.int=T, conf.level=0.9) %>% 
    mutate(ROPE_high= d_to_r(0.2),
           ROPE_low =- d_to_r(0.2)) 
  
  model4a_neg_unadj <- ivreg(get(domain[[1]]) ~ aam_c_14s  | PGS_AaM, data = fulldata) %>% 
    broom::tidy(conf.int=T, conf.level=0.9) %>% 
    mutate(ROPE_high= d_to_r(0.2),
           ROPE_low =- d_to_r(0.2)) 
  
  model4a_neg_adj <- ivreg(get(domain[[1]]) ~ aam_c_14s  +
                             sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                             parent_educ_derived_q1 + parent_cohab_18m +
                             parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                             financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                             epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                             bmi_derived_c_14s  | PGS_AaM  +
                             sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                             parent_educ_derived_q1 + parent_cohab_18m +
                             parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                             financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                             epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                             bmi_derived_c_14s, data = fulldata) %>% 
    broom::tidy(conf.int=T, conf.level=0.9) %>% 
    mutate(ROPE_high= d_to_r(0.2),
           ROPE_low =- d_to_r(0.2)) 
  
  all_res <-list (model4a_unadj,model4a_adj,model4a_neg_unadj,model4a_neg_adj)
  
  domains_res[[names(domains)[[n]]]] <- all_res
  
}

#'
#'## Hypotheses 4.1-3b
#'
# ----

# Get PGS predicted values of AaM from first stage linear models

model4b_stage1_unadj <-lm(aam_c_14s ~ PGS_AaM, data = fulldata)

model4b_stage1_adj <-lm(aam_c_14s ~ PGS_AaM +
                          sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                          parent_educ_derived_q1 + parent_cohab_18m +
                          parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                          financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                          epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                          bmi_derived_c_14s  , data = fulldata)

fulldata <- fulldata%>% 
  add_predictions(model3b_stage1_unadj, var= "pred_unadj") %>% 
  add_predictions(model3b_stage1_adj, var= "pred_adj")

# Run logistic second stage

# List the symptom variables
domains <-list(anx=c("anx_dx_0_8","anx_dx_10_17"),
               adhd = c("adhd_dx_0_8","adhd_dx_10_17"),
               cdodd=c("cdodd_dx_0_8","cdodd_dx_10_17"))


# Run models for all domains and save results in a list
domains_res <- list()
n=0
for(domain in domains){
  n=n+1
  
  model4b_stage2_unadj <- glm(get(domain[[2]]) ~ pred_unadj, 
                              family = binomial(link="logit"), 
                              data=fulldata)
  
  model4b_stage2_adj <- glm(get(domain[[2]]) ~ pred_adj  +
                              sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                              parent_educ_derived_q1 + parent_cohab_18m +
                              parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                              financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                              epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                              bmi_derived_c_14s , 
                            family = binomial(link="logit"), 
                            data=fulldata)
  
  model4b_stage2_8yr_unadj <- glm(get(domain[[1]]) ~ pred_unadj, 
                                  family = binomial(link="logit"), 
                                  data=fulldata)
  
  model4b_stage2_8yr_adj <- glm(get(domain[[1]]) ~ pred_adj  +
                                  sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                                  parent_educ_derived_q1 + parent_cohab_18m +
                                  parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                                  financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                                  epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                                  bmi_derived_c_14s , 
                                family = binomial(link="logit"), 
                                data=fulldata)
  # We also need to adjust the SEs to account for the uncertainty in the first stage:
  # Calculate robust standard errors (from https://stats.stackexchange.com/
  #questions/89999/how-to-replicate-statas-robust-binomial-glm-for-proportion-data-in-r)
  adjSEs <- function(x){
    cov.m1 <- vcovHC(x, type = "HC1") 
    ## This is used as it mirrors statas ivreg2 robust, used in Sequeira et al)
    std.err <- sqrt(diag(cov.m1))
    q.val <- qnorm(0.95)
    comb_res <- cbind(
      Estimate = coef(x)
      , "Robust SE" = std.err
      , z = (coef(x)/std.err)
      , "p"= 2 * pnorm(abs(coef(x)/std.err), lower.tail = FALSE)
      , CI_low = coef(x) - q.val  * std.err
      , CI_high = coef(x) + q.val  * std.err
    ) %>% 
      as.data.frame() %>% 
      mutate(ROPE_low =- d_to_oddsratio(0.25, log=T))
    return(comb_res)
  }
  
  all_res= list  (adjSEs(model4b_stage2_unadj),
                  adjSEs(model4b_stage2_adj),
                  adjSEs(model4b_stage2_8yr_unadj),
                  adjSEs(model4b_stage2_8yrs_adj))
  domains_res[[names(domains)[[n]]]] <- all_res
  
}

# ---- eval=F,include=F
#This is the code to render and knit the script - called automatically by 
#generate_code_walkthrough()
ezknitr::ezspin("scripts/template.R", 
                out_dir="scripts/reports",
                keep_rmd = T, 
                verbose=T)

NA
NA
NA
