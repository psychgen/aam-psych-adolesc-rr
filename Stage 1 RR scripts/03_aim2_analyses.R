# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Aim 2: Investigate to what extent effects of age at menarche extend to other domains of mental health (`03_aim2_analyses.R`) 
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("03_aim2_analyses_codemap.png", error = FALSE)
#'

#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to run the primary analyses corresponding to Aim 2. 
#'Specifically, we will:
#'
#'  * Run a linear regression model with NHS+equivalence tests for effect of aam on all sx variables
#'  * Run a logistic regression with NHS+equivalence tests for effect of aam on on all dx variables 
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
#'
#'## Hypotheses 2.1-4a
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
  model2a_unadj <- lm(get(domain[[2]]) ~ aam_c_14s , 
                      data = fulldata) 
  
  
  model2a_adj <- lm(get(domain[[2]]) ~ aam_c_14s +
                      sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                      parent_educ_derived_q1 + parent_cohab_18m +
                      parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                      financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                      epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                      bmi_derived_c_14s,
                    data = fulldata)
  
  model2a_dep_unadj <- lm(get(domain[[2]]) ~ aam_c_14s + smfq_dep_c_14s , 
                          data = fulldata) 
  
  
  model2a_dep_adj <- lm(get(domain[[2]]) ~ aam_c_14s + smfq_dep_c_14s +
                          sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                          parent_educ_derived_q1 + parent_cohab_18m +
                          parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                          financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                          epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                          bmi_derived_c_14s,
                        data = fulldata)
  
  model2a_dep_8yr_unadj <- lm(get(domain[[2]]) ~ aam_c_14s + smfq_dep_c_14s + get(domain[[1]]), 
                              data = fulldata) 
  
  model2a_dep_8yr_adj <- lm(get(domain[[2]]) ~ aam_c_14s + smfq_dep_c_14s + get(domain[[1]]) +
                              sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                              parent_educ_derived_q1 + parent_cohab_18m +
                              parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                              financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                              epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                              bmi_derived_c_14s,
                            data = fulldata)
  
  models <- list(model2a_unadj, model2a_adj,
                 model2a_dep_unadj,model2a_dep_adj,
                 model2a_dep_8yr_unadj,model2a_dep_8yr_adj)
  
  all_res <- map(models, function(res){
    
    #Extract the results of the NHST
    
    nhst_res <- res %>% 
      summary() %>% 
      .$coefficients %>%
      as.data.frame()%>% 
      rownames_to_column("predictor")
    
    #Run the equivalence test (using d=0.22 as SESOI)
    
    equiv_res <- parameters::equivalence_test(
      res, range= c(-effectsize::d_to_r(0.22),effectsize::d_to_r(0.22)), rule="classic")  %>% 
      as.data.frame() %>% 
      `row.names<-`(NULL)   
    
    #Combine 
    
    comb_res <- nhst_res %>% 
      bind_cols(equiv_res) 
  })
  
  domains_res[[names(domains)[[n]]]] <- all_res
  
}

#'
#'## Hypotheses 2.1-3b
#'
# ----
# List the symptom variables
domains <-list(anx=c("anx_dx_0_8","anx_dx_10_17"),
               adhd = c("adhd_dx_0_8","adhd_dx_10_17"),
               cdodd=c("cdodd_dx_0_8","cdodd_dx_10_17"))


# Run models for all domains and save results in a list
domains_res <- list()
n=0
for(domain in domains){
  n=n+1
  model2b_unadj <- glm(get(domain[[2]]) ~ aam_c_14s ,  
                       family = binomial(link="logit"),
                       data = fulldata) 
  
  
  model2b_adj <- glm(get(domain[[2]]) ~ aam_c_14s +
                       sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                       parent_educ_derived_q1 + parent_cohab_18m +
                       parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                       financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                       epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                       bmi_derived_c_14s, 
                     family = binomial(link="logit"),
                     data = fulldata)
  
  model2b_dep_unadj <- glm(get(domain[[2]]) ~ aam_c_14s + dep_dx_10_17 ,  
                           family = binomial(link="logit"),
                           data = fulldata) 
  
  
  model2b_dep_adj <- glm(get(domain[[2]]) ~ aam_c_14s + dep_dx_10_17 +
                           sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                           parent_educ_derived_q1 + parent_cohab_18m +
                           parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                           financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                           epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                           bmi_derived_c_14s, 
                         family = binomial(link="logit"),
                         data = fulldata)
  
  model2b_dep_8yr_unadj <- glm(get(domain[[2]]) ~ aam_c_14s + dep_dx_10_17 + get(domain[[1]]), 
                               family = binomial(link="logit"), 
                               data = fulldata) 
  
  model2b_dep_8yr_adj <- glm(get(domain[[2]]) ~ aam_c_14s + dep_dx_10_17 + get(domain[[1]]) +
                               sex + age_at_q_ret_8yr + parent_income_derived_q1 +
                               parent_educ_derived_q1 + parent_cohab_18m +
                               parent_cohab_3yr + p_age_at_birth + m_age_at_birth +
                               financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 +
                               epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr +
                               bmi_derived_c_14s, 
                             family = binomial(link="logit"),
                             data = fulldata)
  
  models <- list(model2b_unadj, model2b_adj,
                 model2b_dep_unadj,model2b_dep_adj,
                 model2b_dep_8yr_unadj,model2b_dep_8yr_adj)
  
  all_res <- map(models, function(res){
    
    #Extract the results of the NHST
    
    nhst_res <- res %>% 
      summary() %>% 
      .$coefficients %>%
      as.data.frame()%>% 
      rownames_to_column("predictor")
    
    #Run the equivalence test (using d=0.22 as SESOI)
    
    equiv_res <- parameters::equivalence_test(
      res, range= c(-effectsize::d_to_r(0.22),effectsize::d_to_r(0.22)), rule="classic")  %>% 
      as.data.frame() %>% 
      `row.names<-`(NULL)   
    
    #Combine 
    
    comb_res <- nhst_res %>% 
      bind_cols(equiv_res) 
  })
  
  domains_res[[names(domains)[[n]]]] <- all_res
  
}



# ---- eval=F,include=F
#This is the code to render and knit the script - called automatically by 
#generate_code_walkthrough()
ezknitr::ezspin("scripts/template.R", 
                out_dir="scripts/reports",
                keep_rmd = T, 
                verbose=T)
