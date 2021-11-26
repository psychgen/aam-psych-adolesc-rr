# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Aim 1: Investigate to what extent early age at menarche is associated with elevated rates of adolescent depressive symptoms/depression (`02_aim1_analyses.R`) 
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("02_aim1_analyses_codemap.png", error = FALSE)
#'

#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to run the main analyses corresponding to Aim 1. 
#'Specifically, we will:
#'
#'  * Run a linear regression model with NHS+equivalence tests for aam on dep_sx
#'  * Run a logistic regression with NHS+equivalence tests for aam on dep_DX 
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
#'## Hypothesis 1a
# ----

fulldata <- fulldata %>% 
  mutate(across( c("smfq_dep_c_14s","aam_c_14s","smfq_dep_c_8yr" ), scale ))

model1a_unadj <- lm(smfq_dep_c_14s ~ aam_c_14s , 
                    data = fulldata) 

model1a_adj <- lm(smfq_dep_c_14s ~ aam_c_14s + 
                    sex + age_at_q_ret_8yr + parent_income_derived_q1 + 
                    parent_educ_derived_q1 + parent_cohab_18m + 
                    parent_cohab_3yr + p_age_at_birth + m_age_at_birth + 
                    financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 + 
                    epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr + 
                    bmi_derived_c_14s, 
                  data = fulldata) 

model1a_8yr_unadj <- lm(smfq_dep_c_14s ~ aam_c_14s + smfq_dep_c_8yr, 
                        data = fulldata) 

model1a_8yr_adj <- lm(smfq_dep_c_14s ~ aam_c_14s + smfq_dep_c_8yr + 
                        sex + age_at_q_ret_8yr + parent_income_derived_q1 + 
                        parent_educ_derived_q1 + parent_cohab_18m + 
                        parent_cohab_3yr + p_age_at_birth + m_age_at_birth + 
                        financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 + 
                        epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr + 
                        bmi_derived_c_14s, 
                      data = fulldata) 

models <- list(model1a_unadj, model1a_adj,
               model1a_8yr_unadj,model1a_8yr_adj)


all_res <- map(models, function(res){
  
  #Extract the results of the NHST
  
  nhst_res <- res %>% 
    summary() %>% 
    .$coefficients %>%
    as.data.frame()%>% 
    rownames_to_column("predictor")
  
  #Run the equivalence test (using d=0.23 as SESOI)
  
  equiv_res <- parameters::equivalence_test(
    res, 
    range= c(-effectsize::d_to_r(0.23),effectsize::d_to_r(0.23)), rule="classic")  %>% 
    as.data.frame() %>% 
    `row.names<-`(NULL)   
  
  #Combine 
  
  comb_res <- nhst_res %>% 
    bind_cols(equiv_res) %>% 
    mutate(p_onetailed = ifelse(Estimate<0,`Pr(>|t|)`/2,1)) #One-tailed test p-value 
  
})

#'
#'## Hypothesis 1b
#'
# ----

model1b_unadj <- glm(dep_dx_10_17 ~ aam_c_14s , 
                     family = binomial(link="logit"), 
                     data = fulldata) 

model1b_adj <- glm(dep_dx_10_17 ~ aam_c_14s + 
                     sex + age_at_q_ret_8yr + parent_income_derived_q1 + 
                     parent_educ_derived_q1 + parent_cohab_18m + 
                     parent_cohab_3yr + p_age_at_birth + m_age_at_birth + 
                     financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 + 
                     epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr + 
                     bmi_derived_c_14s,
                   family = binomial(link="logit"), 
                   data = fulldata) 

model1b_8yr_unadj <- glm(dep_dx_10_17 ~ aam_c_14s + dep_dx_0_8, 
                         family = binomial(link="logit"),
                         data = fulldata) 

model1b_8yr_adj <- glm(dep_dx_10_17 ~ aam_c_14s + dep_dx_0_8 + 
                         sex + age_at_q_ret_8yr + parent_income_derived_q1 + 
                         parent_educ_derived_q1 + parent_cohab_18m + 
                         parent_cohab_3yr + p_age_at_birth + m_age_at_birth + 
                         financ_probs_18m + scl_5item_m_q1 + scl_full_m_q3 + 
                         epds_full_m_6m + n_children_8yr + bmi_derived_c_8yr + 
                         bmi_derived_c_14s, 
                       family = binomial(link="logit"),
                       data = fulldata) 

models <- list(model1b_unadj, model1b_adj,
               model1b_8yr_unadj,model1b_8yr_adj)


all_res <- map(models, function(res){
  
  #Extract the results of the NHST
  
  nhst_res <- res %>% 
    summary() %>% 
    .$coefficients %>%
    as.data.frame()%>% 
    rownames_to_column("predictor")
  
  #Run the equivalence test (using d=0.23 as SESOI)
  
  equiv_res <- parameters::equivalence_test(
    res, range= c(log(d_to_oddsratio(-0.23)),log(d_to_oddsratio(0.23))), rule="classic")  %>% 
    as.data.frame() %>% 
    `row.names<-`(NULL)   
  
  #Combine 
  
  comb_res <- nhst_res %>% 
    bind_cols(equiv_res) %>% 
    mutate(p_onetailed = ifelse(Estimate<0,`Pr(>|t|)`/2,1)) #One-tailed test p-value 
  
})
#'
# ---- eval=F,include=F
#This is the code to render and knit the script
ezknitr::ezspin("scripts/02_aim1_analyses.R", 
                out_dir="scripts/reports",
                keep_rmd = T, 
                verbose=T)

