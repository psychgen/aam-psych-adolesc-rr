#07_IPW_sensitivity.R

#This script runs inverse probability weighting sensitivity analyses,
#sourcing the script `07.1_IPW_weights.R` where weights are generated.

#read in the data and load packages:
load(file="data/datlist_tr.RData")
library(tidyverse)
library(effectsize)
library(miceadds)
library(patchwork)
library(AER)

#source script where IP weights are generated (note that this 
#includes multiple imputation which can run for some time)
source("./scripts/07.1_IPW_weights.R")

#load processed data with IP weights
load(file="data/ipw_wdat.RData")

#add weights to the imputed datasets
datlist_ipw <- lapply(datlist, function(dataset) {
  left_join(dataset, select(part_wdat, ind_id, ipw), by = "ind_id")
})

#H1a

#note that there is a bug in lm.cluster when running it inside a larger 
#function ("object 'wgt__' not found"), therefore we use iteration here

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

for(model in names(models)) {
  #initialize a sub-list for the current model
  lm_results[[model]] <- list()
  
  #loop through each dataset in the list of imputed datasets
  for(i in seq_along(datlist_ipw)) {
    lm_results[[model]][[i]] <- miceadds::lm.cluster(
      data = datlist_ipw[[i]],
      formula = models[[model]],
      cluster = datlist_ipw[[i]]$m_id,
      weights = datlist_ipw[[i]]$ipw
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
         CI_low = lci,
         CI_high = uci,
         ROPE_low = - d_to_r(0.23),
         outcome = "14-year symptoms") %>%
  select(Estimate,se,z,p,p_onetailed,CI_low,CI_high,
         ROPE_low,outcome,cov_adj,child_adj)

#H1b

#initialize a list to store the models
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
  for(i in seq_along(datlist_ipw)) {
    glm_results[[model]][[i]] <- miceadds::glm.cluster(
      data = datlist_ipw[[i]],
      formula = models[[model]],
      family=binomial(link="logit"),
      cluster = datlist_ipw[[i]]$m_id,
      weights = datlist_ipw[[i]]$ipw
    )
  }
}

#extract betas and vars
betas <- lapply(glm_results, function(models) {
  lapply(models, function(model) coef(model))
})

variances <- lapply(glm_results, function(models) {
  lapply(models, function(model) vcov(model))
})

#pool the results
model1b_unadj <- summary(pool_mi(qhat=betas$model1b_unadj, u=variances$model1b_unadj))
model1b_adj <- summary(pool_mi(qhat=betas$model1b_adj, u=variances$model1b_adj))

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
         CI_low = exp(lci),
         CI_high = exp(uci),
         ROPE_low = d_to_oddsratio(-0.23),
         outcome = "Diagnoses age 10-17",
         child_adj = "Childhood\nunadjusted") %>%
  select(Estimate,se,z,p,p_onetailed,CI_low,CI_high,
         ROPE_low,outcome,cov_adj,child_adj)

#merge results for symptoms and diagnoses
ests1_ipw <- ests1a %>%
  full_join(ests1b)

#save out
save(ests1_ipw, file="./data/aim1_results_ipw.RData")
load(file="./data/aim1_results_ipw.RData")

#H3a

#This function takes a list of multiply imputed datasets and applies weighted
#instrumental variable regression to each dataset using a specified formula.
#It then calculates robust standard errors for the regression 
#estimates and pools the results across all imputed datasets.
adjust_and_pool <- function(imputed_data_list, formula) {
  
  #placeholder lists for storing estimates and standard errors across imputed datasets
  estimates_list <- vector("list", length(imputed_data_list))
  se_list <- vector("list", length(imputed_data_list))
  
  #loop over each imputed dataset
  for(i in 1:length(imputed_data_list)) {
    
    #run the instrumental variable regression (ivreg) using the specified formula
    model <- ivreg(formula, data = imputed_data_list[[i]], weights = datlist_ipw[[i]]$ipw)
    
    #get robust standard errors accounting for uncertainty in first stage
    robust_se <- sqrt(diag(vcovHC(model, type = "HC1")))
    
    #store the estimates and their robust standard errors for each dataset
    estimates_list[[i]] <- coef(model)
    se_list[[i]] <- robust_se
  }
  
  # Pooling the results from each imputed dataset
  
  #number of imputations
  m <- length(imputed_data_list)
  
  #calculate the pooled estimate for the coefficient of interest (second coefficient in the model)
  pooled_coef <- mean(sapply(estimates_list, function(x) x[2]))
  
  #calculate the within-imputation variance (average of squared standard errors)
  within_variance <- mean(sapply(se_list, function(x) x[2]^2)) * m/(m - 1)
  
  #calculate the between-imputation variance (variance of estimates)
  between_variance <- var(sapply(estimates_list, function(x) x[2])) * (m - 1)/m
  
  #compute the total variance of the pooled estimate
  total_variance <- within_variance + (1 + 1/m) * between_variance
  pooled_se <- sqrt(total_variance)
  
  #calculate degrees of freedom using Barnard and Rubin's adjustment
  df <- (m - 1) * (1 + (within_variance / total_variance))^2
  
  #compute test statistics for the pooled estimate
  z_val <- pooled_coef / pooled_se
  
  #calculate one-tailed and two-tailed p-values
  one_tailed_p <- pt(-abs(z_val), df) 
  two_tailed_p <- 2 * one_tailed_p
  
  #combine the results into a data frame
  comb_res <- data.frame(
    Estimate = pooled_coef,
    se = pooled_se,
    z = z_val,
    two_tailed_p = two_tailed_p,
    one_tailed_p = one_tailed_p,
    CI_low = pooled_coef - qnorm(0.975) * pooled_se, # Lower bound of 95% CI
    CI_high = pooled_coef + qnorm(0.975) * pooled_se  # Upper bound of 95% CI
  )
  
  return(comb_res)
  
}

#run models on all imputed datasets, with robust standard errors, and 
#accounting for participation weights
ests_3a_14_ipw <- adjust_and_pool(datlist_ipw, 
                                  formula = smfq_dep_c_14c ~ aam_c_14c | PGS_AaM) %>%
  mutate(ROPE_low =- d_to_r(0.25),
         outcome="14-year symptoms")

#H3b

# Function to perform two-stage regression and adjust SEs
two_stage_with_adjSEs <- function(data) {
  
  #first stage: linear regression
  lm_model <- lm(aam_c_14c ~ PGS_AaM, data=data, weights = data$ipw)
  data$pred_aam <- predict(lm_model, newdata = data)
  
  #second stage: logistic regression with the predictions from the first stage
  logit_model <- glm(dep_dx_ad ~ pred_aam, family=binomial(link="logit"), data=data, weights = data$ipw)
  
  #extract coefficient for the predictor from the second stage
  estimate_iv <- coef(logit_model)['pred_aam']
  
  #get robust standard errors accounting for uncertainty in first stage
  robust_se <- sqrt(diag(vcovHC(logit_model, type = "HC1"))['pred_aam'])
  
  #calculate z-value, two-tailed p-value, and confidence interval for the estimate
  z_iv <- estimate_iv / robust_se
  pval_iv <- 2 * pnorm(abs(z_iv), lower.tail = FALSE)
  ci_low_iv <- estimate_iv - qnorm(0.975) * robust_se
  ci_high_iv <- estimate_iv + qnorm(0.975) * robust_se
  
  #combine IV results into a data frame
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

#apply this function to each imputed dataset
results_list <- lapply(datlist_ipw, function(data) two_stage_with_adjSEs(data))

#number of imputations
m <- length(datlist_ipw)

#pool the results from each imputed dataset using Rubin's rules 
pooled_estimate <- mean(sapply(results_list, function(x) x$Estimate_IV))
within_var <- mean(sapply(results_list, function(x) x$Robust_SE_IV)^2) * m/(m - 1)
between_var <- var(sapply(results_list, function(x) x$Estimate_IV)) * (m - 1)/m
total_var <- within_var + (1 + 1/m) * between_var
pooled_se <- sqrt(total_var)

#calculate degrees of freedom using Barnard and Rubin's adjustment
df <- (m - 1) * (1 + (within_var / total_var))^2

#calculate adjusted one-tailed and two-tailed p-values
one_tailed_p <- pt(-abs(pooled_estimate / pooled_se), df) 
two_tailed_p <- 2 * one_tailed_p

#constructing pooled results
pooled_results_3b_ipw <- data.frame(
  Estimate = pooled_estimate,
  se = pooled_se,
  z = pooled_estimate / pooled_se,
  two_tailed_p = two_tailed_p,
  one_tailed_p = one_tailed_p,
  CI_low = pooled_estimate - qnorm(0.975) * pooled_se,
  CI_high = pooled_estimate + qnorm(0.975) * pooled_se,
  ROPE_low = d_to_oddsratio(-0.25)
)

#compute additional statistics, exponentiate,
#and rename to enable merging with MR results
ests_3b_14_ipw <- pooled_results_3b_ipw %>%
  mutate(Estimate = exp(Estimate),
         CI_low = exp(CI_low),
         CI_high = exp(CI_high),
         outcome = "Diagnoses age 10-17")

#merge MR results for symptoms and diagnoses
ests3_ipw <- ests_3a_14_ipw %>%
  full_join(ests_3b_14_ipw) %>%
  mutate(adj="1-sample MR")

#load observational results for depression
load(file="./data/aim1_results_ipw.RData")

ests1_ipw <- ests1_ipw %>%
  mutate(two_tailed_p = p,
         one_tailed_p = p_onetailed,
         adj = case_when(outcome=="14-year symptoms" &
                           cov_adj=="Covariate\nunadjusted" &
                           child_adj=="Childhood\nunadjusted" ~ "Unadjusted",
                         outcome=="14-year symptoms" &
                           cov_adj=="Covariate\nadjusted" &
                           child_adj=="Childhood\nadjusted" ~ "Adjusted",
                         outcome=="Diagnoses age 10-17" &
                           cov_adj=="Covariate\nunadjusted" ~ "Unadjusted",
                         outcome=="Diagnoses age 10-17" &
                           cov_adj=="Covariate\nadjusted" ~ "Adjusted")) %>%
  select(-c(cov_adj,child_adj,p,p_onetailed)) %>%
  na.omit()

#merge observational and 
#MR results for plotting
ests_dep_ipw <- ests3_ipw %>%
  full_join(ests1_ipw)

#restrict to 2 decimal points for reporting in text
ests_dep_ipw_2dp <- ests_dep_ipw
ests_dep_ipw_2dp[] <- lapply(ests_dep_ipw_2dp, function(x) if(is.numeric(x)) round(x, 2) else x)
ests_dep_ipw_2dp

#plot results
p3a_ipw <- ests_dep_ipw %>%
  filter(outcome == "14-year symptoms") %>%
  ggplot(mapping=aes(x=Estimate,y="",xmin=CI_low,xmax=CI_high,colour=adj,shape=adj)) +
  geom_vline(aes(xintercept=0), colour = "grey80", linetype = 2) + theme_light(base_size=14) +
  geom_vline(aes(xintercept=ROPE_low), colour = c("#0072B2","#D55E00","#D55E00"), linetype = 2, size=0.5) + 
  geom_errorbar(size = 1.4, alpha=0.6, width = 0, position=position_dodge(0.4)) +
  geom_point(size=4.3, fill="white", alpha=1, position=position_dodge(0.4)) +
  geom_text(aes(label = "Smallest effect size of interest"), x=-0.12, y=1.740, col="black", size=4.9) +
  scale_colour_manual(values=c("#0072B2","#D55E00","#942A12")) +
  scale_shape_manual(values=c(17,15,16)) +
  guides(colour=guide_legend(order=0, reverse=TRUE),shape=guide_legend(order=0, reverse=TRUE)) + 
  theme(axis.text = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14, colour = "black", angle = 0, vjust = 0.5),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        axis.line = element_line(colour="grey80"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin = margin(1.5,1,-1,1,"cm"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.background = element_rect(colour="gray", fill="white")) + 
  xlab(label=expression("Standardised"~beta)) +
  ylab(label="Depressive\nsymptoms") +
  coord_cartesian(xlim = c(-0.2,0.1), clip='off') + 
  scale_x_continuous(breaks = c(-0.2, -0.1, 0.0, 0.1))

p3b_ipw<- ests_dep_ipw %>%
  filter(outcome == "Diagnoses age 10-17") %>%
  ggplot(mapping=aes(x=Estimate,y="",xmin=CI_low,xmax=CI_high,colour=adj,shape=adj)) +
  geom_vline(aes(xintercept=1), colour = "grey80", linetype = 2) + theme_light(base_size=14) +
  geom_vline(aes(xintercept=ROPE_low), colour = c("#0072B2","#D55E00","#D55E00"), linetype = 2, size=0.5) + 
  geom_errorbar(size = 1.4, alpha=0.6, width = 0, position=position_dodge(0.4)) +
  geom_point(size=4.3, fill="white", alpha=1, position=position_dodge(0.4)) +
  scale_colour_manual(values=c("#0072B2","#D55E00","#942A12")) +
  scale_shape_manual(values=c(17,15,16)) +
  guides(shape="none",colour="none") +
  theme(axis.text = element_text(size=14, colour = "black"),
        axis.title.y = element_text(size=14, colour = "black", angle = 0, vjust = 0.5),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        plot.margin = margin(-1,1,1,1,"cm"),
        axis.line = element_line(colour="grey80"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="gray", fill="white")) + 
  xlab(label=expression("Odds ratio (log scale)")) + 
  ylab(label="Depression\ndiagnoses") +
  coord_cartesian(xlim = c(0.482,1.43), clip = 'off') +
  scale_x_log10(breaks = c(0.5,0.7,1,1.4))

patch1 <- p3a_ipw / p3b_ipw

patch_ipw <- patch1 + plot_annotation(tag_levels = "A") &
  plot_layout(guides = "collect") &
  theme(legend.position = "top",
        legend.box.margin = margin(t = 0, r = 0, b = 15, l = 0),
        legend.title = element_blank())

patch_ipw

tiff("figures/dep_MR_ipw_tr.tiff", res=800, compression="lzw", unit="in",
     height=6, width=8)

patch_ipw

dev.off()
