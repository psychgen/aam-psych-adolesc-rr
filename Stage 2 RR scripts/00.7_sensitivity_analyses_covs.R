#00.7_sensitivity_analyses_covs.R

##This script conducts sensitivity analyses regressing the genetic
##instrument for age at menarche on all covariates.

#read in the data and load packages:
load(file="./data/imputed_datasets.RData")
library(tidyverse)
library(mice)
library(miceadds)

covariates <-list("parent_education_q1","parent_income_q1","parent_cohab_18m",
                  "parent_cohab_3yr","p_age_at_birth","m_age_at_birth",
                  "financ_probs_18m","scl_short_m_q1","scl_full_m_q3",
                  "epds_short_m_6m","parity","bmi_derived_c_8yr",
                  "bmi_derived_c_14c","age_8yr","age_14yr")

#standardise the variables for each imputed dataset
standardised_covdat <- lapply(split(imputed_datasets, imputed_datasets$.imp), function(df) {
  df %>% mutate(parent_cohab_18m = as.numeric(parent_cohab_18m),
                parent_cohab_3yr = as.numeric(parent_cohab_3yr)) %>%
    mutate(across(c("PGS_AaM", unlist(covariates)), ~as.vector(scale(.))))
})

#convert the list of data frames back to mids object
standardised_cov <- mice::as.mids(do.call(rbind, standardised_covdat))

#check prediction of covariates in linear/logistic regressions
model1 <- with(standardised_cov, lm(parent_education_q1~PGS_AaM))
res1 <- summary(pool(model1))

model2 <- with(standardised_cov, lm(parent_income_q1~PGS_AaM))
res2 <- summary(pool(model2))

model3 <- with(standardised_cov, lm(parent_cohab_18m~PGS_AaM))
res3 <- summary(pool(model3))

model4 <- with(standardised_cov, lm(parent_cohab_3yr~PGS_AaM))
res4 <- summary(pool(model4))

model5 <- with(standardised_cov, lm(p_age_at_birth~PGS_AaM))
res5 <- summary(pool(model5))

model6 <- with(standardised_cov, lm(m_age_at_birth~PGS_AaM))
res6 <- summary(pool(model6))

model7 <- with(standardised_cov, lm(financ_probs_18m~PGS_AaM))
res7 <- summary(pool(model7))

model8 <- with(standardised_cov, lm(scl_short_m_q1~PGS_AaM))
res8 <- summary(pool(model8))

model9 <- with(standardised_cov, lm(scl_full_m_q3~PGS_AaM))
res9 <- summary(pool(model9))

model10 <- with(standardised_cov, lm(epds_short_m_6m~PGS_AaM))
res10 <- summary(pool(model10))

model11 <- with(standardised_cov, lm(parity~PGS_AaM))
res11 <- summary(pool(model11))

model12 <- with(standardised_cov, lm(bmi_derived_c_8yr~PGS_AaM))
res12 <- summary(pool(model12))

model13 <- with(standardised_cov, lm(bmi_derived_c_14c~PGS_AaM))
res13 <- summary(pool(model13))

model14 <- with(standardised_cov, lm(age_8yr~PGS_AaM))
res14 <- summary(pool(model14))

model15 <- with(standardised_cov, lm(age_14yr~PGS_AaM))
res15 <- summary(pool(model15))

# Function to process a single model summary
process_pooled_models <- function(model_summary) {
  # Extracting the necessary information
  estimates <- model_summary$estimate
  std_errors <- model_summary$std.error
  lower_bounds <- estimates - 1.96 * std_errors
  upper_bounds <- estimates + 1.96 * std_errors
  
  # Creating a data frame with results
  results_df <- data.frame(
    Estimate = estimates,
    Std.Error = std_errors,
    LCI = lower_bounds,
    UCI = upper_bounds,
    P = model_summary$p.value
  )
  
  # Format to 2 decimal places
  results_df <- round(results_df, 2)
  
  return(results_df)
}

# Process each pooled model summary and store the results
all_results <- list(
  res1 = process_pooled_models(res1),
  res2 = process_pooled_models(res2),
  res3 = process_pooled_models(res3),
  res4 = process_pooled_models(res4),
  res5 = process_pooled_models(res5),
  res6 = process_pooled_models(res6),
  res7 = process_pooled_models(res7),
  res8 = process_pooled_models(res8),
  res9 = process_pooled_models(res9),
  res10 = process_pooled_models(res10),
  res11 = process_pooled_models(res11),
  res12 = process_pooled_models(res12),
  res13 = process_pooled_models(res13),
  res14 = process_pooled_models(res14),
  res15 = process_pooled_models(res15)
)

#save out
save(all_results, file="output/sensitivity_cov_results.RData")
