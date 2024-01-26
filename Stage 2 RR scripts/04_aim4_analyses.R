#04_aim4_analyses.R

#Aim 4: Determine whether any associations between earlier age at menarche 
#and symptoms/diagnoses in other domains are likely to be causal

#The purpose of this script is to run primary analyses corresponding to Aim 4. 
#Specifically, we:

# - Run a 2SLS one-sample MR for causal effect of aam on sx 
# - Run a 2SLS one-sample MR for causal effect of aam on childhood sx (negative control)
# - Run a 2SLS (logistic 2nd stage) one-sample MR for causal effect of aam on on all dx variables 

#Read in the data and load packages:
load(file="data/datlist_tr.RData")
library(tidyverse)
library(ggplot2)
library(purrr)
library(forcats)
library(ivreg)
library(effectsize)
library(AER)
library(modelr)
library(broom)
library(sandwich)
library(patchwork)

##Hypotheses 4.1-4a

#This function takes a list of multiple imputed datasets and applies
#instrumental variable regression to each dataset using a specified formula.
#It then calculates robust standard errors for the regression estimates
#and pools the results across all imputed datasets. 
adjust_and_pool <- function(imputed_data_list, formula) {
  
  #placeholder lists for storing estimates and standard errors across datasets
  estimates_list <- vector("list", length(imputed_data_list))
  se_list <- vector("list", length(imputed_data_list))
  
  #loop over each imputed dataset
  for(i in 1:length(imputed_data_list)) {
    
    #run the instrumental variable regression (ivreg) using the specified formula
    model <- ivreg(formula, data = imputed_data_list[[i]])
    
    #calculate robust standard errors using vcovHC function from sandwich package
    robust_se <- sqrt(diag(vcovHC(model, type = "HC1")))
    
    #store the estimates and their robust standard errors for each dataset
    estimates_list[[i]] <- coef(model)
    se_list[[i]] <- robust_se
  }
  
  #pooling the results from each imputed dataset
  
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
  two_tailed_p <- 2*(pt(-abs(z_val), df))
  
  #combine the results into a data frame
  comb_res <- data.frame(
    Estimate = pooled_coef,
    se = pooled_se,
    z = z_val,
    p = two_tailed_p,
    CI_low = pooled_coef - qnorm(0.975) * pooled_se, # Lower bound of 95% CI
    CI_high = pooled_coef + qnorm(0.975) * pooled_se  # Upper bound of 95% CI
  )
  
  return(comb_res)
  
}

#list the symptom variables
domains <-list(anx=c("scared_anx_c_8yr","scared_anx_c_14c"),
               adhd = c("rsdbd_adhd_c_8yr","rsdbd_adhd_c_14m"),
               cd=c("rsdbd_cd_c_8yr","rsdbd_cd_c_14c"),
               odd=c("rsdbd_odd_c_8yr","rsdbd_odd_c_14m"))

#run models for all domains by iterating through the 50 imputed datasets
domains_res <- list()
n = 0

for(domain in domains) {
  n = n + 1
  
  #specify models for 8-year outcomes
  results_8yr <- adjust_and_pool(datlist,
                                 formula = as.formula(paste(domain[[1]], "~ aam_c_14c | PGS_AaM")))
  
  #specify models for 14-year outcomes
  results_14yr <- adjust_and_pool(datlist,
                                  formula = as.formula(paste(domain[[2]], "~ aam_c_14c | PGS_AaM")))
  
  results <- list(results_8yr, results_14yr)
  
  domains_res[[names(domains)[[n]]]] <- results
  
}

ests_8y <- map(domains_res, `[[`, 1)
ests_14y <- map(domains_res, `[[`, 2)

#reformat 8-year results
anx_sx_8y <- ests_8y$anx %>%
  mutate(domain = "ANX")
cd_sx_8y <- ests_8y$cd %>%
  mutate(domain = "CD")
odd_sx_8y <- ests_8y$odd %>%
  mutate(domain = "ODD")
adhd_sx_8y <- ests_8y$adhd %>%
  mutate(domain = "ADHD")

ests_4a_8y <- anx_sx_8y %>%
  full_join(adhd_sx_8y) %>% 
  full_join(cd_sx_8y) %>% 
  full_join(odd_sx_8y) %>%
  mutate(outcome = "Symptom outcomes",
         adj = "1-sample MR")

#reformat 14-year results
anx_sx_14y <- ests_14y$anx %>%
  mutate(domain = "ANX")
cd_sx_14y <- ests_14y$cd %>%
  mutate(domain = "CD")
odd_sx_14y <- ests_14y$odd %>%
  mutate(domain = "ODD")
adhd_sx_14y <- ests_14y$adhd %>%
  mutate(domain = "ADHD")

ests_4a_14y <- anx_sx_14y %>%
  full_join(adhd_sx_14y) %>% 
  full_join(cd_sx_14y) %>% 
  full_join(odd_sx_14y) %>%
  mutate(ROPE_low =- d_to_r(0.20),
         ROPE_high = d_to_r(0.20),
         outcome = "Symptom outcomes",
         adj = "1-sample MR")

##Perform equivalence testing across domains using a function

#prepare a dataframe to store the results
results <- data.frame(
  Domain = character(), Equivalence = logical(), CI_upper_14y = numeric(),
  Equivalence_Bound = numeric(), stringsAsFactors = FALSE
  )

perform_equivalence_test <- function(ests_4a_14y, ests_4a_8y) {
  
  #loop through the domains to perform the equivalence tests
  for (domain in unique(ests_4a_14y$domain)) {
    #extract the estimate and standard error for the current domain
    estimate_14y <- ests_4a_14y$Estimate[ests_4a_14y$domain == domain]
    std_error_14y <- ests_4a_14y$se[ests_4a_14y$domain == domain]
    
    #calculate the upper confidence interval for the 14-year estimate
    t_critical <- qt(1 - 0.1 / 2, df = Inf) #using 90% CI
    
    ci_lower_14y <- estimate_14y - t_critical * std_error_14y
    
    # Extract the equivalence bound for the current domain from the 8-year CI
    equivalence_bound <- ests_4a_8y$CI_low[ests_4a_8y$domain == domain]
    
    # Perform the equivalence test
    is_equivalent <- ci_lower_14y > equivalence_bound
    
    # Append the results to the dataframe
    results <- rbind(results, data.frame(Domain = domain, Equivalence = is_equivalent, 
                                         CI_lower_14y = ci_lower_14y, 
                                         Equivalence_Bound = equivalence_bound))
  }
  
  return(results)
}

#run the function
equivalence_results_4a <- perform_equivalence_test(ests_4a_14y, ests_4a_8y)

#print results
print(equivalence_results_4a)

##Hypotheses 4.1-3b

#Function to perform two-stage regression and adjust SEs for adolescent dx
two_stage_with_adjSEs_ad <- function(data, cluster_id) {
  
  #first stage: linear regression
  lm_model <- lm(aam_c_14c ~ PGS_AaM, data = data)
  predictions <- predict(lm_model, newdata = data)
  data$pred_aam <- predictions

  #second stage: logistic regression with the predictions from the first stage
  logit_model <- glm(as.formula(paste(domain[[2]], "~ pred_aam")), family = binomial(link = "logit"), data = data)
  
  #get robust standard errors accounting for uncertainty in first stage
  robust_se <- sqrt(diag(vcovHC(logit_model, type = "HC1"))['pred_aam'])
  
  #extract coefficient for the predictor from the second stage
  estimate <- coef(logit_model)['pred_aam']
  
  #calculate z-value, two-tailed p-value, and confidence interval for the IV estimate
  z <- estimate / robust_se
  pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
  ci_low <- estimate - qnorm(0.975) * robust_se
  ci_high <- estimate + qnorm(0.975) * robust_se
  
  #combine IV results into a data frame
  results <- data.frame(
    Estimate = estimate,
    Robust_SE = robust_se,
    Z = z,
    P_Value = pval,
    CI_Low = ci_low,
    CI_High = ci_high
  )
  
  return(results)
}

#list the diagnosis variables
domains <-list(anx=c("anx_dx_ch","anx_dx_ad"),
               adhd = c("adhd_dx_ch","adhd_dx_ad"),
               beh=c("beh_dx_ch","beh_dx_ad"))

#run adolescent dx models for all domains across the 50 imputed datasets, 
#adjusting the standard errors based on clustering within mothers
domains_res_ad <- list()
n = 0

for(domain in domains) {
  n = n + 1

#apply this function to each imputed dataset
  results_list <- lapply(datlist, function(data) two_stage_with_adjSEs_ad(data))

  m <- length(datlist)
  
#pool the results from each imputed dataset using Rubin's rules 
  pooled_estimate <- mean(sapply(results_list, function(x) x$Estimate))
  within_var <- mean(sapply(results_list, function(x) x$Robust_SE)^2)*m/(m - 1)
  between_var <- var(sapply(results_list, function(x) x$Estimate))*(m - 1)/m
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

##Note: models with childhood diagnoses were deleted due to non-convergence -
##due to low case numbers, these cannot be used as negative control outcomes

#reformat adolescent dx results
anx_dx_ad <- domains_res_ad$anx %>%
  mutate(domain = "ANX")
dbd_dx_ad <- domains_res_ad$beh %>%
  mutate(domain = "DBD")
adhd_dx_ad <- domains_res_ad$adhd %>%
  mutate(domain = "ADHD")

#exponentiate the dx outcome results
ests_4b_ad <- anx_dx_ad %>%
  full_join(dbd_dx_ad) %>%
  full_join(adhd_dx_ad) %>%
  mutate(Estimate = exp(Estimate),
         CI_low = exp(CI_low),
         CI_high = exp(CI_high),
         outcome = "Diagnoses age 10-17")

#merge MR results for symptoms and diagnoses
ests4_ad <- ests_4a_14y %>%
  full_join(ests_4b_ad) %>%
  mutate(adj="1-sample MR")

#add negative control MR results
ests4 <- ests_4a_8y %>%
  mutate(adj="Negative control MR") %>%
  full_join(ests4_ad)

#load observational results for other domains
load(file="./data/aim2_results.RData")

ests2_ad <- ests2 %>%
  mutate(adj = case_when(outcome=="14-year symptoms" &
                           cov_adj=="Covariate\nunadjusted" &
                           child_adj=="Childhood\nunadjusted" &
                           dep_adj=="Depression\nunadjusted" ~ "Unadjusted",
                         outcome=="14-year symptoms" &
                           cov_adj=="Covariate\nadjusted" &
                           child_adj=="Childhood\nadjusted" & 
                           dep_adj=="Depression\nadjusted" ~ "Adjusted",
                         outcome=="Diagnoses age 10-17" &
                           cov_adj=="Covariate\nunadjusted" & 
                           child_adj=="Childhood\nunadjusted" &
                           dep_adj=="Depression\nunadjusted" ~ "Unadjusted",
                         outcome=="Diagnoses age 10-17" &
                           cov_adj=="Covariate\nadjusted" &
                           child_adj=="Childhood\nadjusted" &
                           dep_adj=="Depression\nadjusted"~ "Adjusted")) %>%
  select(-c(cov_adj,child_adj,dep_adj)) %>%
  mutate(outcome = ifelse(outcome == "14-year symptoms", "Symptom outcomes", outcome)) %>%
  na.omit()

#merge observational with MR results
ests_all_ad <- ests4 %>%
  full_join(ests2_ad) %>%
  select(-model,-z) %>%
  mutate(outcome = ifelse(outcome == "Diagnoses age 10-17", "Diagnostic outcomes", outcome))

#round to 2 decimal places for reporting in main text
ests_all_2dp <- ests_all_ad
ests_all_2dp[] <- lapply(ests_all_2dp, function(x) if(is.numeric(x)) round(x, 2) else x)
ests_all_2dp

#visualise symptom outcome results
p4a<- ests_all_ad %>%
  filter(outcome=="Symptom outcomes") %>%
  mutate(adj=fct_relevel(adj,"Negative control MR","1-sample MR","Adjusted","Unadjusted")) %>%
  mutate(domain = fct_relevel(domain,"ADHD","ODD","CD","ANX")) %>%
  ggplot(mapping=aes(x=Estimate,y=domain,xmin=CI_low,xmax=CI_high,colour=adj,shape=adj)) +
  geom_vline(aes(xintercept=0), colour = "grey80", linetype = 2) + theme_light(base_size=18) +
  geom_vline(aes(xintercept=-0.09950372), colour = "#0072B2", linetype = 2) + 
  geom_vline(aes(xintercept=0.09950372), colour = "#0072B2", linetype = 2) + 
  geom_vline(aes(xintercept=-0.10934048), colour = "#D55E00", linetype = 2) + 
  geom_vline(aes(xintercept=0.10934048), colour = "#D55E00", linetype = 2) +   
  geom_errorbar(size = 1.4, alpha=0.6, width = 0, position=position_dodge(0.85)) +
  geom_point(size=4.1, fill="white", alpha=1, position=position_dodge(0.85)) +
  facet_grid(rows = vars(outcome), scales = "fixed", space = "fixed") +
  scale_shape_manual(values=c(18,17,15,16)) +
  scale_colour_manual(values=c("#C9C9C9","#0072B2","#D55E00","#942A12")) +
  guides(colour=guide_legend(order=0, reverse=TRUE),shape=guide_legend(order=0, reverse=TRUE)) + 
  geom_text(label="Region of practical equivalence to zero", x=0, y=4.9, col="black", size=5.3) +
  theme(axis.text = element_text(size=15, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(0,0,-2,0,"cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        strip.text = element_text(size=15, colour = 'black'),
        strip.background = element_rect(colour="gray", fill="white")) + 
  xlab(label=expression("Standardised"~beta)) + 
  coord_cartesian(xlim = c(-0.2,0.2), clip = 'off') +
  scale_x_continuous(breaks = c(-0.2, -0.1, 0.0, 0.1, 0.2))

#visualise diagnostic outcome results
p4b <- ests_all_ad %>%
  filter(outcome=="Diagnostic outcomes") %>%
  mutate(domain = fct_relevel(domain,"ADHD","DBD","ANX")) %>%
  ggplot(mapping=aes(x=Estimate,y=domain,xmin=CI_low,xmax=CI_high,colour=adj,shape=adj)) +
  geom_vline(aes(xintercept=1), colour = "grey80", linetype = 2) + theme_light(base_size=18) +
  geom_vline(aes(xintercept=0.69575348), colour = "#0072B2", linetype = 2) + 
  geom_vline(aes(xintercept=1.43729069), colour = "#0072B2", linetype = 2) + 
  geom_vline(aes(xintercept=0.67096664), colour = "#D55E00", linetype = 2) + 
  geom_vline(aes(xintercept=1.49038706), colour = "#D55E00", linetype = 2) + 
  geom_errorbar(size = 1.4, alpha=0.6, width = 0, position=position_dodge(0.8)) +
  geom_point(size=4.1, fill="white", alpha=1, position=position_dodge(0.8)) +
  facet_grid(rows = vars(outcome), scales = "free", space = "fixed") +
  scale_colour_manual(values=c("#0072B2","#D55E00","#942A12")) +
  scale_shape_manual(values=c(17,15,16)) +
  guides(shape="none",colour="none") +
  theme(axis.text = element_text(size=15, colour = "black"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill="transparent"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.margin = margin(0.75, 0, 0, 0, "cm"),
        strip.text = element_text(size=15, colour = 'black'),
        strip.background = element_rect(colour="gray", fill="white")) + 
  xlab(label=expression("Odds ratio (log scale)")) + 
  coord_cartesian(xlim = c(0.48,2.07), clip = 'off') +
  scale_x_log10(breaks = c(0.5,0.7,1,1.4,2))

patch3 <- p4a / p4b

patch4 <- patch3 + plot_annotation(tag_levels = "A") &
  plot_layout(guides = "collect", nrow = 2, heights = c(4.5,3.2)) &
  theme(legend.position = "top",
        legend.box.margin = margin(t = 0, r = 0, b = 20, l = 0))

patch4

tiff("figures/all_domains_tr.tiff", res=800, compression="lzw", unit="in",
     height=8, width=8)

patch4

dev.off()
