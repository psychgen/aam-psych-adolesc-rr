#04_aim3_analyses.R

#Aim 3: Determine whether any associations between earlier age at menarche and 
#adolescent depressive symptoms/depression are likely to be causal

#The purpose of this script is to run primary analyses corresponding to Aim 4. 
#Specifically, we:

# - Run a 2SLS one-sample MR for causal effect of aam on dep sx 
# - Run a 2SLS one-sample MR for causal effect of aam on earlier dep sx (negative control)
# - Run a 2SLS (logistic 2nd stage) one-sample MR for causal effect of aam on on dep dx  

#Read in the data and load packages:
load(file="data/datlist_tr.RData")
load(file="./data/standardised_MI_tr.RData")
library(tidyverse)
library(ggplot2)
library(ivreg)
library(parameters)
library(lmtest)
library(effectsize)
library(AER)
library(modelr)
library(broom)
library(sandwich)
library(mice)

##Hypothesis 3a

#This function takes a list of multiply imputed datasets and applies
#instrumental variable regression to each dataset using a specified formula.
#It then calculates robust standard errors for the 2SLS (using option = HC1)
#estimates and pools the results across all imputed datasets.
adjust_and_pool <- function(imputed_data_list, formula) {
  
  #placeholder lists for storing estimates and standard errors across imputed datasets
  estimates_list <- vector("list", length(imputed_data_list))
  se_list <- vector("list", length(imputed_data_list))
  
  #loop over each imputed dataset
  for(i in 1:length(imputed_data_list)) {
    
    #run the instrumental variable regression (ivreg) using the specified formula
    model <- ivreg(formula, data = imputed_data_list[[i]])
    
    #calculate robust standard errors using vcovHC function from the sandwich package
    robust_se <- sqrt(diag(vcovHC(model, type = "HC1")))
    
    #store the coefficient estimates and robust standard errors for each dataset
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
    CI90_low = pooled_coef - qnorm(0.95) * pooled_se, # Lower bound of 90% CI
    CI90_high = pooled_coef + qnorm(0.95) * pooled_se, # Upper bound of 90% CI
    CI95_low = pooled_coef - qnorm(0.975) * pooled_se, # Lower bound of 95% CI
    CI95_high = pooled_coef + qnorm(0.975) * pooled_se # Upper bound of 95% CI
  )
  
  return(comb_res)
  
}

#run models on all imputed datasets, with robust standard errors
ests_3a_14 <- adjust_and_pool(datlist, 
                              formula = smfq_dep_c_14c ~ aam_c_14c | PGS_AaM) %>%
  mutate(ROPE_low =- d_to_r(0.25),
         outcome="14-year symptoms")

ests_3a_8y <- adjust_and_pool(datlist, 
                              formula = smfq_dep_c_8yr ~ aam_c_14c | PGS_AaM) %>%
  mutate(ROPE_low =- d_to_r(0.25),
         outcome="8-year symptoms")

#sanity check: does the results of this custom function match mice results?
#
#  #run ivreg using with()
#  model <- with(standardised_MI, ivreg(smfq_dep_c_14c ~ aam_c_14c | PGS_AaM))
#  
#  #pool the results of the analyses
#  pooled_results <- pool(model)
#  
#  #summary of pooled results
#  pooled <- summary(pooled_results) %>% #se = 0.03538
#    mutate(one_tailed_p = ifelse(estimate < 0, p.value / 2, 1 - p.value / 2))

##Perform equivalence testing using the lower bound of negative control estimate
perform_one_sided_equivalence_test <- function(ests_3a_14, ests_3a_8y) {
  
  #perform the one-sided equivalence test
  is_equivalent <- ests_3a_14$CI90_low > ests_3a_8y$CI95_low
  
  # Construct the result
  result <- data.frame(
    Equivalence = is_equivalent,
    Estimate_14y = ests_3a_14$Estimate,
    CI_lower_14y = ests_3a_14$CI90_low,
    Equivalence_Bound = ests_3a_8y$CI95_low
  )
  
  return(result)
}

#run the function
equivalence_results_3a <- perform_one_sided_equivalence_test(ests_3a_14,ests_3a_8y)

#print results
print(equivalence_results_3a)

##Hypothesis 3b
  
#Function to perform two-stage regression and adjust SEs
two_stage_with_adjSEs <- function(data) {
  
  #first stage: linear regression
  lm_model <- lm(aam_c_14c ~ PGS_AaM, data=data)
  data$pred_aam <- predict(lm_model, newdata = data)
  
  #second stage: logistic regression with the predictions from the first stage
  logit_model <- glm(dep_dx_ad ~ pred_aam, family=binomial(link="logit"), data=data)
  
  #extract coefficient for the predictor from the second stage
  estimate <- coef(logit_model)['pred_aam']
  
  #get robust standard errors accounting for uncertainty in the first stage
  robust_se <- sqrt(diag(vcovHC(logit_model, type = "HC1"))['pred_aam'])
  
  #calculate z-value, p-value, and confidence interval for the estimate
  z <- estimate / robust_se
  pval <- 2 * pnorm(abs(z), lower.tail = FALSE)
  ci90_low <- estimate - qnorm(0.95) * robust_se
  ci90_high <- estimate + qnorm(0.95) * robust_se
  ci95_low <- estimate - qnorm(0.975) * robust_se
  ci95_high <- estimate + qnorm(0.975) * robust_se
  
  #combine IV results into a data frame
  results <- data.frame(
    Estimate = estimate,
    Robust_SE = robust_se,
    Z = z,
    P_Value = pval,
    CI90_Low = ci90_low,
    CI90_High = ci90_high,
    CI95_Low = ci95_low,
    CI95_High = ci95_high
  )
  
  return(results)
}

#apply this function to each imputed dataset
results_list <- lapply(datlist, function(data) two_stage_with_adjSEs(data))

#number of imputations
m <- length(datlist)

#pool the results from each imputed dataset using Rubin's rules 
pooled_estimate <- mean(sapply(results_list, function(x) x$Estimate))
within_var <- mean(sapply(results_list, function(x) x$Robust_SE)^2) * m/(m - 1)
between_var <- var(sapply(results_list, function(x) x$Estimate)) * (m - 1)/m
total_var <- within_var + (1 + 1/m) * between_var
pooled_se <- sqrt(total_var)

#calculate degrees of freedom using Barnard and Rubin's adjustment
df <- (m - 1) * (1 + (within_var / total_var))^2

#calculate adjusted one-tailed and two-tailed p-values
one_tailed_p <- pt(-abs(pooled_estimate / pooled_se), df) 
two_tailed_p <- 2 * one_tailed_p

#constructing pooled results
pooled_results <- data.frame(
  Estimate = pooled_estimate,
  se = pooled_se,
  z = pooled_estimate / pooled_se,
  two_tailed_p = two_tailed_p,
  one_tailed_p = one_tailed_p,
  CI90_low = pooled_estimate - qnorm(0.95) * pooled_se,
  CI90_high = pooled_estimate + qnorm(0.95) * pooled_se,
  CI95_low = pooled_estimate - qnorm(0.975) * pooled_se,
  CI95_high = pooled_estimate + qnorm(0.975) * pooled_se,
  ROPE_low = d_to_oddsratio(-0.25)
)

#compute additional statistics, exponentiate,
#and rename to enable merging with MR results
ests_3b_14 <- pooled_results %>%
  mutate(Estimate = exp(Estimate),
         CI90_low = exp(CI90_low),
         CI90_high = exp(CI90_high),
         CI95_low = exp(CI95_low),
         CI95_high = exp(CI95_high),
         outcome = "Diagnoses age 10-17")

#merge MR results for symptoms and diagnoses
ests3_14 <- ests_3a_14 %>%
  full_join(ests_3b_14) %>%
  mutate(adj="1-sample MR")

#add MR negative control results
ests3 <- ests_3a_8y %>%
  mutate(adj="Negative control MR") %>%
  full_join(ests3_14) 

#load observational results for depression
load(file="./data/aim1_results_tr.RData")

ests1 <- ests1 %>%
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
ests_dep <- ests3 %>%
  full_join(ests1)

#restrict to 2 decimal points for reporting in text
ests_dep_2dp <- ests_dep
ests_dep_2dp[] <- lapply(ests_dep_2dp, function(x) if(is.numeric(x)) round(x, 2) else x)
ests_dep_2dp

#plot results
p3a <- ests_dep %>%
  filter(outcome == "14-year symptoms" | outcome == "8-year symptoms") %>%
  mutate(adj=fct_relevel(adj,"Negative control MR","1-sample MR","Adjusted","Unadjusted")) %>%
  ggplot(mapping=aes(x=Estimate,y="",xmin=CI95_low,xmax=CI95_high,colour=adj,shape=adj)) +
  geom_vline(aes(xintercept=0), colour = "grey80", linetype = 2) + theme_light(base_size=14) +
  geom_vline(aes(xintercept=ROPE_low), colour = c("#D55E00","#0072B2","#0072B2","#D55E00"), linetype=2, size=0.5) + 
  geom_errorbar(size = 1.4, alpha=0.6, width = 0, position=position_dodge(0.55)) +
  geom_point(size=4.3, fill="white", alpha=1, position=position_dodge(0.55)) +
  geom_text(aes(label = "Smallest effect sizes of interest"), x=-0.12, y=1.740, col="black", size=4.9) +
  scale_colour_manual(values=c("#C9C9C9","#0072B2","#D55E00","#942A12")) +
  scale_shape_manual(values=c(18,17,15,16)) +
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
        legend.text = element_text(size = 13),
        strip.background = element_rect(colour="gray", fill="white")) + 
  xlab(label=expression("Standardised"~beta)) +
  ylab(label="Depressive\nsymptoms") +
  coord_cartesian(xlim = c(-0.2,0.1), clip='off') + 
  scale_x_continuous(breaks = c(-0.2, -0.1, 0.0, 0.1))

p3b<- ests_dep %>%
  filter(outcome == "Diagnoses age 10-17") %>%
  ggplot(mapping=aes(x=Estimate,y="",xmin=CI95_low,xmax=CI95_high,colour=adj,shape=adj)) +
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

patch1 <- p3a / p3b

patch2 <- patch1 + plot_annotation(tag_levels = "A") &
  plot_layout(guides = "collect") &
  theme(legend.position = "top",
        legend.box.margin = margin(t = 0, r = 0, b = 15, l = 0),
        legend.title = element_blank())

patch2

tiff("figures/dep_MR_tr.tiff", res=800, compression="lzw", unit="in",
     height=6, width=8)

patch2

dev.off()
