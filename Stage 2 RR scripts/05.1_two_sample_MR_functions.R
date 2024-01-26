#05.1_two_sample_MR_functions.R

#make function to read in data for outcomes:
read_my_outcome_data <- function(filename) {
  read_outcome_data(
    snps = exp_dat$SNP,
    filename = filename,
    sep = " ",
    snp_col = "ID",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ",
    pval_col = "P",
    samplesize_col = "N"
  )
}

#make function to harmonise exposure and outcome data
harmonise_outcomes <- function(outcome_data) {
  harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_data, action = 1)
}

#I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

#regression dilution function
regression_dilution <- function(data) {
  
  #calculate the regression dilution Isq. 
  BetaXG   = data$beta.exposure
  seBetaXG = data$se.exposure 
  seBetaYG <- data$se.outcome
  BXG  = abs(BetaXG) 
  
  Isq_unweighted <- Isq(BXG,seBetaXG) 
  Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG))
  
  #print unweighted and weighted
  print(Isq_unweighted)
  print(Isq_weighted)
}

simex_correction <- function(data) {
  #ensure all gene-exposure estimates are positive
  BYG = data$beta.outcome * sign(data$beta.exposure) 
  BXG = abs(data$beta.exposure)
  
  #MR-Egger regression (unweighted)
  Fit <- lm(BYG ~ BXG, x = TRUE, y = TRUE)
  
  #simulation extrapolation
  mod.sim <- simex(Fit, B = 10000, measurement.error = data$se.exposure, 
                   SIMEXvariable = "BXG", fitting.method = "quad", asymptotic = FALSE)
  summary_mod <- summary(mod.sim)
  
  #return corrected MR-Egger intercept and slope estimates
  estimates_df <- as.data.frame(summary_mod$coefficients$jackknife[c("(Intercept)", "BXG"), ])
  rownames(estimates_df) <- NULL
  estimates_df$LCI <- estimates_df$Estimate - 1.96*estimates_df$`Std. Error`
  estimates_df$UCI <- estimates_df$Estimate + 1.96*estimates_df$`Std. Error`
  estimates_df$Egger_estimate <- c("Intercept", "Slope")
  
  return(estimates_df)
}

#make function to run across outcomes
conmix <- function(data, outcome){
  
  #select data
  conmix_data <- data %>%
    dplyr::select(snps=SNP, betaX=beta.exposure, betaXse=se.exposure, exposure=id.exposure,
                  betaY=beta.outcome, betaYse=se.outcome, outcome=id.outcome,
                  effect_allele=effect_allele.exposure, other_allele=other_allele.exposure)
  
  #create MR object
  mr_input <- mr_input(bx = conmix_data$betaX, bxse = conmix_data$betaXse,
                       by = conmix_data$betaY, byse = conmix_data$betaYse,
                       exposure = "aam", outcome = outcome, snps = conmix_data$snps,
                       effect_allele = conmix_data$effect_allele,
                       other_allele = conmix_data$other_allele)
  
  #run contamination mixture method
  conmix_result <- mr_conmix(mr_input,psi=0,CIMin=NA,CIMax=NA,CIStep=0.01,alpha=0.05)
  
  return(conmix_result)
  
}

#make function for Steiger filtering
steiger_filter <- function(data){
  
  #sample size/units of exposures/outcomes
  data$samplesize.outcome <- 9832
  data$units.outcome <- "SD"
  data$samplesize.exposure <- 252000
  data$units.exposure <- "SD"
  
  #run Steiger filtering
  analysis_df <- subset(data, data$mr_keep == TRUE)
  steiger <- TwoSampleMR::steiger_filtering(analysis_df)
  
  #filter data to include only true variants
  filtered_data <- analysis_df[steiger$steiger_dir, ]
  
  #report the proportion of SNPs that are true
  false <- length(steiger$steiger_dir[steiger$steiger_dir == FALSE])
  true <- length(steiger$steiger_dir[steiger$steiger_dir == TRUE])
  percent <- (true / (true + false)) * 100
  
  #print statistics
  print(data.frame(true, false, percent))
  
  #return the filtered data
  return(filtered_data)
}

#Create function to update MR-Egger estimates with unweighted SIMEX estimates, 
#and include results from the contamination mixture method (conmix)
update_results <- function(results_df, simex_results, conmix_results) {
  #update MR-Egger estimates
  results_df <- results_df %>% 
    mutate(b = ifelse(method == "MR Egger", simex_results$Estimate[2],b),
           se = ifelse(method == "MR Egger", simex_results$`Std. Error`[2],se),
           pval = ifelse(method == "MR Egger", simex_results$`Pr(>|t|)`[2],pval),
           one_tailed_pval = pval/2,
           lci95 = b - se * 1.96,  #calculate 95% confidence intervals
           uci95 = b + se * 1.96)  
  
  #include results from the contamination mixture method
  conmix_row <- data.frame(
    id.exposure = "aam",
    id.outcome = results_df$id.outcome,
    outcome = results_df$outcome,
    exposure = "Age at menarche",
    method = "Contamination mixture method",
    nsnp = conmix_results@SNPs,
    b = conmix_results@Estimate,     
    se = NA,  
    lci95 = conmix_results@CILower,
    uci95 = conmix_results@CIUpper,
    pval = conmix_results@Pvalue,
    one_tailed_pval = conmix_results@Pvalue/2
  )
  
  #add conmix to results and return
  return(rbind(results_df, conmix_row))
}

#Create function for plotting results for different outcomes
plot_scatter <- function(results, data, base_size, legend_size) {
  p <- mr_scatter_plot(results, data)[[1]] +
    theme_minimal(base_size = base_size) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = legend_size),
          legend.position = "top",
          plot.background = element_rect(fill = "white", color = "white"))
  
  return(p)
}

#Make function for leave-one-out analyses
perform_loo_analysis <- function(data, filename_loo) {
  res_loo <- mr_leaveoneout(data)
  p2 <- mr_leaveoneout_plot(res_loo)
  ggsave(p2[[1]], file = paste0("./2-sample plots/",filename_loo), width = 7, height = 10)
}

#Make function for single-SNP analyses
perform_single_snp_analysis <- function(data, filename_singlesnp) {
  res_single <- mr_singlesnp(data, all_method = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  p3 <- mr_forest_plot(res_single)
  ggsave(p3[[1]], file = paste0("./2-sample plots/",filename_singlesnp), width = 7, height = 10)
}
