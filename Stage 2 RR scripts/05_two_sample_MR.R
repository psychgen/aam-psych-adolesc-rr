#05_two_sample_MR.R

##This script conducts two-sample MR analyses and sensitivity checks (outside of 
##TSD), sourcing the script '05.1_two_sample_MR_functions.R'. 
source("05.1_two_sample_MR_functions.R")

################################################################################
## Contents
# Reading in data
# 1. Clumping exposure data
# 2. Harmonising exposure and outcome data 
# 3. 2-sample MR analysis
# 4. MR Egger sensitivity checks
# 5. Tests for pleiotropy 
# 6. Test for reverse causation
# 7. Visual inspection
################################################################################

##set working directory
setwd("~/Desktop/RR")

##load packages
library(devtools)
library(TwoSampleMR)
library(ieugwasr)
library(knitr)
library(tidyverse)
library(data.table)
library(simex) 
library(remotes)
library(MendelianRandomization)
library(glmnet)
library(MASS)
library(quantreg)
library(patchwork)

################################################################################
##########                    Reading in data                      #############
################################################################################

## Example of the required column headers in the two-sample MR package: 
# Phenotype 	  SNP			    CHR 	BP 			    effect_allele 	other_allele 	eaf	 beta 	 se 		  pval
# Neuroticism 	rs4653651 	1 		225862060 	A 				      G 					  0.25 -0.091  0.0202 	6.443e-06

#read in exposure sumstats and make CSV
#sstats <- fread("Menarche_1KG_NatGen2017_WebsiteUpload_subset.txt",data.table=FALSE)
#write_csv(sstats, "aam.csv")

#read in data for exposure:
exp_dat <- read_exposure_data(
  filename = "aam.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "Pvalue"
  )

#calculate SE based on beta and p-value and rename exposure
exp_dat <- exp_dat %>%
  mutate(se.exposure = get_se(beta.exposure, pval.exposure),
         exposure = "age at menarche",
         id.exposure = "aam",
         mr_keep.exposure = "TRUE")

#apply function to read in data for different summary statistics
out_dep <- read_my_outcome_data("smfq_dep_c_14c.txt") %>%
  mutate(outcome="14-year depressive symptoms", id.outcome="smfq")

out_anx <- read_my_outcome_data("scared_anx_c_14c.txt") %>%
  mutate(outcome="anxiety symptoms", id.outcome="scared")

out_cd <- read_my_outcome_data("rsdbd_cd_c_14c.txt") %>%
  mutate(outcome="conduct symptoms", id.outcome="rsdbd")

out_odd <- read_my_outcome_data("rsdbd_odd_c_14m.txt") %>%
  mutate(outcome="oppositional symptoms", id.outcome="rsdbd")

out_adhd <- read_my_outcome_data("rsdbd_adhd_c_14m.txt") %>%
  mutate(outcome="ADHD traits", id.outcome="rsdbd")

################################################################################
###############         1. Clumping exposure data             ##################
################################################################################

#perform clumping locally using plink, and downloaded EUR reference panel
clump_exp_dat <- ld_clump(
  dplyr::tibble(rsid=exp_dat$SNP, pval=exp_dat$pval.exposure, id=exp_dat$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/Users/joadriandahlaskelund/Desktop/RR/EUR",
  clump_kb = 10000, 
  clump_r2 = 0.001,
  clump_p = 5e-8
)

#subset exp_dat to independent SNPs
exp_dat <- exp_dat %>%
  filter(SNP %in% clump_exp_dat$rsid)

################################################################################
#########      2. Harmonising exposure and outcome data             ############
################################################################################

#apply harmonisation function to different outcomes (action = 1)
dep_data <- harmonise_outcomes(out_dep)
anx_data <- harmonise_outcomes(out_anx)
cd_data <- harmonise_outcomes(out_cd)
odd_data <- harmonise_outcomes(out_odd)
adhd_data <- harmonise_outcomes(out_adhd)

################################################################################
###############              3. 2-sample MR analysis            ################
################################################################################

dep_results <- mr(dep_data,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
anx_results <- mr(anx_data,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
cd_results <- mr(cd_data,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
odd_results <- mr(odd_data,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
adhd_results <- mr(adhd_data,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))

################################################################################
#########              4. MR Egger sensitivity checks              #############
################################################################################

#MR Egger Intercept
#estimating the degree of bias due to directional horizontal pleiotropy
mr_egger_int_dep <- mr_pleiotropy_test(dep_data)#small but significant
mr_egger_int_anx <- mr_pleiotropy_test(anx_data)#not sig.
mr_egger_int_cd <- mr_pleiotropy_test(cd_data)#not sig.
mr_egger_int_odd <- mr_pleiotropy_test(odd_data)#not sig.
mr_egger_int_adhd <- mr_pleiotropy_test(adhd_data)#small but significant

#regression dilution I-squared:
#are the MR Egger estimates valid? 
#0.6 <Isq> 0.9 = perform SIMEX correction

#apply functions to different outcomes
regression_dilution(dep_data)
regression_dilution(anx_data)
regression_dilution(cd_data)
regression_dilution(odd_data)
regression_dilution(adhd_data)

##conclusion: we need to perform SIMEX corrections (unweighted)

#apply SIMEX correction function to different outcomes
egger_dep <- simex_correction(dep_data)
egger_anx <- simex_correction(anx_data)
egger_cd <- simex_correction(cd_data)
egger_odd <- simex_correction(odd_data)
egger_adhd <- simex_correction(adhd_data)

################################################################################
##########                5. Tests for Pleiotropy                   ############
################################################################################

##Estimating heterogeneity
mr_het_dep <- mr_heterogeneity(dep_data)
mr_het_anx <- mr_heterogeneity(anx_data)
mr_het_cd <- mr_heterogeneity(cd_data)
mr_het_odd <- mr_heterogeneity(odd_data)
mr_het_adhd <- mr_heterogeneity(adhd_data)

##MR PRESSO
run_mr_presso(dep_data, NbDistribution=1000, SignifThreshold=0.05)#no outliers
run_mr_presso(anx_data, NbDistribution=1000, SignifThreshold=0.05)#no outliers
run_mr_presso(cd_data, NbDistribution=1000, SignifThreshold=0.05)#no outliers
run_mr_presso(odd_data, NbDistribution=1000, SignifThreshold=0.05)#no outliers
run_mr_presso(adhd_data, NbDistribution=1000, SignifThreshold=0.05)#no outliers

##Contamination mixture method

#apply conmix functions to different outcomes
dep_conmix <- conmix(dep_data, dep_data$id.outcome[1])
anx_conmix <- conmix(anx_data, anx_data$id.outcome[1])
cd_conmix <- conmix(cd_data, cd_data$id.outcome[1])
odd_conmix <- conmix(odd_data, odd_data$id.outcome[1])
adhd_conmix <- conmix(adhd_data, adhd_data$id.outcome[1])

################################################################################
##########             6. Test for reverse causation                ############
################################################################################

#apply Steiger filtering function 
dep_steiger <- steiger_filter(dep_data)
anx_steiger <- steiger_filter(anx_data)
cd_steiger <- steiger_filter(cd_data)
odd_steiger <- steiger_filter(odd_data)
adhd_steiger <- steiger_filter(adhd_data)

#re-run 2-sample MR sensitivity analyses after Steiger filtering
dep_steiger_res <- mr(dep_steiger,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")) %>%
  mutate(lci = b - se * 1.96, uci = b + se * 1.96, one_tailed_p = pval/2) #95% CIs and one-tailed pval
anx_steiger_res <- mr(anx_steiger,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")) %>%
  mutate(lci = b - se * 1.96, uci = b + se * 1.96) #calculate 95% confidence intervals
cd_steiger_res <- mr(cd_steiger,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")) %>%
  mutate(lci = b - se * 1.96, uci = b + se * 1.96) #calculate 95% confidence intervals
odd_steiger_res <- mr(odd_steiger,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")) %>%
  mutate(lci = b - se * 1.96, uci = b + se * 1.96) #calculate 95% confidence intervals
adhd_steiger_res <- mr(adhd_steiger,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")) %>%
  mutate(lci = b - se * 1.96, uci = b + se * 1.96) #calculate 95% confidence intervals

#re-run contamination mixture method after Steiger filtering
dep_steiger_conmix <- conmix(dep_steiger, dep_steiger$id.outcome[1])
anx_steiger_conmix <- conmix(anx_steiger, anx_steiger$id.outcome[1])
cd_steiger_conmix <- conmix(cd_steiger, cd_steiger$id.outcome[1])
odd_steiger_conmix <- conmix(odd_steiger, odd_steiger$id.outcome[1])
adhd_steiger_conmix <- conmix(adhd_steiger, adhd_steiger$id.outcome[1])

################################################################################
##########                  7. Visual Inspection                    ############
################################################################################

#create lists of datasets and their corresponding names and results
result_sets <- list(dep_results, anx_results, cd_results, odd_results, adhd_results)
simex_sets <- list(egger_dep, egger_anx, egger_cd, egger_odd, egger_adhd)
conmix_sets <- list(dep_conmix, anx_conmix, cd_conmix, odd_conmix, adhd_conmix)
data_sets <- list(dep_data, anx_data, cd_data, odd_data, adhd_data)
names_sets <- c("dep", "anx", "cd", "odd", "adhd")

#iterate over the datasets, updating the results, and making plots
other_plots <- list()

for (i in seq_along(result_sets)) {
  updated_results <- update_results(result_sets[[i]], simex_sets[[i]], conmix_sets[[i]])
  
  print(updated_results)
  
  base_size <- if (names_sets[i] != "dep") 14 else 15
  legend_size <- if (names_sets[i] != "dep") 20 else 16
  
  #generate scatterplot
  p <- plot_scatter(updated_results, data_sets[[i]], base_size, legend_size)
  
  #save the "dep" plot separately
  if (names_sets[i] == "dep") {
    ggsave(p, file = paste0("./2-sample plots/scatter_", names_sets[i], ".png"), width = 7, height = 7)
  } else {
    #collect the other plots
    other_plots[[names_sets[i]]] <- p
  }
  
  #perform and save leave-one-out analysis plots
  loo_filename <- paste0("loo_", names_sets[i], ".png")
  perform_loo_analysis(data_sets[[i]], loo_filename)
  
  #perform and save single SNP analysis plots
  singlesnp_filename <- paste0("single_forest_", names_sets[i], ".png")
  perform_single_snp_analysis(data_sets[[i]], singlesnp_filename)
}

#combine plots for other symptom domains
combined_plot <- wrap_plots(other_plots, ncol = 2, guides = "collect") &
  theme(legend.position = "top")

#save out combined plot
ggsave(combined_plot, filename = "./2-sample plots/combined_scatter_plots.png", width = 10, height = 9)
