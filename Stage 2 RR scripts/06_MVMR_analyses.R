#06_MVMR_analyses.R

##This script conducts MVMR analyses and sensitivity checks (outside of TSD 
#based only on summary statistics), sourcing the script '06.1_MVMR_functions.R'. 
source("06.1_MVMR_functions.R")

##set working directory
setwd("~/Desktop/RR")

##load packages
library(TwoSampleMR)
library(MVMR)
library(tidyverse)
library(data.table)
library(ieugwasr)

################################################################################
##########                 MVMR (childhood body size)              #############
################################################################################

##read in and process adult BMI sumstats
#sstats_adult_BMI <- fread("BMI_adult_female.txt",data.table=FALSE) %>%
# filter(nchar(ALLELE1) == 1 & ALLELE1 %in% c("A", "T", "C", "G")) %>%
# filter(nchar(ALLELE0) == 1 & ALLELE0 %in% c("A", "T", "C", "G"))
#
##write out
#write_csv(sstats_adult_BMI, "BMI_adult_female.csv")

##read in and process childhood body size sumstats
#sstats_child_BMI <- fread("BMI_10_female.txt",data.table=FALSE) %>%
#  filter(nchar(ALLELE1) == 1 & ALLELE1 %in% c("A", "T", "C", "G")) %>%
#  filter(nchar(ALLELE0) == 1 & ALLELE0 %in% c("A", "T", "C", "G"))
#
##write out
#write_csv(sstats_child_BMI, "BMI_10_female.csv")

#read in data for the main exposure:
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

#read in a second exposure (child BMI)
exp2_dat <- read_exposure_data(
  filename = "BMI_10_female.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_LINREG")

#label exposures 1 and 2
exp_dat$id.exposure<- 1
exp2_dat$id.exposure<- 2
exp2_dat$exposure<- "Child BMI"

#merge the dataframes
exposures <- rbind(exp_dat, exp2_dat)

#check for overlapping SNPs
n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq >1, ]
exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq >1], ]

#extract both sets of instruments from the exposure 1 GWAS
exp1_MV <- read_outcome_data(snps = exposures$SNP, 
                             filename = "aam.csv", 
                             sep = ",",
                             snp_col = "SNP",
                             beta_col = "Effect",
                             effect_allele_col = "Allele1",
                             other_allele_col = "Allele2",
                             pval_col = "Pvalue")
exp1_MV$id.outcome <- 1
exp1_MV$Phenotype<- "exposure1"

#calculate SE based on beta and p-value
exp1_MV <- exp1_MV %>%
  mutate(se.outcome = get_se(beta.outcome, pval.outcome),
         mr_keep.outcome = "TRUE")

#extract the instruments from the second exposure
exp2_MV<- read_outcome_data(snps = exposures$SNP, 
                            filename = "BMI_10_female.csv", 
                            sep = ",",
                            snp_col = "SNP",
                            beta_col = "BETA",
                            effect_allele_col = "ALLELE1",
                            other_allele_col = "ALLELE0",
                            pval_col = "P_LINREG")
exp2_MV$id.outcome<- 2
exp2_MV$Phenotype<- "exposure2"

#calculate SE based on beta and p-value
exp2_MV <- exp2_MV %>%
  mutate(se.outcome = get_se(beta.outcome, pval.outcome),
         mr_keep.outcome = "TRUE")

#merge them to extract from outcome
exposures <- rbind(exp1_MV, exp2_MV)

#turn outcome to exposure to clump
names(exposures) <- gsub("outcome", "exposure", names(exposures))

exposures$id.exposure[exposures$id.exposure == "2"] <- "1"

#perform clumping locally using plink, and downloaded EUR reference panel
clump_exposures <- ld_clump(
  dplyr::tibble(rsid=exposures$SNP, pval=exposures$pval.exposure, id=exposures$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/Users/joadriandahlaskelund/Desktop/RR/EUR",
  clump_kb = 10000, 
  clump_r2 = 0.001
)

#subset exposures to independent SNPs
exposures <- exposures %>%
  filter(SNP %in% clump_exposures$rsid)

#keep only snps that are present across both exposures
n_occur <- data.frame(table(exposures$SNP))
n_occur[n_occur$Freq == 2,]
exposures<- exposures[exposures$SNP %in% n_occur$Var1[n_occur$Freq == 2],]

#split again to harmonise based on exposure id
exposures$id.exposure[exposures$Phenotype=="exposure1"] <- 1
exposures$id.exposure[exposures$Phenotype=="exposure2"] <- 2
exp1 = split(exposures, exposures$id.exposure)[['1']]
exp2 = split(exposures, exposures$id.exposure)[['2']]

#harmonise exposure1 on exposure 2
#to ensure that the SNPs are aligned correctly
names(exp2) <- gsub("exposure", "outcome", names(exp2))
exposures<- harmonise_data(exp1, exp2, action = 1)
#keep only snps MrKeep= TRUE
exposures<- exposures[exposures$mr_keep== TRUE, ]

#split the tables
exp1_H<- subset(exposures, id.exposure== id.exposure[1], 
                select= c(SNP, exposure, id.exposure, effect_allele.exposure, 
                          other_allele.exposure, beta.exposure, se.exposure, pval.exposure))
#split the tables
exp2_H<- subset(exposures, id.outcome== id.outcome[2], 
                select= c(SNP, outcome, id.outcome, effect_allele.outcome, 
                          other_allele.outcome, beta.outcome, se.outcome, pval.outcome))

#turn to exposure to merge the datasets
names(exp2_H) <- gsub("outcome", "exposure", names(exp2_H))

#merge the harmonised datasets
Exposures_H<- rbind(exp1_H, exp2_H)
Exposures_H["Phenotype"]<- NA
Exposures_H$Phenotype[Exposures_H$id.exposure == 1] <- "exposure1"
Exposures_H$Phenotype[Exposures_H$id.exposure == 2] <- "exposure2"
Exposures_H["eaf.exposure"]<- NA

#extract from the outcome
outcome<- read_outcome_data(snps = Exposures_H$SNP,
                            filename = "smfq_dep_c_14c.txt",
                            sep = " ",
                            snp_col = "ID",
                            beta_col = "BETA",
                            se_col = "SE",
                            effect_allele_col = "ALLELE1",
                            other_allele_col = "ALLELE0",
                            pval_col = "P")

#harmonise the exposures and outcome
mv_data <- harmonise_data(Exposures_H, outcome, action = 1)

#keep only mrkeep= TRUE
mv_data<- mv_data[mv_data$mr_keep== TRUE, ]

#format to analyse with MVMR package
bX1<- c(mv_data$beta.exposure[mv_data$id.exposure== 1])
bX2<- c(mv_data$beta.exposure[mv_data$id.exposure== 2])
bY<- c(mv_data$beta.outcome[mv_data$id.exposure== 1])
bYse<- c(mv_data$se.outcome[mv_data$id.exposure== 1])
bXse1<- c(mv_data$se.exposure[mv_data$id.exposure== 1])
bXse2<- c(mv_data$se.exposure[mv_data$id.exposure== 2])

#merge them into one dataframe
df<- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)
df_mvmr<- format_mvmr(df[, c(1, 3)], df[,5], df[,c(2, 4)], df[, 6])

#run the MVMR analysis
res<- ivw_mvmr(df_mvmr, gencov=0) %>%
  as.data.frame() %>%
  mutate(lci95 = Estimate - `Std. Error`*1.96,
         uci95 = Estimate + `Std. Error`*1.96,
         p_one_tailed=`Pr(>|t|)`/2)

#print mvmr results
print(res)

#F-statistic estimation
Fstat<- strength_mvmr(df_mvmr, gencov=0)

##############################################################################
#########       MVMR sensitivity analyses (childhood body size)   ############
##############################################################################

#modified Q-statistic for horisontal pleiotropy
ptr<-pleiotropy_mvmr(df_mvmr, gencov=0)

#reformat data for MVMR-Egger sensitivity analysis
mvmr_egger_input <- mr_mvinput(
  bx = cbind(bX1, bX2),
  by = bY,
  bxse = cbind(bXse1, bXse2),
  byse = bYse
)

#MVMR-Egger
mvmr_egger <- mr_mvegger(
  mvmr_egger_input,
  orientate = 1,
  alpha = 0.05
)

#print results
print(mvmr_egger)

#reformat data for further MVMR sensitivity analyses
bx = as.matrix(df_mvmr[, c("betaX1", "betaX2")])  #exposure effect sizes
sebx = as.matrix(df_mvmr[, c("sebetaX1", "sebetaX2")])  #exposure ses
by = df_mvmr$betaYG  #outcome effect sizes
seby = df_mvmr$sebetaYG  #outcome standard errors

# Identify rows with INF values
inf_rows = apply(sebx, 1, function(x) any(is.infinite(x)))

# Filter out these rows from all inputs
bx <- bx[!inf_rows, ]
sebx <- sebx[!inf_rows, ]
by <- by[!inf_rows]
seby <- seby[!inf_rows]

#use mvmr_median function to obtain results
mvmr_median_results = mvmr_median(bx, sebx, by, seby, boot_it = 1000)

#format results and calculate CIs
mvmr_median <- mvmr_median_results %>%
  as.data.frame() %>%
  mutate(lci95 = coefficients - se*1.96,
         uci95 = coefficients + se*1.96)
#print
print(mvmr_median)

#use mvmr_lasso function to run lasso
mvmr_lasso_results = mvmr_lasso(bx, by, seby)

#format results and calculate CIs
mvmr_lasso_ests <- as.data.frame(mvmr_lasso_results$th_post)
colnames(mvmr_lasso_ests)[1] <- "coefficients"

mvmr_lasso_ses <- as.data.frame(mvmr_lasso_results$se_post)
colnames(mvmr_lasso_ses)[1] <- "se"

mvmr_lasso <- mvmr_lasso_ests %>%
  cbind(mvmr_lasso_ses) %>%
  mutate(lci95 = coefficients - se*1.96,
         uci95 = coefficients + se*1.96)
#print
print(mvmr_lasso)

################################################################################
#############                  MVMR (adult BMI)                #################
################################################################################

#read in the other exposure (adult BMI)
exp4_dat <- read_exposure_data(
  filename = "BMI_adult_female.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P_LINREG")

#label exposures 1 and 2
exp_dat$id.exposure<- 1
exp4_dat$id.exposure<- 2
exp4_dat$exposure<- "Adult BMI"

#merge the dataframes
exposures1 <- rbind(exp_dat, exp4_dat)

#check for overlapping SNPs
n_occur1 <- data.frame(table(exposures1$SNP))
n_occur1[n_occur1$Freq >1, ]
exposures[exposures$SNP %in% n_occur1$Var1[n_occur1$Freq >1], ]

#extract both sets of instruments from the exposure 1 GWAS
exp3_MV <- read_outcome_data(snps = exposures1$SNP, 
                             filename = "aam.csv", 
                             sep = ",",
                             snp_col = "SNP",
                             beta_col = "Effect",
                             effect_allele_col = "Allele1",
                             other_allele_col = "Allele2",
                             pval_col = "Pvalue")
exp3_MV$id.outcome <- 1
exp3_MV$Phenotype<- "exposure1"

#calculate SE based on beta and p-value
exp3_MV <- exp3_MV %>%
  mutate(se.outcome = get_se(beta.outcome, pval.outcome),
         mr_keep.outcome = "TRUE")

#extract the instruments from the second exposure
exp4_MV<- read_outcome_data(snps = exposures1$SNP, 
                            filename = "BMI_adult_female.csv", 
                            sep = ",",
                            snp_col = "SNP",
                            beta_col = "BETA",
                            effect_allele_col = "ALLELE1",
                            other_allele_col = "ALLELE0",
                            pval_col = "P_LINREG")
exp4_MV$id.outcome<- 2
exp4_MV$Phenotype<- "exposure2"

#calculate SE based on beta and p-value
exp4_MV <- exp4_MV %>%
  mutate(se.outcome = get_se(beta.outcome, pval.outcome),
         mr_keep.outcome = "TRUE")

#merge them to extract from outcome
exposures2 <- rbind(exp3_MV, exp4_MV)

#turn outcome to exposure to clump
names(exposures2) <- gsub("outcome", "exposure", names(exposures2))

#clump the data
exposures2$id.exposure[exposures2$id.exposure == "2"] <- "1"

#perform clumping locally using plink, and EUR reference panel
clump_exposures2 <- ld_clump(
  dplyr::tibble(rsid=exposures2$SNP, pval=exposures2$pval.exposure, id=exposures2$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/Users/joadriandahlaskelund/Desktop/RR/EUR",
  clump_kb = 10000, 
  clump_r2 = 0.001
)

#subset exposures to independent SNPs
exposures2 <- exposures2 %>%
  filter(SNP %in% clump_exposures2$rsid)

#keep only snps that are present across both exposures
n_occur2 <- data.frame(table(exposures2$SNP))
n_occur2[n_occur2$Freq == 2,]
exposures2<- exposures2[exposures2$SNP %in% n_occur2$Var1[n_occur2$Freq == 2],]

#split again to harmonise based on exposure id
exposures2$id.exposure[exposures2$Phenotype=="exposure1"] <- 1
exposures2$id.exposure[exposures2$Phenotype=="exposure2"] <- 2
exp3 = split(exposures2, exposures2$id.exposure)[['1']]
exp4 = split(exposures2, exposures2$id.exposure)[['2']]

#Harmonise exposure 3 on exposure 4,
#to ensure that the SNPs are aligned correctly
names(exp4) <- gsub("exposure", "outcome", names(exp4))
exposures2<- harmonise_data(exp3, exp4, action = 1)
#keep only snps MrKeep=TRUE
exposures2<- exposures2[exposures2$mr_keep== TRUE, ]

#split the tables
exp3_H<- subset(exposures2, id.exposure== id.exposure[1], 
                select= c(SNP, exposure, id.exposure, effect_allele.exposure, 
                          other_allele.exposure, beta.exposure, se.exposure, pval.exposure))
#split the tables
exp4_H<- subset(exposures2, id.outcome== id.outcome[2], 
                select= c(SNP, outcome, id.outcome, effect_allele.outcome, 
                          other_allele.outcome, beta.outcome, se.outcome, pval.outcome))

#turn to exposure to merge the datasets
names(exp4_H) <- gsub("outcome", "exposure", names(exp4_H))

#merge the harmonised datasets
Exposures_H2<- rbind(exp3_H, exp4_H)
Exposures_H2["Phenotype"]<- NA
Exposures_H2$Phenotype[Exposures_H2$id.exposure == 1] <- "exposure1"
Exposures_H2$Phenotype[Exposures_H2$id.exposure == 2] <- "exposure2"
Exposures_H2["eaf.exposure"]<- NA

#extract from the outcome
outcome2<- read_outcome_data(snps = Exposures_H2$SNP,
                             filename = "smfq_dep_c_14c.txt",
                             sep = " ",
                             snp_col = "ID",
                             beta_col = "BETA",
                             se_col = "SE",
                             effect_allele_col = "ALLELE1",
                             other_allele_col = "ALLELE0",
                             pval_col = "P")

#harmonise the exposures and outcome
mv_data2 <- harmonise_data(Exposures_H2, outcome2, action = 1)

#keep only mrkeep= TRUE
mv_data2<- mv_data2[mv_data2$mr_keep== TRUE, ]

#format to analyse with MVMR package
bX1_2<- c(mv_data2$beta.exposure[mv_data2$id.exposure== 1])
bX2_2<- c(mv_data2$beta.exposure[mv_data2$id.exposure== 2])
bY_2<- c(mv_data2$beta.outcome[mv_data2$id.exposure== 1])
bYse_2<- c(mv_data2$se.outcome[mv_data2$id.exposure== 1])
bXse1_2<- c(mv_data2$se.exposure[mv_data2$id.exposure== 1])
bXse2_2<- c(mv_data2$se.exposure[mv_data2$id.exposure== 2])

#merge them into one dataframe
df2<- data.frame(bX1_2, bXse1_2, bX2_2, bXse2_2, bY_2, bYse_2)
df_mvmr2<- format_mvmr(df2[, c(1, 3)], df2[,5], df2[,c(2, 4)], df2[, 6])

#run the analysis
res2<- ivw_mvmr(df_mvmr2) %>%
  as.data.frame() %>%
  mutate(lci95 = Estimate - `Std. Error`*1.96,
         uci95 = Estimate + `Std. Error`*1.96,
         p_one_tailed=`Pr(>|t|)`/2)

#print results
print(res2)

#conditional F-statistic estimation
Fstat2<- strength_mvmr(df_mvmr2, gencov=0)

##############################################################################
##############       MVMR sensitivity analyses (adult BMI)   #################
##############################################################################

#modified Q-statistic for horisontal pleiotropy
ptr2<-pleiotropy_mvmr(df_mvmr2, gencov=0)

#reformat data for MVMR-Egger sensitivity analysis
mvmr_egger_input2 <- mr_mvinput(
  bx = cbind(bX1_2, bX2_2),
  by = bY_2,
  bxse = cbind(bXse1_2, bXse2_2),
  byse = bYse_2
)

#MVMR-Egger
mvmr_egger2 <- mr_mvegger(
  mvmr_egger_input2,
  orientate = 1,
  alpha = 0.05
)

#print results
print(mvmr_egger2)

#reformat data for further MVMR sensitivity analyses
bx2 = as.matrix(df_mvmr2[, c("betaX1", "betaX2")])  # Exposure effect sizes
sebx2 = as.matrix(df_mvmr2[, c("sebetaX1", "sebetaX2")])  # Exposure standard errors
by2 = df_mvmr2$betaYG  # Outcome effect sizes
seby2 = df_mvmr2$sebetaYG  # Outcome standard errors

# Identify rows with INF values in any column of sebx2
inf_rows = apply(sebx2, 1, function(x) any(is.infinite(x)))

# Filter out these rows from all inputs
bx2 <- bx2[!inf_rows, ]
sebx2 <- sebx2[!inf_rows, ]
by2 <- by2[!inf_rows]
seby2 <- seby2[!inf_rows]

#use mvmr_median function to obtain results
mvmr_median_results2 = mvmr_median2(bx2, sebx2, by2, seby2, boot_it = 1000)

#format results and calculate CIs
mvmr_median2 <- mvmr_median_results2 %>%
  as.data.frame() %>%
  mutate(lci95 = coefficients - se*1.96,
         uci95 = coefficients + se*1.96)
#print
print(mvmr_median2)

#use mvmr_lasso function to run lasso
mvmr_lasso_results2 = mvmr_lasso2(bx2, by2, seby2)

#format results (estimates and standard errors)
mvmr_lasso_ests2 <- as.data.frame(mvmr_lasso_results2$th_post)
colnames(mvmr_lasso_ests2)[1] <- "coefficients"

mvmr_lasso_ses2 <- as.data.frame(mvmr_lasso_results2$se_post)
colnames(mvmr_lasso_ses2)[1] <- "se"

#merge estimates and SEs, and calculate CIs
mvmr_lasso2 <- mvmr_lasso_ests2 %>%
  cbind(mvmr_lasso_ses2) %>%
  mutate(lci95 = coefficients - se*1.96,
         uci95 = coefficients + se*1.96)

#print
print(mvmr_lasso2)

################################################################################
##########                    MVMR (estradiol UKB)                 #############
################################################################################

#read in and process sumstats, removing INDELS
#sstats_estradiol <- fread("GCST90020092_buildGRCh38.tsv", data.table=FALSE) %>%
#  select(-c(odds_ration,ci_lower,ci_upper)) %>%
#  filter(nchar(effect_allele) == 1 & effect_allele %in% c("A", "T", "C", "G")) %>%
#  filter(nchar(other_allele) == 1 & other_allele %in% c("A", "T", "C", "G"))
#
#write_csv(sstats_estradiol, file="estradiol_female.csv")

#Read in additional exposure (estradiol UKB females)
exp5_dat <- read_exposure_data(
  filename = "estradiol_female.csv",
  sep = ",",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value")

#label exposures 1 and 2
exp_dat$id.exposure<- 1
exp5_dat$id.exposure<- 2
exp5_dat$exposure<- "Estradiol"

#Merge the dataframes
exposures3 <- rbind(exp_dat, exp5_dat)

#Check for overlapping SNPs
n_occur3 <- data.frame(table(exposures3$SNP))
n_occur3[n_occur3$Freq >1, ]
exposures3[exposures3$SNP %in% n_occur3$Var1[n_occur3$Freq >1], ]

#Extract both sets of instruments from the exposure 1 GWAS
exp6_MV <- read_outcome_data(snps = exposures3$SNP, 
                             filename = "aam.csv", 
                             sep = ",",
                             snp_col = "SNP",
                             beta_col = "Effect",
                             effect_allele_col = "Allele1",
                             other_allele_col = "Allele2",
                             pval_col = "Pvalue")
exp6_MV$id.outcome <- 1
exp6_MV$Phenotype<- "exposure1"

#calculate SE based on beta and p-value
exp6_MV <- exp6_MV %>%
  mutate(se.outcome = get_se(beta.outcome, pval.outcome),
         mr_keep.outcome = "TRUE")

#Extract the instruments from the second exposure
exp7_MV<- read_outcome_data(snps = exposures3$SNP, 
                            filename = "estradiol_female.csv",
                            sep = ",",
                            snp_col = "variant_id",
                            beta_col = "beta",
                            se_col = "standard_error",
                            effect_allele_col = "effect_allele",
                            other_allele_col = "other_allele",
                            eaf_col = "effect_allele_frequency",
                            pval_col = "p_value")
exp7_MV$id.outcome<- 2
exp7_MV$Phenotype<- "exposure2"

#calculate SE based on beta and p-value
exp7_MV <- exp7_MV %>%
  mutate(se.outcome = get_se(beta.outcome, pval.outcome),
         mr_keep.outcome = "TRUE")

#Merge them to extract from outcome
exposures4 <- rbind(exp6_MV, exp7_MV)

#Turn outcome to exposure to clump
names(exposures4) <- gsub("outcome", "exposure", names(exposures4))

#Clump the data
exposures4$id.exposure[exposures4$id.exposure == "2"] <- "1"

#perform clumping locally using plink, and downloaded EUR reference panel
clump_exposures4 <- ld_clump(
  dplyr::tibble(rsid=exposures4$SNP, pval=exposures4$pval.exposure, id=exposures4$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/Users/joadriandahlaskelund/Desktop/RR/EUR",
  clump_kb = 10000, 
  clump_r2 = 0.001
)

#subset exposures to independent SNPs
exposures4 <- exposures4 %>%
  filter(SNP %in% clump_exposures4$rsid)

#Keep only snps that are present across both exposures
n_occur4 <- data.frame(table(exposures4$SNP))
n_occur4[n_occur4$Freq == 2,]
exposures4<- exposures4[exposures4$SNP %in% n_occur4$Var1[n_occur4$Freq == 2],]

#Split again to harmonise based on exposure id
exposures4$id.exposure[exposures4$Phenotype=="exposure1"] <- 1
exposures4$id.exposure[exposures4$Phenotype=="exposure2"] <- 2
exp5 = split(exposures4, exposures4$id.exposure)[['1']]
exp6 = split(exposures4, exposures4$id.exposure)[['2']]

#Harmonise exposure 5 on exposure 6
#To ensure that the SNPs are aligned correctly
names(exp6) <- gsub("exposure", "outcome", names(exp6))
exposures4<- harmonise_data(exp5, exp6, action = 1)
#Keep only snps MrKeep= TRUE
exposures4<- exposures4[exposures4$mr_keep== TRUE, ]

#Split the tables
exp5_H<- subset(exposures4, id.exposure== id.exposure[1], 
                select= c(SNP, exposure, id.exposure, effect_allele.exposure, 
                          other_allele.exposure, beta.exposure, se.exposure, pval.exposure))
#Split the tables
exp6_H<- subset(exposures4, id.outcome== id.outcome[2], 
                select= c(SNP, outcome, id.outcome, effect_allele.outcome, 
                          other_allele.outcome, beta.outcome, se.outcome, pval.outcome))

#Turn to exposure to merge the datasets
names(exp6_H) <- gsub("outcome", "exposure", names(exp6_H))

#Merge the harmonised datasets
Exposures_H3<- rbind(exp5_H, exp6_H)
Exposures_H3["Phenotype"]<- NA
Exposures_H3$Phenotype[Exposures_H3$id.exposure == 1] <- "exposure1"
Exposures_H3$Phenotype[Exposures_H3$id.exposure == 2] <- "exposure2"
Exposures_H3["eaf.exposure"]<- NA

#Extract from the outcome
outcome3<- read_outcome_data(snps = Exposures_H3$SNP,
                             filename = "smfq_dep_c_14c.txt",
                             sep = " ",
                             snp_col = "ID",
                             beta_col = "BETA",
                             se_col = "SE",
                             effect_allele_col = "ALLELE1",
                             other_allele_col = "ALLELE0",
                             pval_col = "P")

#harmonise the exposures and outcome
mv_data3 <- harmonise_data(Exposures_H3, outcome3, action = 1)

#Keep only mrkeep= TRUE
mv_data3<- mv_data3[mv_data3$mr_keep== TRUE, ]

#format to analyse with MVMR package
bX1_3<- c(mv_data3$beta.exposure[mv_data3$id.exposure== 1])
bX2_3<- c(mv_data3$beta.exposure[mv_data3$id.exposure== 2])
bY_3<- c(mv_data3$beta.outcome[mv_data3$id.exposure== 1])
bYse_3<- c(mv_data3$se.outcome[mv_data3$id.exposure== 1])
bXse1_3<- c(mv_data3$se.exposure[mv_data3$id.exposure== 1])
bXse2_3<- c(mv_data3$se.exposure[mv_data3$id.exposure== 2])

#Merge them into one dataframe
df3<- data.frame(bX1_3, bXse1_3, bX2_3, bXse2_3, bY_3, bYse_3)

df_mvmr3<- format_mvmr(df3[, c(1, 3)], df3[,5], df3[,c(2, 4)], df3[, 6])

#run the analysis
res3<- ivw_mvmr(df_mvmr3) %>%
  as.data.frame() %>%
  mutate(lci95 = Estimate - `Std. Error`*1.96,
         uci95 = Estimate + `Std. Error`*1.96,
         p_one_tailed = `Pr(>|t|)`/2)

#print the results
print(res3)

#F-statistic estimation
Fstat3<- strength_mvmr(df_mvmr3, gencov=0)

#test for horizontal pleiotropy
ptr3<-pleiotropy_mvmr(df_mvmr3, gencov=0)

##############################################################################
############       MVMR sensitivity analyses (estradiol UKB)   ###############
##############################################################################

#reformat data for MVMR-Egger sensitivity analysis
mvmr_egger_input3 <- mr_mvinput(
  bx = cbind(bX1_3, bX2_3),
  by = bY_3,
  bxse = cbind(bXse1_3, bXse2_3),
  byse = bYse_3
)

#MVMR-Egger
mvmr_egger3 <- mr_mvegger(
  mvmr_egger_input3,
  orientate = 1,
  alpha = 0.05
)

#print results
print(mvmr_egger3)

#reformat data for further MVMR sensitivity analyses
bx3 = as.matrix(df_mvmr3[, c("betaX1", "betaX2")])  # Exposure effect sizes
sebx3 = as.matrix(df_mvmr3[, c("sebetaX1", "sebetaX2")])  # Exposure standard errors
by3 = df_mvmr3$betaYG  # Outcome effect sizes
seby3 = df_mvmr3$sebetaYG  # Outcome standard errors

#use mvmr_median function to obtain results
mvmr_median_results3 = mvmr_median(bx3, sebx3, by3, seby3, boot_it = 1000)

#format results and calculate CIs
mvmr_median3 <- mvmr_median_results3 %>%
  as.data.frame() %>%
  mutate(lci95 = coefficients - se*1.96,
         uci95 = coefficients + se*1.96)
#print
print(mvmr_median3)

#use mvmr_lasso function to run lasso
mvmr_lasso_results3 = mvmr_lasso(bx3, by3, seby3)

#format results (estimates and standard errors)
mvmr_lasso_ests3 <- as.data.frame(mvmr_lasso_results3$th_post)
colnames(mvmr_lasso_ests3)[1] <- "coefficients"

mvmr_lasso_ses3 <- as.data.frame(mvmr_lasso_results3$se_post)
colnames(mvmr_lasso_ses3)[1] <- "se"

#merge estimates and SEs, and calculate CIs
mvmr_lasso3 <- mvmr_lasso_ests3 %>%
  cbind(mvmr_lasso_ses3) %>%
  mutate(lci95 = coefficients - se*1.96,
         uci95 = coefficients + se*1.96)

#print
print(mvmr_lasso3)

################################################################################
##########                  MVMR (major depression)               ##############
################################################################################

#read in and process sumstats
#sstats_dep <- fread("PGC_UKB_depression_genome-wide.txt",data.table=FALSE)
#write_csv(sstats_dep, "dep.csv")

#Read in second exposure (depression PGC + UKB)
exp7_dat <- read_exposure_data(
  filename = "dep.csv",
  sep = ",",
  snp_col = "MarkerName",
  beta_col = "LogOR",
  se_col = "StdErrLogOR",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "Freq",
  pval_col = "P")

#label exposures 1 and 2
exp_dat$id.exposure<- 1
exp7_dat$id.exposure<- 2
exp7_dat$exposure<- "Depression"

#merge exposure dataframes
exposures5 <- rbind(exp_dat, exp7_dat)

#check for overlapping SNPs
n_occur5 <- data.frame(table(exposures5$SNP))
n_occur5[n_occur5$Freq >1, ]
exposures5[exposures5$SNP %in% n_occur5$Var1[n_occur5$Freq >1], ]

#extract instruments from the first exposure GWAS
exp8_MV <- read_outcome_data(snps = exposures5$SNP, 
                             filename = "aam.csv", 
                             sep = ",",
                             snp_col = "SNP",
                             beta_col = "Effect",
                             effect_allele_col = "Allele1",
                             other_allele_col = "Allele2",
                             pval_col = "Pvalue")
exp8_MV$id.outcome <- 1
exp8_MV$Phenotype<- "exposure1"

#calculate SE based on beta and p-value
exp8_MV <- exp8_MV %>%
  mutate(se.outcome = get_se(beta.outcome, pval.outcome),
         mr_keep.outcome = "TRUE")

#extract instruments from the second exposure
exp9_MV<- read_outcome_data(snps = exposures5$SNP, 
                            filename = "dep.csv",
                            sep = ",",
                            snp_col = "MarkerName",
                            beta_col = "LogOR",
                            se_col = "StdErrLogOR",
                            effect_allele_col = "A1",
                            other_allele_col = "A2",
                            eaf_col = "Freq",
                            pval_col = "P")
exp9_MV$id.outcome<- 2
exp9_MV$Phenotype<- "exposure2"

#calculate SE based on beta and p-value
exp9_MV <- exp9_MV %>%
  mutate(se.outcome = get_se(beta.outcome, pval.outcome),
         mr_keep.outcome = "TRUE")

#merge to extract from the outcome
exposures6 <- rbind(exp8_MV, exp9_MV)

#turn outcome to exposure to clump
names(exposures6) <- gsub("outcome", "exposure", names(exposures6))

exposures6$id.exposure[exposures6$id.exposure == "2"] <- "1"

#perform clumping using plink
clump_exposures6 <- ld_clump(
  dplyr::tibble(rsid=exposures6$SNP, pval=exposures6$pval.exposure, id=exposures6$id.exposure),
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = "/Users/joadriandahlaskelund/Desktop/RR/EUR",
  clump_kb = 10000, 
  clump_r2 = 0.001
)

#subset exposures to independent SNPs
exposures6 <- exposures6 %>%
  filter(SNP %in% clump_exposures6$rsid)

#keep only snps that are present across both exposures
n_occur6 <- data.frame(table(exposures6$SNP))
n_occur6[n_occur6$Freq == 2,]
exposures6<- exposures6[exposures6$SNP %in% n_occur6$Var1[n_occur6$Freq == 2],]

#split again to harmonise based on exposure id
exposures6$id.exposure[exposures6$Phenotype=="exposure1"] <- 1
exposures6$id.exposure[exposures6$Phenotype=="exposure2"] <- 2
exp7 = split(exposures6, exposures6$id.exposure)[['1']]
exp8 = split(exposures6, exposures6$id.exposure)[['2']]

#harmonize exposure 1 on exposure 2
names(exp8) <- gsub("exposure", "outcome", names(exp8))
exposures6<- harmonise_data(exp7, exp8, action = 1)
#Keep only snps MrKeep= TRUE
exposures6<- exposures6[exposures6$mr_keep== TRUE, ]

#split the tables
exp7_H<- subset(exposures6, id.exposure== id.exposure[1], 
                select= c(SNP, exposure, id.exposure, effect_allele.exposure, 
                          other_allele.exposure, beta.exposure, se.exposure, pval.exposure))
#split the tables
exp8_H<- subset(exposures6, id.outcome== id.outcome[2], 
                select= c(SNP, outcome, id.outcome, effect_allele.outcome, 
                          other_allele.outcome, beta.outcome, se.outcome, pval.outcome))

#turn to exposure to merge the datasets
names(exp8_H) <- gsub("outcome", "exposure", names(exp8_H))

#merge the harmonised datasets
Exposures_H4<- rbind(exp7_H, exp8_H)
Exposures_H4["Phenotype"]<- NA
Exposures_H4$Phenotype[Exposures_H4$id.exposure == 1] <- "exposure1"
Exposures_H4$Phenotype[Exposures_H4$id.exposure == 2] <- "exposure2"
Exposures_H4["eaf.exposure"]<- NA

#specify symptom outcomes for the depression MVMR analyses
outcome_filenames <- c("scared_anx_c_14c.txt", "rsdbd_cd_c_14c.txt", 
                       "rsdbd_odd_c_14m.txt", "rsdbd_adhd_c_14m.txt")

#apply function to run MVMR on different symptom outcomes accounting for
#depression, including MVMR sensitivity analyses (MVMR Egger, Median, and Lasso)
results <- lapply(outcome_filenames, perform_mvmr_analysis)

#print or process results for each outcome
for (i in seq_along(results)) {
  print(paste("Results for:", outcome_filenames[i]))
  
  #extract the main MVMR results
  mvmr_res <- results[[i]]$mvmr_results
  #extract the F and Q statistics
  f_stat <- results[[i]]$F_statistic
  pleiotropy_res <- results[[i]]$pleiotropy_test
  #extract MVMR Egger results
  mvmr_egger <- results[[i]]$mvmr_egger
  #extract MVMR Median results
  mvmr_median <- results[[i]]$mvmr_median
  #extract MVMR Lasso results
  mvmr_lasso <- results[[i]]$mvmr_lasso
  
  #combine F and Q statistics with MVMR results + sensitivity output
  combined_results <- cbind(mvmr_res, F_stat = f_stat, Pleiotropy_stat = pleiotropy_res,
                            mvmr_median = mvmr_median, mvmr_lasso = mvmr_lasso)
  
  #round the results for presentation
  rounded_results <- round(combined_results, 2)
  
  print(rounded_results)
  print(mvmr_egger)
  
}
