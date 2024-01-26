#00.3_multiple_imputation.R

##This script conducts multiple imputation of all variables (including 
##covariates to avoid dropping a large number of observations), using the 
##Multivariate Imputation by Chained Equations (mice) package. It also conducts 
##a sensitivity analysis related to the imputation of age at menarche. 

#load data and required packages
load(file="data/analysis_dataset.RData")
library(dplyr)
library(mice)
           
#reorder and transform variables to ensure they are handled correctly by mice
fulldata <- fulldata %>%
  select(ind_id,preg_id,BARN_NR,m_id,parity,m_age_at_birth,
         p_age_at_birth,parent_education_q1,parent_income_q1,
         scl_short_m_q1,scl_full_m_q3,epds_short_m_6m,menarche_c_14c,
         bmi_derived_c_14c,breast_c_14c,skin_c_14c,hair_c_14c,
         growth_c_14c,aam_c_14c,age_8yr,age_14yr,financ_probs_18m,
         parent_cohab_18m,parent_cohab_3yr,smfq_dep_c_14c,
         rsdbd_cd_c_14c,scared_anx_c_14c,scared_anx_c_8yr,
         smfq_dep_c_8yr,rsdbd_cd_c_8yr,rsdbd_odd_c_8yr,
         rsdbd_adhd_c_8yr,rsdbd_odd_c_14m,rsdbd_adhd_c_14m,
         bmi_derived_c_8yr,PGS_AaM,anx_dx_ch,anx_dx_ad,adhd_dx_ch,
         adhd_dx_ad,dep_dx_ad,beh_dx_ch,beh_dx_ad) %>%
  mutate(parent_cohab_18m = as.factor(parent_cohab_18m),
         parent_cohab_3yr = as.factor(parent_cohab_3yr),
         aam_c_14c = as.numeric(aam_c_14c),
         menarche_c_14c = as.numeric(menarche_c_14c),
         breast_c_14c = as.numeric(breast_c_14c),
         growth_c_14c = as.numeric(growth_c_14c),
         hair_c_14c = as.numeric(hair_c_14c),
         skin_c_14c = as.numeric(skin_c_14c)) %>%
  mutate(across(matches("dx"), as.factor))

#check number of missing
pMiss <- function(x) {
  sum(is.na(x)) / length(x) * 100
}
apply(fulldata, 2, pMiss)

#create the condition vector outside of mice
menarche_condition <- !is.na(fulldata$menarche_c_14c) & fulldata$menarche_c_14c == 0 & is.na(fulldata$aam_c_14c)
menarche_condition_subset <- menarche_condition[which(is.na(fulldata$aam_c_14c))]

#settings for multiple imputation
init = mice(fulldata, maxit = 1, m = 1)
predictorMatrix = init$predictorMatrix
meth = init$method
post = init$post

#remove ID vars as they provide little info, and menarche which is collinear
predictorMatrix[, c("preg_id", "BARN_NR", "m_id", "ind_id","menarche_c_14c")]=0

#specify which method to use for imputation
meth[c("menarche_c_14c", "aam_c_14c")] = c("", "norm")

post["aam_c_14c"] <- "imp[[j]][, i] <- ifelse(!is.na(menarche_condition_subset) & menarche_condition_subset, squeeze(imp[[j]][, i], c(9, 12)), imp[[j]][, i])"

#run imputation
imp <-
  mice(
    fulldata,
    method = meth,
    maxit = 15,
    m = 50,
    predictorMatrix = predictorMatrix,
    post = post,
    seed = 37
  )

#check logged events
imp$loggedEvents

#get summary
summary(imp)

#get all multiply imputed datasets
imputed_datasets <- complete(imp, action = 'long', include = TRUE)

#save out/load imputed datasets
save(imp, file = "./data/imputations.RData")
save(imputed_datasets, file = "./data/imputed_datasets.RData")

#to validate the imputation of age at menarche we set an equivalent %
#of randomly selected individuals to missing, imputing their aam, and
#comparing the imputed and actual values (sensitivity analysis)

#load imputed datasets
load(file = "./data/imputations.RData")

#take aam column from one imputed dataset
aam_mi <- complete(imp, 5)$aam_c_14c

#append the imputed aam variable to original dataset
fulldata_mi <- fulldata %>% cbind(aam_mi)

#set seed to ensure reproducibility
set.seed(10)

#identify rows where aam (non-imputed) is non-missing
selected_rows <- which(!is.na(fulldata_mi$aam_c_14c))

#randomly sample 972 rows from this selection
na_indices <- sample(selected_rows, 972)

#set aam_mi to NA in these rows
fulldata_mi$aam_mi[na_indices] <- NA

#remove original aam
fulldata_mi <- fulldata_mi %>% select(-aam_c_14c)

#settings for multiple imputation
init_sens = mice(fulldata_mi, maxit = 1, m = 1)
predictorMatrix_sens = init_sens$predictorMatrix
meth_sens = init_sens$method

#remove ID vars as predictors as they provide little info
predictorMatrix_sens[, c("preg_id", "BARN_NR", "m_id", "ind_id", "menarche_c_14c")] =0

#specify which method to use for imputation
meth_sens[c("menarche_c_14c", "aam_mi")] = c("", "norm")

#run imputation on sensitivity dataset
imp_sens <-
  mice(
    fulldata_mi,
    method = meth_sens,
    maxit = 15,
    m = 50,
    predictorMatrix = predictorMatrix_sens,
    seed = 3
  )

#save out/load imputed sensitivity datasets
save(imp_sens, file = "./data/imputations_sensitivity.RData")

#number of imputed datasets
m <- 50

#extract the completed datasets
imputed_datasets_sens <-
  lapply(1:m, function(i)
    complete(imp_sens, action = i))

#calculate the mean imputed values across all datasets
mean_imputed_values <-
  Reduce("+", lapply(imputed_datasets_sens, `[[`, "aam_mi")) / m

#now 'mean_imputed_values' has the average imputed values for the variable 'aam_mi'

#compare imputed and actual values for cases set to missing
comparison <- data.frame(actual = fulldata$aam_c_14c[na_indices],
                         imputed = mean_imputed_values[na_indices])

#check means
(mean(comparison$actual))
(mean(comparison$imputed))

#compute Mean Absolute Error (MAE)
(mae <- mean(abs(comparison$actual - comparison$imputed)))

#assess the relationship between the actual and imputed values
(cor <- cor(comparison$actual, comparison$imputed, use = "complete.obs"))
