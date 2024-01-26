#00.4_process_imputed_data.R

##This script does processing of the multiply imputed datasets from 
##'00.3_multiple_imputation.R', including transformation and 
##standardisation of variables for main analyses and categorisation
##of age at menarche/dichotomisation of depressive symptoms outcome
##for sensitivity analyses

#load data and required packages
load(file = "./data/imputed_datasets.RData")
library(dplyr)
library(mice)
library(moments)
library(MASS)
library(purrr)

#split imputed datasets into a list to enable processing
imputed_dataset_list <- split(imputed_datasets, imputed_datasets$.imp)

#function to categorise age at menarche and smfq
categorise_vars <- function(df, aam, smfq_14yr, smfq_8yr) {
  mean_aam <- mean(df[[aam]], na.rm = TRUE)
  sd_aam <- sd(df[[aam]], na.rm = TRUE)
  
  #categorise age at menarche
  df$aam_cat <- ifelse(df[[aam]] <= mean_aam - sd_aam, 'early',
                       ifelse(df[[aam]] >= mean_aam + sd_aam, 'late', 'average'))
  
  #make factor
  df$aam_cat <- as.factor(df$aam_cat)
  
  #dichotomise smfq
  df$smfq_14yr_dic <- as.factor(ifelse(df[[smfq_14yr]] > 16, 1, 0))
  df$smfq_8yr_dic <- as.factor(ifelse(df[[smfq_8yr]] > 11, 1, 0))
  
  return(df)
}

#apply the function to each dataset in the list
categorised_data <-
  lapply(
    imputed_dataset_list,
    categorise_vars,
    aam = "aam_c_14c",
    smfq_14yr = "smfq_dep_c_14c",
    smfq_8yr = "smfq_dep_c_8yr"
  )

#calculate the proportion of cases in smfq_14yr
proportion_1s <- function(df, smfq) {
  prop_1 <- mean(df[[smfq]] == 1, na.rm = TRUE)
  return(prop_1)
}

#apply the function to each imputed dataset and store the proportions
proportions <- sapply(categorised_data, proportion_1s, smfq = "smfq_14yr_dic")

#compute the mean proportion across all imputed datasets
mean_proportion <- mean(proportions)

#print mean proportion
print(mean_proportion)#0.160 vs 0.155 in ALSPAC

#convert the list of data frames back to mids object
MI_cat <- mice::as.mids(do.call(rbind, categorised_data))

#convert from mids object to list of datasets
datlist_cat <- miceadds::mids2datlist(MI_cat)

#save out
save(MI_cat, file = "./data/MI_categorised.RData")
save(datlist_cat, file = "./data/datlist_categorised.RData")

##Calculate descriptive statistics of age at menarche and depression  

#initialize a matrix to store mean and sd values
num_datasets <- length(imputed_dataset_list) - 1
mean_matrix <- matrix(nrow = num_datasets, ncol = 2)
sd_matrix <- matrix(nrow = num_datasets, ncol = 2)

for (i in 1:num_datasets) {
  dataset <- imputed_dataset_list[[i]]
  mean_values <- sapply(dataset[c("smfq_dep_c_14c", "aam_c_14c")], mean)
  sd_values <- sapply(dataset[c("smfq_dep_c_14c", "aam_c_14c")], sd)
  mean_matrix[i,] <- mean_values
  sd_matrix[i,] <- sd_values
}

#calculate the mean and sd values for each variable
mean_vals <- colMeans(mean_matrix, na.rm = TRUE)
sd_vals <- colMeans(sd_matrix, na.rm = TRUE)

#print the mean and sd values
print(mean_vals)
print(sd_vals)

##check skewness and kurtosis of outcome vars before transformation
skewness_values_8y <- sapply(imputed_dataset_list$`1`[c("smfq_dep_c_8yr",
                                                        "scared_anx_c_8yr",
                                                        "rsdbd_cd_c_8yr",
                                                        "rsdbd_odd_c_8yr",
                                                        "rsdbd_adhd_c_8yr")], skewness)
kurtosis_values_8y <- sapply(imputed_dataset_list$`1`[c("smfq_dep_c_8yr",
                                                        "scared_anx_c_8yr",
                                                        "rsdbd_cd_c_8yr",
                                                        "rsdbd_odd_c_8yr",
                                                        "rsdbd_adhd_c_8yr")], kurtosis)
skewness_values_14y <- sapply(imputed_dataset_list$`1`[c("smfq_dep_c_14c",
                                                        "scared_anx_c_14c",
                                                        "rsdbd_cd_c_14c",
                                                        "rsdbd_odd_c_14m",
                                                        "rsdbd_adhd_c_14m")], skewness)
kurtosis_values_14y <- sapply(imputed_dataset_list$`1`[c("smfq_dep_c_14c",
                                                        "scared_anx_c_14c",
                                                        "rsdbd_cd_c_14c",
                                                        "rsdbd_odd_c_14m",
                                                        "rsdbd_adhd_c_14m")], kurtosis)
#print results
print(skewness_values_8y)
print(kurtosis_values_8y)
print(skewness_values_14y)
print(kurtosis_values_14y)

#Distributional properties vary a lot between outcomes - therefore,
#apply a box cox transformation to decide the appropriate transformation
dep14<-boxcox(lm(imputed_dataset_list$`1`$smfq_dep_c_14c + 1 ~ 1)) 
anx14<-boxcox(lm(imputed_dataset_list$`1`$scared_anx_c_14c + 1 ~ 1)) 
cd14<-boxcox(lm(imputed_dataset_list$`1`$rsdbd_cd_c_14c + 1 ~ 1)) 
adhd14<-boxcox(lm(imputed_dataset_list$`1`$rsdbd_adhd_c_14m + 1 ~ 1)) 
odd14<-boxcox(lm(imputed_dataset_list$`1`$rsdbd_odd_c_14m + 1 ~ 1)) 

dep8<-boxcox(lm(imputed_dataset_list$`1`$smfq_dep_c_8yr + 1 ~ 1)) 
anx8<-boxcox(lm(imputed_dataset_list$`1`$scared_anx_c_8yr + 1 ~ 1)) 
cd8<-boxcox(lm(imputed_dataset_list$`1`$rsdbd_cd_c_8yr + 1 ~ 1)) 
adhd8<-boxcox(lm(imputed_dataset_list$`1`$rsdbd_adhd_c_8yr + 1 ~ 1)) 
odd8<-boxcox(lm(imputed_dataset_list$`1`$rsdbd_odd_c_8yr + 1 ~ 1)) 

#calculate exact lambdas
(dep14_lambda <- dep14$x[which.max(dep14$y)])     #lambda = 0.34
(anx14_lambda <- anx14$x[which.max(anx14$y)])     #lambda = 0.14
(cd14_lambda <- cd14$x[which.max(cd14$y)])        #lambda = -2
(adhd14_lambda <- adhd14$x[which.max(adhd14$y)])  #lambda = -0.10
(odd14_lambda <- odd14$x[which.max(odd14$y)])     #lambda = -0.14

(dep8_lambda <- dep8$x[which.max(dep8$y)])     #lambda = -0.42
(anx8_lambda <- anx8$x[which.max(anx8$y)])     #lambda = -0.30
(cd8_lambda <- cd8$x[which.max(cd8$y)])        #lambda = -2
(adhd8_lambda <- adhd8$x[which.max(adhd8$y)])  #lambda = 0.22
(odd8_lambda <- odd8$x[which.max(odd8$y)])     #lambda = 0.14

#create list of symptom domains
domains <- list("smfq_dep_c_8yr","smfq_dep_c_14c",
                "scared_anx_c_8yr","scared_anx_c_14c",
                "rsdbd_adhd_c_8yr","rsdbd_adhd_c_14m",
                "rsdbd_cd_c_8yr","rsdbd_cd_c_14c",
                "rsdbd_odd_c_8yr","rsdbd_odd_c_14m")

transf <- list("sqrt", "sqrt",
               "sqrt", "log",
               "log", "log",
               "log", "log", 
               "log", "log")

#create a named list of transformation functions
transf_functions <- list(
  "sqrt" = function(x) sqrt(x),
  "log" = function(x) log1p(x)
)

# Apply transformations
transformed_data <- lapply(split(imputed_datasets, imputed_datasets$.imp), function(df) {
  for(i in seq_along(domains)) {
    # Get the transformation function name for the current domain
    trans_fun_name <- transf[[i]]
    # Get the actual function from the transf_list
    trans_fun <- transf_functions[[trans_fun_name]]
    # Apply the transformation
    df[[domains[[i]]]] <- trans_fun(df[[domains[[i]]]])
  }
  df
})

##check skewness and kurtosis of outcome vars after transformation
skewness_values_8y <- sapply(transformed_data$`1`[c("smfq_dep_c_8yr",
                                                    "scared_anx_c_8yr",
                                                    "rsdbd_cd_c_8yr",
                                                    "rsdbd_odd_c_8yr",
                                                    "rsdbd_adhd_c_8yr")], skewness)
kurtosis_values_8y <- sapply(transformed_data$`1`[c("smfq_dep_c_8yr",
                                                    "scared_anx_c_8yr",
                                                    "rsdbd_cd_c_8yr",
                                                    "rsdbd_odd_c_8yr",
                                                    "rsdbd_adhd_c_8yr")], kurtosis)
skewness_values_14y <- sapply(transformed_data$`1`[c("smfq_dep_c_14c",
                                                     "scared_anx_c_14c",
                                                     "rsdbd_cd_c_14c",
                                                     "rsdbd_odd_c_14m",
                                                     "rsdbd_adhd_c_14m")], skewness)
kurtosis_values_14y <- sapply(transformed_data$`1`[c("smfq_dep_c_14c",
                                                     "scared_anx_c_14c",
                                                     "rsdbd_cd_c_14c",
                                                     "rsdbd_odd_c_14m",
                                                     "rsdbd_adhd_c_14m")], kurtosis)
#print results
print(skewness_values_8y)
print(kurtosis_values_8y)
print(skewness_values_14y)
print(kurtosis_values_14y)

#standardise the variables for each imputed dataset
standardised_data <- lapply(transformed_data, function(df){
  df %>%
    mutate(across(c("aam_c_14c", "PGS_AaM", unlist(domains)), ~ as.vector(scale(.))))
})

#convert the list of data frames back to mids object
standardised_MI <- as.mids(do.call(rbind, standardised_data))

#convert from mids object to list of datasets
datlist <- miceadds::mids2datlist(standardised_MI)

#save out
save(standardised_MI, file = "./data/standardised_MI_tr.RData")
save(datlist, file = "./data/datlist_tr.RData")

