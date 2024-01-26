#00.1_process_dx_data.R

##This script takes the phenotools-generated dataset from '00_prepare_data.R',
##reshapes it to a long format and creates dx variables for each year.
##It then accounts for the censoring in the data by setting unobserved
##years to missing, saving out the processed dx data.

#load required packages
library(tidyverse)
library(phenotools)

#load phenotools-generated dataset
load(file="data/pheno_data.RData")

#get long-form registry data
my_long_npr_data <- pivot_curated_npr(pheno_data,large_file=TRUE)
my_long_kuhr_data <- pivot_curated_kuhr(pheno_data,large_file=TRUE)

#process long npr data
long_npr <- my_long_npr_data %>% 
  select(preg_id, BARN_NR, birth_yr, dx_group, hc_contact, info, value) %>%
  mutate(childageyrs = case_when(info == "childageyrs" ~ value)) %>%
  filter(info=="childageyrs") %>%
  mutate(childageyrs = as.numeric(childageyrs))

hist(long_npr$childageyrs)
table(long_npr$childageyrs)

#loop through each year and each diagnostic group,
#creating a variable with 1 if a person got an NPR 
#diagnosis in that group that year, or 0 if not
years <- rep(1:17, each = 1)
diags <- unique(long_npr$dx_group)

for (year in years){
  for (diag in diags){
    
    name <- paste0(diag,"_",year)
    
    print(name)
    
    long_npr[[name]] <- ifelse(
      long_npr$dx_group == diag & is.na(long_npr$childageyrs), 0,
      ifelse(long_npr$dx_group == diag & long_npr$childageyrs == year, 
             1, 
             0)
    )
  }
}

npr <- long_npr %>%
  select(preg_id, BARN_NR, birth_yr, matches("npr")) %>%
  mutate(ind_id = paste0(preg_id,"_",BARN_NR)) %>%
  distinct() %>%
  group_by(ind_id,birth_yr) %>%
  summarise(across(matches("npr"), sum, na.rm = TRUE))

#repeat the same steps for kuhr data
long_kuhr <- my_long_kuhr_data %>% 
  select(preg_id, BARN_NR, birth_yr, dx_group, hc_contact, info, value) %>%
  mutate(childageyrs = case_when(info == "childageyrs" ~ value)) %>%
  mutate(childageyrs = as.numeric(childageyrs))

hist(long_kuhr$childageyrs)
table(long_kuhr$childageyrs)

#loop through each year and each diagnostic group,
#creating a variable with 1 if a person got a KUHR 
#diagnosis in that group that year, or 0 if not
years <- rep(1:17, each = 1)
diags <- unique(long_kuhr$dx_group)

for (year in years){
  for (diag in diags){
    
    name <- paste0(diag,"_",year)
    
    print(name)
    
    long_kuhr[[name]] <- ifelse(
      long_kuhr$dx_group == diag & is.na(long_kuhr$childageyrs), 0,
      ifelse(long_kuhr$dx_group == diag & long_kuhr$childageyrs == year, 
             1, 
             0)
    )
  }
}

kuhr <- long_kuhr %>%
  select(preg_id, BARN_NR, birth_yr, matches("kuhr")) %>%
  mutate(ind_id = paste0(preg_id,"_",BARN_NR)) %>%
  distinct() %>%
  group_by(ind_id, birth_yr) %>%
  summarise(across(matches("kuhr"), sum, na.rm = TRUE))

#filter kuhr and npr to ages 0-17, the 
#ages that are eligible for inclusion
long_kuhr <- long_kuhr %>%
  filter(childageyrs <= 17) %>%
  filter(childageyrs > 0)

long_npr <- long_npr %>%
  filter(childageyrs <= 17) %>%
  filter(childageyrs > 0)

#plot birth year against child age when receiving diagnoses, to see what 
#combinations of birth year and age of diagnosis are censored in the registries
plot(long_npr$birth_yr,long_npr$childageyrs)
plot(long_kuhr$birth_yr,long_kuhr$childageyrs)

#merge npr and kuhr data
dx_data <- npr %>%
  left_join(kuhr)

#apply ifelse() function to each pair of npr and kuhr columns from 1-17,
#to create new columns including 1 if either npr or kuhr is 1, otherwise 0
for (domain in c("adhd","beh","dep","anx")) {
  for (i in 1:17) {
    dx_data[paste0(domain,"_dx_",i)] <- ifelse(
      dx_data[paste0(domain,"_npr_",i)] | dx_data[paste0(domain,"_kuhr_",i)] > 0,
      1,
      0
    )
  }
  # drop the original columns
  dx_data <- dx_data %>% 
    select(-matches(paste0("^", domain, "_npr_.*|^", domain, "_kuhr_.*")))
}

#set censored observations to missing based on birth year and time of diagnosis:
#first, apply criteria on columns corresponding to age 13-17 with a nested loop,
#where values are censored if birth_yr is higher than 2008 for age 13, or higher 
#than 2007 for age 14 and so on until if birth_yr is higher than 2004 for age 17 

domains <- c("adhd","beh","dep","anx")
for (i in 13:17) {
  for (j in domains) {
    
    name <- paste0(j,"_dx_",i)
    print(name)
    
    if (grepl(i, name)) {
      dx_data[[name]][dx_data$birth_yr > 2008-(i-13)] <- NA
      
    }
  }
}

#then apply criteria across columns corresponding to age 1-9:
#here, values are censored if birth_yr is lower than 2007 for age 1, or 2006 for 
#age 2, 2005 for age 3, and so on until if birth_yr is lower than 2000 for age 8

for (i in 1:9) {
  for (j in domains) {
    
    name <- paste0(j,"_dx_",i)
    print(name)
    
    if (grepl(i, name)) {
      dx_data[[name]][dx_data$birth_yr < 2008-i] <- NA
      
    }
  }
}

#save out processed child dx data
save(dx_data, file="./data/dx_data.RData")
