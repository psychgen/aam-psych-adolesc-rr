#XX_project_diagnoses.R

library(tidyverse)
library(phenotools)
set.seed(23904)

## Curate datatset with variables from MoBa phenotypic data
load(file="data/pheno_data.RData")
diagnoses = c("dep","anx","adhd","beh")
all_case_rates <- data.frame()
for( diag in diagnoses){
  
  npr <- pivot_curated_npr(pheno_data %>% select(preg_id:dx_owner, matches(diag)))
  
  kuhr <- pivot_curated_kuhr(pheno_data %>% select(preg_id:dx_owner, matches(diag)))
  
  #First get case rate per age for all complete
  
  #denom is all girls with full data age 10-17
  denom <-pheno_data %>% 
    select(preg_id:birth_yr) %>% 
    distinct() %>% 
    mutate(years_on_rec = 2018-birth_yr) %>% 
    filter(years_on_rec >=17)
  
  
  #rates per year 
  case_rates <- npr %>% 
    bind_rows(kuhr) %>% 
    select(preg_id:birth_yr, info, value) %>% 
    filter(info=="childageyrs") %>% #Just focus on how old when they got the dx
    mutate(years_on_rec = 2018-birth_yr) %>% 
    filter(years_on_rec >=17, #Only take those in the denominator set
           value  <=17) %>%  #Only take those years where all in the denominator can contribute
    mutate(iid= paste0(preg_id, BARN_NR)) %>% 
    group_by(value,iid) %>%
    slice_sample(n=1) %>% #Important so that only one hc_contact per individual per year
    ungroup() %>% 
    group_by(value) %>% 
    summarise(rate= n()/nrow(denom)) %>% 
    mutate(diagnosis = diag)
  
  all_case_rates <- bind_rows(all_case_rates,case_rates)
}


# Now we sample from our 8year dataset to get possible versions of our analytic
# dataset, simulate an AaM variable, simulate diagnoses per individual per year,
# based on the rates above, and calculate a yes/no diagnosis outcome for each
# category based on projected data availability as of 2022 and time window

sample_size=c(12500,10000)
timewindows=c(100)
niterations=50
bounds=c(FALSE)
impute_post_2021=c(TRUE)
all_exp_rates <- data.frame()
for(impute in impute_post_2021){
  for(lbound in bounds){
    for (timewindow in timewindows){
      for(iter in 1:niterations){
        for(sample_size in sample_size){
        print(iter)
        dat <- pheno_data %>%
          select(preg_id:smfq_dep_c_8yr) %>% 
          filter(!is.na(smfq_dep_c_8yr)) %>%
          slice_sample(n= sample_size) %>% 
          mutate(AaM = round (rnorm(n = sample_size, 
                                    mean = 151.52,
                                    sd = 14.11)/12),
                 null=1,
                 window=9+timewindow) 
        
        
        if(impute==FALSE){
          dat <- dat %>% filter((birth_yr+window)<2022 )  #restrict to only those who we will have data for up to end 2021
        }
        
        
        ess<- nrow(dat) #this is our effective sample size for the analysis
        
        dat <- dat %>%
          left_join(all_case_rates %>% 
                      unite(diag_age, c("diagnosis","value") ) %>% 
                      mutate(null=1)) 
        
        dat$dx <-  rbinom(n=nrow(dat),size=1,prob=dat$rate )
        
        exp_rates<-dat %>% 
          separate(diag_age, into = c("diag","age")) %>% 
          select(preg_id, BARN_NR,diag,window,diag,AaM,age,dx) %>% 
          mutate(age=as.numeric(age)) %>% 
          group_by(preg_id,BARN_NR,diag) %>% 
          filter(age<=window)
        
        if(lbound==TRUE){
          exp_rates <- exp_rates %>% filter(age>=AaM) 
        }
        exp_rates <- exp_rates %>% 
          summarise(n_dx=sum(dx)) %>% 
          ungroup() %>% 
          mutate(dx=ifelse(n_dx>0,1,0)) %>% 
          group_by(diag) %>% 
          summarise(rate=sum(dx)/sample_size) %>% 
          mutate(repl=iter,
                 tw=timewindow,
                 ess=ess,
                 lbound_at_aam=lbound,
                 impute_post_2021=imput,
                 sample_size=sample_size)
        
        all_exp_rates <- bind_rows(all_exp_rates,exp_rates)
        
      }
    }
  }
  }
}

save(all_exp_rates, file= "./output/projected_diagnosis_rates.RData")

all_exp_rates %>% group_by(diag,tw,lbound_at_aam, impute_post_2021) %>% 
  summarise(se=sd(rate)/sqrt(100),
            rate=mean(rate)) %>% 
  mutate(lci=rate+1.96*se,
         uci=rate-1.96*se)
