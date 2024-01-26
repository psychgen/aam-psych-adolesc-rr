#07.1_IPW_weights.R

##This script creates inverse probability weights, which are added to the
##processed data for use in IPW sensitivity analyses in '07_IPW_sensitivty.R'.

#read in the data and load packages:
load(file="data/pheno_data.RData")
library(tidyverse)
library(mice)
library(weights)
library(tictoc)

#source script with modified functions
source("./scripts/07.2_utils.R")

#create function to convert dx vars from y/n to numeric
convert_to_numeric <- function(x) {
  ifelse(x == "yes", 1, ifelse(x == "no", 0, NA))
}

#select/rename phenotypic variables for IPW and create 
#individual identifier ("ind_id") for merging
ipw_pheno_data <- pheno_data %>%
  select(preg_id,BARN_NR,m_id,sex,birth_yr,parity=PARITET_5_raw,
         birthweight=VEKT_raw,gest_age=SVLEN_DG_raw,
         m_age_at_birth=MORS_ALDER_raw,p_age_at_birth=FARS_ALDER_raw,
         AA1315_raw,AA1316_raw,AA1124_raw,AA1126_raw,scl_short_m_q1,
         smfq_dep_c_14c,matches("received_dx_")) %>%
  select(-matches("received_dx_2x")) %>%
  mutate(preg_id = as.numeric(preg_id)) %>%
  mutate(ind_id = paste0(preg_id,"_",BARN_NR)) %>%
  mutate_at(vars(matches("received_dx")), convert_to_numeric) 

##manipulate variables

#compute and change response ordering of education variable (combine '2' & '3')
ipw_pheno_data <- ipw_pheno_data %>%
  mutate(education_m_q1 = case_when(AA1124_raw==1~1,
                                    AA1124_raw==2~2, 
                                    AA1124_raw==3~2,
                                    AA1124_raw==4~3,
                                    AA1124_raw==5~4,
                                    AA1124_raw==6~5,
                                    AA1124_raw==0~NA_real_)) %>%
  mutate(education_f_q1 = case_when(AA1126_raw==1~1,
                                    AA1126_raw==2~2, 
                                    AA1126_raw==3~2,
                                    AA1126_raw==4~3,
                                    AA1126_raw==5~4,
                                    AA1126_raw==6~5,
                                    AA1126_raw==0~NA_real_)) %>%
  rowwise() %>%
  mutate(parent_education_q1 = mean(c(education_m_q1,education_f_q1),na.rm=TRUE)) %>%
  ungroup() %>% 
  mutate_all(~ifelse(is.nan(.), NA, .))

#compute and recode income variable (set '0' & '8' to NA)
ipw_pheno_data <- ipw_pheno_data %>%
  mutate(income_m_q1 = case_when(AA1315_raw==1~1,
                                 AA1315_raw==2~2, 
                                 AA1315_raw==3~3,
                                 AA1315_raw==4~4,
                                 AA1315_raw==5~5,
                                 AA1315_raw==6~6,
                                 AA1315_raw==7~7,
                                 AA1315_raw==0~NA_real_)) %>%
  mutate(income_f_q1 = case_when(AA1316_raw==1~1,
                                 AA1316_raw==2~2, 
                                 AA1316_raw==3~3,
                                 AA1316_raw==4~4,
                                 AA1316_raw==5~5,
                                 AA1316_raw==6~6,
                                 AA1316_raw==7~7, 
                                 AA1316_raw==8~NA_real_,
                                 AA1316_raw==0~NA_real_)) %>%
  rowwise() %>%
  mutate(parent_income_q1 = mean(c(income_m_q1,income_f_q1),na.rm=TRUE)) %>%
  ungroup() %>% 
  mutate_all(~ifelse(is.nan(.), NA, .))

#create 14-year response variable for IPW (based on child-reported SMFQ)
ipw_pheno_data <- ipw_pheno_data %>%
  mutate(participation = ifelse(!is.na(smfq_dep_c_14c), 1, 0)) %>%
  select(-c(sex,smfq_dep_c_14c,education_m_q1,education_f_q1,income_m_q1,
            income_f_q1, matches(c("_raw"))))

#reorder variables
ipw_data <- ipw_pheno_data %>%
  mutate_at(vars(matches("received_dx")), as.factor) %>%
  mutate(participation = as.factor(participation),
         preg_id = as.character(preg_id),
         BARN_NR = as.character(BARN_NR)) %>%
  select(ind_id,preg_id,BARN_NR,m_id,birth_yr,participation,
         parity,birthweight,gest_age,m_age_at_birth,p_age_at_birth,
         parent_education_q1,parent_income_q1,scl_short_m_q1,
         matches("received_dx_"))

##first, check the extent of differences due to participation

#separate out participants and non-participants
participants = ipw_data %>% 
  filter(participation == 1)

non_participants = ipw_data %>% 
  filter(participation == 0)

#function to compare numeric variables with t-tests
compare_numeric <- function(participants, non_participants, vars) {
  numeric_results <- data.frame()
  
  for (var in vars) {
    if (is.numeric(participants[[var]])) {
      test_result <- t.test(participants[[var]], non_participants[[var]])
      temp <- data.frame(
        variable = var,
        estimate = test_result$estimate[1],
        estimate2 = ifelse(length(test_result$estimate) > 1, test_result$estimate[2], NA),
        statistic = test_result$statistic,
        p.value = test_result$p.value,
        parameter = test_result$parameter,
        conf.low = test_result$conf.int[1],
        conf.high = test_result$conf.int[2],
        method = "t-test",
        type = "numeric"
      )
      numeric_results <- rbind(numeric_results, temp)
    }
  }
  
  return(numeric_results)
}

#function to compare binary variables with fisher.test
compare_binary <- function(participants, non_participants, vars) {
  binary_results <- data.frame()
  
  for (var in vars) {
    if (is.factor(participants[[var]]) && all(levels(participants[[var]]) %in% c("0", "1"))) {
      yes_part <- sum(participants[[var]] == "1", na.rm = TRUE)
      no_part <- sum(participants[[var]] == "0", na.rm = TRUE)
      yes_non_part <- sum(non_participants[[var]] == "1", na.rm = TRUE)
      no_non_part <- sum(non_participants[[var]] == "0", na.rm = TRUE)
      contingency_table <- matrix(c(yes_part, no_part, yes_non_part, no_non_part), nrow = 2)
      test_result <- fisher.test(contingency_table)
      temp <- data.frame(
        variable = var,
        p.value = test_result$p.value,
        OR = test_result$estimate,
        method = "Fisher's Exact Test",
        type = "binary"
      )
      binary_results <- rbind(binary_results, temp)
    }
  }
  
  return(binary_results)
}

# Define the range of columns to apply the functions to
numeric_vars <- names(participants)[which(names(participants) == "parity"):which(names(participants) == "scl_short_m_q1")]
binary_vars <- names(participants)[which(names(participants) == "received_dx_dep_npr"):which(names(participants) == "received_dx_beh_kuhr")]

# Run the functions
numeric_comparison <- compare_numeric(participants, non_participants, numeric_vars) %>%
  mutate(diff=estimate-estimate2)
binary_comparison <- compare_binary(participants, non_participants, binary_vars)

# Print the comparison results
print(numeric_comparison)
print(binary_comparison)

ggplot(binary_comparison %>% 
         filter(variable%in%names(ipw_data %>% select(received_dx_dep_npr:received_dx_beh_kuhr)))
       , aes(x=OR, y=variable))+
  geom_vline(aes(xintercept=1),colour="grey60", linetype=2,size=1.1)+
  geom_point(shape=21, size=4, stroke=1.1, fill="blue")+
  theme_minimal()+
  scale_x_continuous("Difference participators and non-participators (OR)")+
  scale_y_discrete("Predictor")

###############################################################################  
################### 2 - MI and IPW to set up adjustment #######################
###############################################################################

### Multiple imputation

MI <- mice(ipw_data, maxit=15, m=50, seed=30)

### IPW
all = complete(MI, "all")
all_weights = ipw_data %>% 
  select(ind_id)

for (i in 1:length(all)){
  
  comp_dat <-  all[[i]] 
  
  # use to generate ipw (stabilized)
  ipws <- myipwpoint(exposure=participation,
                     family="binomial",
                     link="logit",
                     numerator = ~ 1,
                     denominator = ~ 1 + birthweight + gest_age + m_age_at_birth + p_age_at_birth + parent_education_q1 + parent_income_q1 + scl_short_m_q1 + received_dx_dep_kuhr + received_dx_dep_npr + received_dx_anx_kuhr + received_dx_anx_npr + received_dx_beh_kuhr + received_dx_beh_npr + received_dx_adhd_kuhr + received_dx_adhd_npr,
                     data=comp_dat)
  
  range(ipws$ipw.weights)
  
  # check accuracy of predictions
  test <- comp_dat %>% 
    mutate(pred = round(predict(ipws$den.mod, newdata=., type="response"))) 
  
  accuracy <- table(test$pred, test$participation)
  print(sum(diag(accuracy))/sum(accuracy))
  
  all_weights =cbind(all_weights, ipws$ipw.weights)
}

# rename and save raw weights
colnames(all_weights)<-c("ind_id", paste0("ipw",seq(1,50)))

save(all_weights, file="./data/all_weights.RData")
load(file="./data/all_weights.RData")

##############################################################################
############### Test out adding weights to each imputed dataset ##############
##############################################################################

full <- complete(MI, action = "long", include = TRUE)

long_weights= ipw_data %>% 
  select(ind_id) %>% 
  mutate(ipw=NA) %>% 
  bind_rows(all_weights %>% 
            pivot_longer(-ind_id, values_to = "ipw") %>% 
            arrange(name,ind_id) %>% 
            select(-name)) %>% 
  mutate(`.imp`=rep(seq(0,50,1), each=nrow(ipw_data)))

full_weights = full %>%
  left_join(long_weights)

midat = as.mids(full_weights)

##################################################################
####### Smooth weights by averaging across all imputations #######
##################################################################

smoothed_weights <- all_weights %>% 
  rowwise() %>% 
  group_by(ind_id) %>% 
  summarise(ipw= mean(c_across(ipw1:ipw50)))

save(smoothed_weights, file="./data/smoothed_weights.RData")
load(file="./data/smoothed_weights.RData")

# check performance with weighted t-tests
wdat <- ipw_data %>% 
  left_join(smoothed_weights) %>%
  mutate(ipw = as.character(ipw),
         birth_yr = as.character(birth_yr)) %>%
  mutate_if(is.numeric, scale) %>%
  mutate_if(is.numeric, as.numeric) %>%
  mutate(ipw = as.numeric(ipw),
         birth_yr = as.numeric(birth_yr))

part_wdat = wdat %>% 
  filter(participation==1)

nonpart_wdat = wdat %>% 
  filter(participation==0)

# make function to compare numeric vars using weighted t-tests
wtcompare = function(part, nonpart, vars, wts){
  all=data.frame()
  for(var in vars){
    temp=weights::wtd.t.test(x=part[[var]], y= nonpart[[var]], weight=part[[wts]]) 
    temp = c(temp$additional,temp$coefficients) %>% 
      as_tibble_row() %>% 
      rename_all(~paste0("wtd.",.))
    base= weights::wtd.t.test(x=part[[var]], y= nonpart[[var]]) 
    base = c(base$additional,base$coefficients) %>% 
      as_tibble_row() %>% 
      rename_all(~paste0("unw.",.))
    comp= cbind(base,temp) %>% 
      mutate(pred=var) %>% 
      select(pred,everything())
    all=rbind(all,comp)
  }
  return(as_tibble(all))
}

# run function
wtcomparison = wtcompare(part_wdat,nonpart_wdat,names(ipw_data %>%
                                                        select(parity:scl_short_m_q1)), "ipw")
wtcomparison %>% filter(pred%in%names(ipw_data %>% select(parity:scl_short_m_q1))) %>%  print(n=Inf)

wtcomplong = wtcomparison %>% 
  filter(pred%in%names(ipw_data %>% select(parity:scl_short_m_q1))) %>% 
  select(pred, unw.Difference,`unw.Std. Err`,wtd.Difference,`wtd.Std. Err`) %>% 
  pivot_longer(cols=-pred) %>% 
  mutate(name = str_replace_all(name,"Std\\. ","Std")) %>% 
  separate(name, into=c("model","param"), sep="\\." ) %>% 
  pivot_wider(names_from = param, values_from = value) %>% 
  mutate(conf.low=Difference-1.96*StdErr,
         conf.high=Difference+1.96*StdErr,
         pred=as.factor(pred),
         model=as.factor(model))


# create plot
p<-ggplot(wtcomplong, aes(x=Difference, y=pred, fill=factor(model)))+
  geom_vline(aes(xintercept=0),colour="black", linetype=2,size=1)+
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), height=0, size=1.1, position=position_dodge(0.4))+
  geom_point(shape=21, size=4, stroke=1.1, position=position_dodge(0.4))+
  theme_minimal(base_size=16)+
  scale_x_continuous("Mean difference between participants and non-participants (SD)")+
  scale_fill_manual(values=c("blue","red"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=18),
        legend.position = "top",
        axis.title.y = element_blank(),
        axis.text = element_text(size=16))
p

wtcompare = function(part, nonpart, vars, wts){
  all=data.frame()
  for(var in vars){
    temp=weights::wtd.t.test(x=part[[var]], y= nonpart[[var]], weight=part[[wts]]) 
    temp = c(temp$additional,temp$coefficients) %>% 
      as_tibble_row() %>% 
      rename_all(~paste0("wtd.",.))
    base= weights::wtd.t.test(x=part[[var]], y= nonpart[[var]]) 
    base = c(base$additional,base$coefficients) %>% 
      as_tibble_row() %>% 
      rename_all(~paste0("unw.",.))
    comp= cbind(base,temp) %>% 
      mutate(pred=var) %>% 
      select(pred,everything())
    all=rbind(all,comp)
  }
  return(as_tibble(all))
}

# run function
wtcomparison = wtcompare(part_wdat,nonpart_wdat,names(ipw_data %>%
                                                        select(parity:scl_short_m_q1)), "ipw")
wtcomparison %>% filter(pred%in%names(ipw_data %>% select(parity:scl_short_m_q1))) %>%  print(n=Inf)

#save out processed data with IP weights
save(part_wdat, file="data/ipw_wdat.RData")
