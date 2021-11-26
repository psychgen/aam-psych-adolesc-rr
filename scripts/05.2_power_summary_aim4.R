# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Summary of power analyses for aim 4 hypothesis/equivalence tests (`05.1_power_summary_aim4.R`) 
#'
# ---- echo=F, out.width='33%', out.extra = 'style="float:right;"'
knitr::include_graphics("052_power_summary_aim4_codemap.png", error = FALSE)
#'
#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to summarise the power analyses corresponding to 
#'Aim 4 hypothesis/equivalence tests. 
#'
#'  
#'
#'&nbsp;
#'&nbsp;
#'
#'## Load required packages
# ---- 

library(tidyverse)
library(effectsize)
library(parameters)

load(file="./output/fullsim_power34a.RData")
#'### Power 
#'Our power in NHST in each scenario is the proportion of iterations in which a true 
#'effect is detected, given the effect exists in the population
#'
#'Our power in equivalence/inferiority each scenario is the proportion of iterations 
#'in which an effect is declared equivalent to zero, given that in the population it is 
#'less than the SESOI
#'
# ----


#Summarise at the level of scenario (i.e., across iterations)
power4a <- fullsim %>% 
  group_by(causal_eff_d,R2_inst, sample_size) %>% 
  summarise(power_nhst= mean(p.value < 0.05),
            se_nhst=sd(p.value < 0.05)/sqrt(n()),
            power_equiv= mean(conf.low  > ROPE_low_aim4 & conf.high < ROPE_high_aim4 ),
            se_equiv=sd(conf.low  > ROPE_low_aim4 & conf.high < ROPE_high_aim4)/sqrt(n())) %>% 
  rename(effect_size= causal_eff_d,N = sample_size) %>% 
  pivot_longer(cols=matches("power|se_"), names_to = "type", values_to="power") %>% 
  separate(type, into=c("est","test")) %>% 
  pivot_wider(names_from=est,values_from=power) %>% 
  mutate( test=factor(test,
                      levels=c("nhst","equiv"),
                      labels=c("NHST\n(H1:Effect > 0)","Equivalence\n(H1:Effect = 0)")),
          N=factor(N))%>% 
  filter(N%in%c("9500","10500"))


r2labs <- c("R2 inst. 5%","R2 inst. 7.5%","R2 inst. 10%")
names(r2labs)<- c("0.05","0.075", "0.1")
#Plot the power curves 
ggplot(data=power4a, aes(x= effect_size, y= power, fill=N,
                         colour=N,group=N))+
  geom_errorbar(aes(ymin=power-se,ymax=power+se), width=0, size=1.2, alpha=0.4)+
  geom_point(size=2.6,
             alpha=0.4)+
  geom_line(size=0.8)+
  geom_hline(yintercept=0.95, colour="red",size=1.4, linetype=2)+
  geom_hline(yintercept=0.80, colour="orange",size=1.4, linetype=3, alpha=0.7)+
  facet_grid(test~R2_inst, labeller = labeller(R2_inst=r2labs))+
  coord_cartesian(ylim=c(0,1))+
  theme(text=element_text(size=11),
        strip.text.y = element_text(angle=0, face="bold"),
        legend.position = "bottom",
        legend.direction= "horizontal")+
  scale_x_continuous("Effect size (d)")+
  scale_y_continuous("Power") 
ggsave("./output/power4a.tiff")


#'
#'We can summarise by giving a) the smallest effect size for which we have
#'80 and 95% power in NHST and b) the largest effect size for which we have 
#'80 and 95% power to declare equivalence to zero:
# ----

# Smallest effect size for which 80% power (NHST):
power4a %>% 
  filter(power>.8,
         str_detect(test,"NHST")) %>%
  group_by(N, R2_inst) %>% 
  summarise(`min_es_80%power` = max(effect_size))%>% 
  ungroup()
print(paste0("sesoi = ",-r_to_d(unique(fullsim$ROPE_high_aim4))))
# Smallest effect size for which 95% power (NHST):
power4a %>% 
  filter(power>.95,
         str_detect(test,"NHST")) %>%
  group_by(N, R2_inst) %>% 
  summarise(`min_es_95%power` = max(effect_size))%>% 
  ungroup()
print(paste0("sesoi = ",-r_to_d(unique(fullsim$ROPE_high_aim4))))

# Largest effect size for which 80% power (equiv):
power4a %>% 
  filter(power>.8,
         str_detect(test,"NHST",negate=T)) %>%
  group_by(N, R2_inst) %>% 
  summarise(`max_es_80%power` = min(effect_size))%>% 
  ungroup()
print(paste0("sesoi = ",-r_to_d(unique(fullsim$ROPE_high_aim4))))
# Largest effect size for which 95% power (equiv):
power4a %>% 
  filter(power>.95,
         str_detect(test,"NHST",negate=T)) %>%
  group_by(N, R2_inst) %>% 
  summarise(`max_es_95%power` = min(effect_size))%>% 
  ungroup()
print(paste0("sesoi = ",-r_to_d(unique(fullsim$ROPE_high_aim4))))

#

#'## Power analysis 3/4b
#'
load(file="./output/fullsim_power34b.RData")

#'### Power 
#'Our power in NHST in each scenario is the proportion of iterations in which a true 
#'effect is detected, given the effect exists in the population
#'
#'Our power in equivalence/inferiority each scenario is the proportion of iterations 
#'in which an effect is declared equivalent to zero, given that in the population it is 
#'less than the SESOI
#'
# ----
#First, compare the two MR approaches
methodcomp <- fullsim %>% 
  group_by(effect_size,r2_inst, sample_size, case_rate, method) %>%
  summarise(bias = mean(est_COR-effect_size, na.rm=T),   
            mse = mean(est_COR-effect_size,na.rm=T)^2,   
            coverage = mean(effect_size > est_COR_lci90 & effect_size < est_COR_uci90 ,na.rm=T),
            precision = mean(est_COR_uci90 - est_COR_lci90, na.rm =T))

methodcomp %>% 
  group_by(method,case_rate) %>% 
  summarise(bias = mean(bias),  
            mse = mean(mse),  
            coverage = mean(coverage),
            precision = mean(precision) )

## logistic 2nd stage seems to be less biased, have lower mse and equivalent coverage
methodcomp %>% 
  group_by(method,sample_size) %>% 
  summarise(bias = mean(bias),  
            mse = mean(mse),  
            coverage = mean(coverage),
            precision = mean(precision) )
methodcomp %>% 
  group_by(method,effect_size) %>% 
  summarise(bias = mean(bias),  
            mse = mean(mse),  
            coverage = mean(coverage),
            precision = mean(precision) )

#Summarise at the level of scenario (i.e., across iterations)
power4b <- fullsim %>% 
  filter(method=="logistic_2nd_stage") %>% 
  group_by(effect_size,r2_inst, sample_size, case_rate) %>% 
  summarise(power_nhst= mean(p < 0.05),
            se_nhst=sd(p < 0.05)/sqrt(n()),
            power_equiv= mean(CI_low  > ROPE_low_aim4 & CI_high < ROPE_high_aim4, na.rm=T),
            se_equiv=sd(CI_low  > ROPE_low_aim4 & CI_high < ROPE_high_aim4, na.rm=T)/sqrt(n())) %>% 
  mutate(effect_size=oddsratio_to_d(effect_size )) %>% 
  rename(N = sample_size, prev=case_rate) %>% 
  pivot_longer(cols=matches("power|se_"), names_to = "type", values_to="power") %>% 
  separate(type, into=c("est","test")) %>% 
  pivot_wider(names_from=est,values_from=power) %>% 
  mutate( test=factor(test,
                      levels=c("nhst","equiv"),
                      labels=c("NHST\n(H1:Effect > 0)","Equivalence\n(H1:Effect = 0)")),
          N=factor(N),
          r2_inst=as.character(r2_inst))

r2labs <- c("R2 inst. 5%","R2 inst. 7.5%","R2 inst. 10%")
names(r2labs)<- c("0.05","0.075","0.1")

#Plot the power curves 
ggplot(data=power4b, aes(x= effect_size, y= power, fill=factor(prev),
                         colour=factor(prev),group=interaction(N,prev), shape=N))+
  geom_errorbar(aes(ymin=power-se,ymax=power+se), width=0, size=1.2, alpha=0.4)+
  geom_point(size=2.6,
             alpha=0.4)+
  geom_line(size=0.8)+
  geom_hline(yintercept=0.95, colour="red",size=1.4, linetype=2)+
  geom_hline(yintercept=0.80, colour="orange",size=1.4, linetype=3, alpha=0.7)+
  facet_grid(test~r2_inst, labeller = labeller(r2_inst= r2labs) )+
  coord_cartesian(ylim=c(0,1))+
  theme(text=element_text(size=11),
        strip.text.y = element_text(angle=0, face="bold"),
        legend.position = "bottom",
        legend.direction= "horizontal")+
  scale_x_continuous("Effect size (d)")+
  scale_y_continuous("Power") 
ggsave("./output/power4b.tiff")

#'We can summarise by giving a) the smallest effect size for which we have
#'80 and 95% power in NHST and b) the largest effect size for which we have 
#'80 and 95% power to declare equivalence to zero:
# ----

# Smallest effect size for which 80% power (NHST):
power4b %>% 
  filter(power>.8,
         str_detect(test,"NHST")) %>%
  group_by(N, r2_inst,prev) %>% 
  summarise(`min_es_80%power` = max(effect_size))%>% 
  ungroup()
print(paste0("sesoi = ",-oddsratio_to_d(exp(unique(fullsim$ROPE_high_aim4)))))
# Smallest effect size for which 95% power (NHST):
power4b %>% 
  filter(power>.95,
         str_detect(test,"NHST")) %>%
  group_by(N, r2_inst,prev) %>% 
  summarise(`min_es_95%power` = max(effect_size))%>% 
  ungroup()
print(paste0("sesoi = ",-oddsratio_to_d(exp(unique(fullsim$ROPE_high_aim4)))))

# Largest effect size for which 80% power (equiv):
power4b %>% 
  filter(power>.8,
         str_detect(test,"NHST",negate=T)) %>%
  group_by(N, r2_inst,prev) %>% 
  summarise(`max_es_80%power` = min(effect_size))%>% 
  ungroup() 
print(paste0("sesoi = ",-oddsratio_to_d(exp(unique(fullsim$ROPE_high_aim4)))))
# Largest effect size for which 95% power (equiv):
power4b %>% 
  filter(power>.95,
         str_detect(test,"NHST",negate=T)) %>%
  group_by(N, r2_inst,prev) %>% 
  summarise(`max_es_95%power` = min(effect_size))%>% 
  ungroup()

print(paste0("sesoi = ",-oddsratio_to_d(exp(unique(fullsim$ROPE_high_aim4)))))

# ---- eval=F,include=F
#This is the code to render and knit the script - called automatically by 
#generate_code_walkthrough()
# ezknitr::ezspin("scripts/template.R", 
#                 out_dir="scripts/reports",
#                 keep_rmd = T, 
#                 verbose=T)

NA
NA
NA
NA
NA
NA
NA
NA
NA
