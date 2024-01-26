# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Power analyses for hypothesis/equivalence tests in Aim 1 (`02.1_power_analyses_aim1.R`) 
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("{021_power_analyses_aim1_codemap.png}", error = FALSE)
#'
#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to summarise the power analyses for the 
#'hypothesis/equivalence tests in Aim 1. 
#'
#'
#'&nbsp;
#'&nbsp;
# ----
load(file="./output/fullsim_power1a.RData")
library(tidyverse)
library(effectsize)
#'
#'### Power 
#'Our power in NHST in each scenario is the proportion of iterations in which a true 
#'effect is detected, given the effect exists in the population
#'
#'Our power in equivalence/inferiority each scenario is the proportion of iterations 
#'in which an effect is declared equivalent to zero, given that in the population it is 
#'less than the SESOI
#'
# ----
#Compute one-sided p-value
fullsim <- fullsim %>% 
  mutate(p = ifelse(Estimate<0,`Pr(>|t|)`/2,1))
#Summarise at the level of scenario (i.e., across iterations)
power1a <- fullsim %>% 
  group_by(corr, sample_size) %>% 
  summarise(power_nhst= mean(p < 0.05),
            se_nhst=sd(p < 0.05)/sqrt(n()),
            power_equiv= mean(CI_low  > ROPE_low),
            se_equiv=sd(CI_low  > ROPE_low)/sqrt(n())) %>% 
  rename(effect_size= corr,N = sample_size) %>% 
  mutate(effect_size=r_to_d(effect_size )) %>% 
  pivot_longer(cols=matches("power|se_"), names_to = "type", values_to="power") %>% 
  separate(type, into=c("est","test")) %>% 
  pivot_wider(names_from=est,values_from=power) %>% 
  mutate( test=factor(test,
                     levels=c("nhst","equiv"),
                     labels=c("NHST\n(H1:Effect > 0)","Equivalence\n(H1:Effect = 0)")),
         N=factor(N))%>% 
  filter(N%in%c("12000","13000"))
#'
# ----
#Plot the power curves 
ggplot(data=power1a, aes(x= effect_size, y= power, fill=N,
                         colour=N,group=N))+
  geom_errorbar(aes(ymin=power-se,ymax=power+se), width=0, size=1.2, alpha=0.4)+
  geom_point(size=2.6,
             alpha=0.4)+
  geom_line(size=0.8)+
  geom_hline(yintercept=0.95, colour="red",size=1.4, linetype=2)+
  geom_hline(yintercept=0.80, colour="orange",size=1.4, linetype=3, alpha=0.7)+
  facet_grid(test~.)+
  coord_cartesian(ylim=c(0,1))+
  theme(text=element_text(size=11),
        strip.text.y = element_text(angle=0, face="bold"),
        legend.position = "bottom",
        legend.direction= "horizontal")+
  scale_x_continuous("Effect size (d)")+
  scale_y_continuous("Power")
ggsave("./output/power1a.tiff")

#'
#'We can summarise by giving a) the smallest effect size for which we have
#'80 and 95% power in NHST and b) the largest effect size for which we have 
#'80 and 95% power to declare equivalence to zero:
# ----

# Smallest effect size for which 80% power (NHST):
power1a %>% 
  filter(power>.8,
         str_detect(test,"NHST")) %>%
  group_by(N) %>% 
  summarise(`min_es_80%power` = max(effect_size))%>% 
  mutate(sesoi = rep(-oddsratio_to_d(1.4),3))
# Smallest effect size for which 95% power (NHST):
power1a %>% 
  filter(power>.95,
         str_detect(test,"NHST")) %>%
  group_by(N) %>% 
  summarise(`min_es_95%power` = max(effect_size))%>% 
  mutate(sesoi = rep(-oddsratio_to_d(1.4),3))

# Largest effect size for which 80% power (equiv):
power1a %>% 
  filter(power>.8,
         str_detect(test,"NHST",negate = T)) %>%
  group_by(N) %>% 
  summarise(`max_es_80%power` = min(effect_size))%>% 
  mutate(sesoi = rep(-oddsratio_to_d(1.4),3))
# Largest effect size for which 95% power (equiv):
power1a %>% 
  filter(power>.95,
         str_detect(test,"NHST",negate=T)) %>%
  group_by(N) %>% 
  summarise(`max_es_95%power` = min(effect_size))%>% 
  mutate(sesoi = rep(-oddsratio_to_d(1.4),3))

#
#'## Power analysis 1b
# ----
load(file="./output/fullsim_power1b.RData")
#'
#'### Power 
#'Our power in NHST in each scenario is the proportion of iterations in which a true 
#'effect is detected, given the effect exists in the population
#'
#'Our power in equivalence/inferiority each scenario is the proportion of iterations 
#'in which an effect is declared equivalent to zero, given that in the population it is 
#'less than the SESOI
#'
# ----

fullsim <- fullsim %>% 
  mutate(p = ifelse(Estimate<0,`Pr(>|z|)`/2,1))
#Summarise at the level of scenario (i.e., across iterations)
power1b <- fullsim %>% 
  group_by(effect_size, sample_size,case_rate) %>% 
  summarise(power_nhst= mean(p < 0.05),
            se_nhst=sd(p < 0.05)/sqrt(n()),
            power_equiv= mean(CI_low  > ROPE_low),
            se_equiv=sd(CI_low  > ROPE_low)/sqrt(n())) %>% 
  rename(N = sample_size, prevalence = case_rate) %>% 
  mutate(effect_size=oddsratio_to_d(effect_size )) %>% 
  pivot_longer(cols=matches("power|se_"), names_to = "type", values_to="power") %>% 
  separate(type, into=c("est","test")) %>% 
  pivot_wider(names_from=est,values_from=power) %>% 
  mutate( test=factor(test,
                      levels=c("nhst","equiv"),
                      labels=c("NHST\n(H1:Effect > 0)","Equivalence\n(H1:Effect = 0)")),
          N=factor(N)) %>%  
  filter(N%in%c("12000","13000"))
#'
# ----
#Plot the power curves 

prevlabs <- c("2% avg. \nprevalence","6% avg. \nprevalence","10% avg. \nprevalence")
names(prevlabs)<- c("0.02","0.06", "0.1")
ggplot(data=power1b, aes(x= effect_size, y= power, fill=N,
                         colour=N,group=N))+
  geom_errorbar(aes(ymin=power-se,ymax=power+se), width=0, size=1.2, alpha=0.4)+
  geom_point(size=2.6,
             alpha=0.4)+
  geom_line(size=0.8)+
  geom_hline(yintercept=0.95, colour="red",size=1.4, linetype=2)+
  geom_hline(yintercept=0.80, colour="orange",size=1.4, linetype=3, alpha=0.7)+
  facet_grid(test~prevalence, labeller=labeller(prevalence=prevlabs))+
  coord_cartesian(ylim=c(0,1))+
  theme(text=element_text(size=11),
        strip.text.y = element_text(angle=0, face="bold"),
        legend.position = "bottom",
        legend.direction= "horizontal")+
  scale_x_continuous("Effect size (d)")+
  scale_y_continuous("Power")
ggsave("./output/power1b.tiff")


#'We can summarise by giving a) the smallest effect size for which we have
#'80 and 95% power in NHST and b) the largest effect size for which we have 
#'80 and 95% power to declare equivalence to zero:
# ----

# Smallest effect size for which 80% power (NHST):
power1b %>% 
  filter(power>.8,
         str_detect(test,"NHST")) %>%
  group_by(N, prevalence) %>% 
  summarise(`min_es_80%power` = max(effect_size))%>% 
  mutate(sesoi = rep(-oddsratio_to_d(1.4),3))
# Smallest effect size for which 95% power (NHST):
power1b %>% 
  filter(power>.95,
         str_detect(test,"NHST")) %>%
  group_by(N, prevalence) %>% 
  summarise(`min_es_95%power` = max(effect_size))%>% 
  mutate(sesoi = rep(-oddsratio_to_d(1.4),3))

# Largest effect size for which 80% power (equiv):
power1b %>% 
  filter(power>.8,
         str_detect(test,"NHST",negate = T)) %>%
  group_by(N, prevalence) %>% 
  summarise(`max_es_80%power` = min(effect_size))%>% 
  mutate(sesoi = rep(-oddsratio_to_d(1.4),3))
# Largest effect size for which 95% power (equiv):
power1b %>% 
  filter(power>.95,
         str_detect(test,"NHST",negate=T)) %>%
  group_by(N, prevalence) %>% 
  summarise(`max_es_95%power` = min(effect_size))%>% 
  mutate(sesoi = rep(-oddsratio_to_d(1.4),3))

