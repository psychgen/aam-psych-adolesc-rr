# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Power analyses for hypothesis/equivalence tests in Aim 1 (`02.2_power_simulation_aim1b.R`) 
#'
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("{022_power_simulation_aim1b_codemap}.png", error = FALSE)
#'

#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to run power analyses for the hypothesis/equivalence
#'tests in Aim 1. Specifically, we will:
#'  
#'  * Simulate a depressive symptoms outcome variable based on a range of 
#'  scenarios in which the association with (simulated) age at menarche increases
#'  from 0.01 to 0.1 in increments of 0.01, and in which availability of 14-year
#'  data varies between 12k, 13k, and 14k
#'  
#'  * Run linear models on the simulated data for 1000 iterations of each
#'  scenario, deriving power empirically as the proportion of significant 
#'  (#'  at alpha = 5%) effects detected in each scenario  (power analysis 1a)
#'  
#'  * Simulate equivalent data and run equivalent power analyses for the 
#'  generalized linear models of depression diagnoses, additionally varying
#'  prevalence of depression  (power analysis 1b)
#'
#'
#'&nbsp;
#'&nbsp;
#'
#'## Load required packages
# ---- 

library(tidyverse)
library(faux)
library(effectsize)
library(parameters)
#'
#'## Power analysis 1b
#'
#'### Simulation 
#'
#'Here we create a function which simulates both the 14 yr age at menarche (AAM)
#'and depression cases, based on:
#'1. The effect size
#'2. 14-yr data availability (as a proportion of 8-year availability)
#'3. Prevalence of depression in sample
#'...and runs a linear model and accompanying equivalence test, saving 
#'standardised betas and inferences based on decision rules in a data.frame
# ----

set.seed(82928)


sim_dep_dx <- function(effect_size_d, 
                       sample_size, 
                       case_rate, 
                       iteration, 
                       equiv_d){
  
  # Mean and SD values for AaM are from  10.1192/bjp.bp.115.168617 Sequeira et al.
  # Supp table DS1
  
  or <- d_to_oddsratio(effect_size_d)
  
  simdat <- tibble("AaM" =  rnorm(n = sample_size,
                                  m = 151.52,
                                  sd = 14.11))
  
  # The new way, the downside of which is that the only "error" is in the measurement
  # of AaM, but since we specify beta1 for AaM as it is "measured", this is not propagated
  # to the models - the net result being that we always make the same inference
  # in each iteration of a given scenario (not realistic)
  
  # To avoid this, we will add noise to the beta1 parameter
  
  
  beta0 <- log((1/(1-case_rate))-1)
  beta1 <- log(or) + rnorm(1, mean=0,sd=0.10)
  pi_x <- exp(beta0 + beta1 * scale(simdat$AaM)) / (1 + exp(beta0 + beta1 * scale(simdat$AaM)))
  simdat$"DepDx" <- rbinom(n=length(simdat$AaM), size=1, prob=pi_x)
  
  obs_case_rate <- prop.table(table(simdat$DepDx))[[2]]
  
  
  # Note that the possibility of AaM values > 168 represents the situation
  # after we have imputed right censored values
  
  
  #Run the model
  
  res <-glm(DepDx ~ scale(AaM), family = binomial(link="logit"), data=simdat)
  
  #Extract the results of the NHST
  
  nhst_res <- res %>% 
    summary() %>% 
    .$coefficients %>%
    as.data.frame()%>% 
    `row.names<-`(NULL) %>% 
    .[2,] %>% 
    mutate(oddsr = exp(Estimate),  # Odds ratio/gradient
           var.diag = diag(vcov(res))[2],  # Variance of each coefficient
           or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 
  
  #Run the equivalence test (using OR = 1.4 as SESOI)
  
  equiv_res <- parameters::equivalence_test(
    res, range= c(log(d_to_oddsratio(-equiv_d)),log(d_to_oddsratio(equiv_d))), rule="classic")  %>% 
    as.data.frame() %>% 
    `row.names<-`(NULL) %>% 
    .[2,]
  
  #Combine and apply decision rules - inferiority, i.e., only use high bound
  
  comb_res <- nhst_res %>% 
    bind_cols(equiv_res) %>% 
    mutate(p = ifelse(Estimate<0,`Pr(>|z|)`/2,1), #One-tailed test p-value computation
           iteration = iteration,
           effect_size = or,
           sample_size=sample_size,
           case_rate=case_rate,
           obs_case_rate=obs_case_rate) # We monitor this because, at a certain effect size,
  # the case rate gets inflated for a given sample size
  
  
  simres_full<- comb_res
  
  
  
  row.names(simres_full)<- NULL
  return(simres_full)
}
#'
#'We apply this function with different parameters for the second power analysis
#'
# ----

# Set range of parameters
es <- seq(0.0,-0.26,-0.02)
ns <- c(12000,13000)
rates <- c(0.02,0.06,.10,.14)

equiv=oddsratio_to_d(1.4)

nsims=1000

## Load the checkpointing file, if it exists:
checkFile <- "checkpoint_aim1b.RData"
tempFile <- "tempCheckpoint_aim1b.RData"
if (file.exists(checkFile)) {
  load(checkFile)
  cat("Resuming after iteration", iter, "\n")
  iter <- iter + 1 
}else{
  # Create some objects to make the timer work
  iter<- 1
  durs <- vector()
  # Create an empty DF for the output
  fullsim<- data.frame()
}

if(iter <= nsims){
  
  for(iter in iter:nsims){
    
    
    message(paste0(
      "\nReplicate ",iter, " of ", nsims
    ))
    ptm <- proc.time()
    
    # Loop over the function
    for(e in es){
      for(n in ns){
        for(c in rates){
          scen_temp <- suppressMessages(sim_dep_dx(effect_size_d=e,
                                                   sample_size=n, 
                                                   case_rate = c,
                                                   equiv_d=equiv,
                                                   iteration=iter))
          fullsim <- rbind(fullsim, scen_temp)
          ## Save the results of the iteration.  (By first saving to a temporary file
          ## and then renaming it into the checkpointing file, we guard against being
          ## interrupted while saving.)
          save.image(tempFile)
          file.rename(tempFile, checkFile)
        }
      }
    }
    
    dur <- proc.time()-ptm
    durs <-c(durs,dur[3])
    message(paste0("\nProjected finish in ", 
                   round( mean(durs)*(nsims-iter)/60,2), " mins" ))
  }
}

save(fullsim, file="./output/fullsim_power1b.RData")
## Clean up (_after_ saving the results)
if (file.exists(checkFile)) file.remove(checkFile)
NA
NA
NA
NA
NA
NA
