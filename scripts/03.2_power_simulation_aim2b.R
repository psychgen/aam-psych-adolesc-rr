# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Power analyses for aim 2 hypothesis/equivalence tests (`03.2_power_simulation_aim2b.R`) 
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("{032_power_simulation_aim2b_codemap.png}", error = FALSE)
#'
#'
#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to run the power analyses corresponding to 
#'Aim 2 hypothesis tests. 
#'Specifically, we will:
#'
#'  * Simulate outcomes in four new symptom domains (anx, cd, odd, adhd) in
#'  which associations with aam are randomly drawn from a normal distribution
#'  with an incrementally increasing mean, and in which availability of 14-year
#'  data (that is, sample size for the analyses) varies (as a percentage of the
#'  observed 8-year sample size) from 55-65% in increments of 5%
#'  
#'  * Run linear models on simulated data for 1000 iterations per scenario, 
#'  deriving power empirically for the composite NHS+Equivalence test of a null
#'  hypothesis that *all* effects are practically equivalent to zero (alpha = 5%)
#'  
#'  * Run a logistic regression model with aam on other diagnoses (
#'  anx, cd, odd, adhd) and use an omnibus test in conjunction with equivalence
#'  testing to test of a null hypothesis that *all* effects are practically 
#'  equivalent to zero (alpha = 5%)
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
library(faux)

#'
#'## Power analysis 2b
#'
#'### Simulation 
#'
#'Here we create a function which simulates both the 14 yr age at menarche (AAM)
#'and anx/cd/odd/adhd cases, based on:
#'1. The effect size
#'2. 14-yr data availability (as a proportion of 8-year availability)
#'3. Prevalence of anx/cd/odd/adhd in sample
#'...and runs a generalized linear model and accompanying equivalence test, saving 
#'standardised betas and inferences based on decision rules in a data.frame
# ----

set.seed(82928)

#First, load the real data:

sim_oth_dx <- function(effect_size_d=-0.10, 
                       sample_size=13000, 
                       case_rate=0.06, 
                       iteration=1, 
                       equiv_d=oddsratio_to_d(1.49)){
  
  # Mean and SD values for AaM are from  10.1192/bjp.bp.115.168617 Sequeira et al.
  # Supp table DS1
  
  # Convert OR to r as we are simulating the underlying continuous liability
  #corr <- effectsize::convert_oddsratio_to_r(or, log=F)
  
  or <- d_to_oddsratio(effect_size_d)  
  
  simdat <- tibble("AaM" =  rnorm(n = sample_size,
                                  m = 151.52,
                                  sd = 14.11))
  
  # The new way, the downside of which is that the only "error" is in the measurement
  # of AaM, but since we specify beta1 for AaM as it is "measured", this is not propagated
  # to the models - the net result being that we always make the same inference
  # in each iteration of a given scenario (not realistic)
  
  # To avoid this, we will add noise to the beta1 parameter
  
  # Here, we also vary the case rate, with a lower bound of 1%
  cr_anx<- case_rate+ rnorm(1, mean=0,sd=0.02)
  cr_anx<-ifelse(cr_anx<0.01,0.01,cr_anx)
  cr_cdodd<- case_rate+ rnorm(1, mean=0,sd=0.02)
  cr_cdodd<-ifelse(cr_cdodd<0.01,0.01,cr_cdodd)
  cr_adhd <-case_rate+ rnorm(1, mean=0,sd=0.02)
  cr_adhd <-ifelse(cr_adhd<0.01,0.01,cr_adhd)
  
  #Anx
  beta0 <- log((1/(1-cr_anx))-1)
  beta1 <- log(or) + rnorm(1, mean=0,sd=0.05)
  pi_x <- exp(beta0 + beta1 * scale(simdat$AaM)) / (1 + exp(beta0 + beta1 * scale(simdat$AaM)))
  simdat$"AnxDx" <- rbinom(n=length(simdat$AaM), size=1, prob=pi_x)
  #Conduct/ODD
  beta0 <- log((1/(1-cr_cdodd))-1)
  beta1 <- log(or) + rnorm(1, mean=0,sd=0.05)
  pi_x <- exp(beta0 + beta1 * scale(simdat$AaM)) / (1 + exp(beta0 + beta1 * scale(simdat$AaM)))
  simdat$"CdOddDx" <- rbinom(n=length(simdat$AaM), size=1, prob=pi_x)
  #ADHD
  beta0 <- log((1/(1-cr_adhd))-1)
  beta1 <- log(or) + rnorm(1, mean=0,sd=0.05)
  pi_x <- exp(beta0 + beta1 * scale(simdat$AaM)) / (1 + exp(beta0 + beta1 * scale(simdat$AaM)))
  simdat$"AdhdDx" <- rbinom(n=length(simdat$AaM), size=1, prob=pi_x)
  
  simdat$AaM <- scale(simdat$AaM)
  
  #Run the models
  
  res <-
    list(glm(AnxDx~ AaM, family = binomial(link="logit"), data=simdat),
         glm(CdOddDx~ AaM, family = binomial(link="logit"), data=simdat),
         glm(AdhdDx~ AaM, family = binomial(link="logit"), data=simdat))
  
  #Extract the results of the NHST
  
  nhst_res <- map(res,function(x){
    x %>% 
      summary() %>% 
      .$coefficients %>%
      as.data.frame()%>% 
      `row.names<-`(NULL) %>% 
      .[2,] %>% 
      mutate(oddsr = exp(Estimate),  # Odds ratio/gradient
             var.diag = diag(vcov(x))[2],  # Variance of each coefficient
             or.se = sqrt(or^2 * var.diag))  # Odds-ratio adjusted 
  })%>% 
    reduce(bind_rows) %>% 
    rename(p=`Pr(>|z|)`)
  
  #Run the equivalence test 
  
  equiv_res <- map(res, function(x){
    x %>% 
      parameters::equivalence_test( 
        range= c(d_to_oddsratio(-equiv_d, log=T),d_to_oddsratio(equiv_d, log=T)), 
        rule="classic",
        ci=1-(0.05/nrow(nhst_res)),
        p_values=T )  %>% 
      as.data.frame() %>% 
      `row.names<-`(NULL) %>% 
      .[2,]
  }) %>% 
    reduce(bind_rows)
  
  #Combine and apply decision rules
  
  comb_res <- nhst_res %>% 
    bind_cols(equiv_res) %>% 
    mutate(iteration = iteration,
           effect_size = or,
           sample_size=sample_size,
           case_rate=case_rate)
  most_extreme <- max(abs(c(equiv_res$CI_low,equiv_res$CI_high)))
  
  comb_res <- comb_res %>% 
    filter(abs(CI_high)==most_extreme|
             abs(CI_low)==most_extreme)  
  # We only need take the result of the equivalence test for the most extreme parameter estimate
  
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
ns <- c(12000,13000,14000)
rates <- c(0.02,0.06,.10)

equiv=oddsratio_to_d(1.49)


nsims=1000

## Load the checkpointing file, if it exists:
checkFile <- "checkpoint_aim2b.RData"
tempFile <- "tempCheckpoint_aim2b.RData"
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
  # Loop over the function
  for(iter in iter:nsims){
    
    
    message(paste0(
      "\nReplicate ",iter, " of ", nsims
    ))
    ptm <- proc.time()
    
    
    # Loop over the function
    for(e in es){
      for(n in ns){
        for(c in rates){
          
          scen_temp <- suppressMessages(sim_oth_dx(effect_size_d=e,
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
    message(paste0("\nProjected finish in ", round(mean(durs)*(nsims-iter)/60,2), " mins" ))
  }
}

save(fullsim, file="./output/fullsim_power2b.RData")
## Clean up (_after_ saving the results)
if (file.exists(checkFile)) file.remove(checkFile)
NA
NA
NA
