# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Power analyses for aim 2 hypothesis/equivalence tests (`03.1_power_simulation_aim2a.R`) 
#'
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("{031_power_simulation_aim2a_codemap}.png", error = FALSE)
#'

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
#'## Power analysis 2a
#'
#'### Simulation
#'
#'Here we create a function which simulates both the 14 yr symptoms
#'variables and age at menarche (AAM) variable based on two inputs:
#'1. The average correlation (across all domains) between with aam
#'2. The 14-yr data availability (as a proportion of 8-year availability)
#'...and runs a linear model and accompanying equivalence test, saving 
#'standardised betas and inferences based on decision rules in a data.frame
# ----


set.seed(91734)


sim_all_sx <- function(corr=0.1,  
                       sample_size=13000, 
                       iteration=1, 
                       equiv_d=0.1){
  
  # Mean and SD values are from 10.1192/bjp.bp.115.168617 Sequeira et al.
  # Supp table DS1 (aam) and observed 8year data in MoBa:
  
  #anx m= 1.03 sd= 1.20
  #cd m= 0.78 sd= 1.51
  #odd m= 3.42 sd= 3.16
  #adhd m= 8.54  sd= 7.23
  
simdat <- rnorm_multi(n = sample_size, 
                          mu = c(151.52, 1.03),
                          sd = c(14.11, 1.2),
                          r = rnorm(1, corr, 0.02), 
                          varnames = c("AaM", "AnxSx14yr"),
                          empirical = FALSE) %>% 
      mutate(AnxSx14yr = ifelse(AnxSx14yr<0,0,round(AnxSx14yr)),
             #Assume a floor effect; this attenuates the correlation slightly
             CDsx14yr = rnorm_pre(AaM, mu= 0.78,sd=1.51 ,
                                  r =rnorm(1, corr, 0.02), empirical=FALSE ),
             ODDsx14yr = rnorm_pre(AaM, mu= 3.42,sd= 3.16,
                                   r =rnorm(1, corr, 0.02), empirical=FALSE ),
             ADHDsx14yr = rnorm_pre(AaM, mu= 8.54,sd= 7.23,
                                    r =rnorm(1, corr, 0.02), empirical=FALSE ),
             CDsx14yr = ifelse(CDsx14yr<0,0,round(CDsx14yr)),
             #Assume a floor effect; this attenuates the correlation slightly
             ODDsx14yr = ifelse(ODDsx14yr<0,0,round(ODDsx14yr)),
             #Assume a floor effect; this attenuates the correlation slightly
             ADHDsx14yr = ifelse(ADHDsx14yr<0,0,round(ADHDsx14yr)),
             #Assume a floor effect; this attenuates the correlation slightly
             AaM = round(AaM/12) )
    
    
    # Note that the possibility of AaM values > 168 represents the situation
    # after we have imputed right censored values
    
    #Run the models 
    
    res <-list(
      lm(scale(AnxSx14yr)~ scale(AaM), data = simdat), 
      lm(scale(CDsx14yr)~ scale(AaM), data = simdat),
      lm(scale(ODDsx14yr)~ scale(AaM), data = simdat),
      lm(scale(ADHDsx14yr) ~ scale(AaM), data = simdat) )
    
    #Extract the results of the NHST
    
    nhst_res <- map(res,function(x){
      x %>% 
        summary() %>% 
        .$coefficients %>%
        as.data.frame()%>% 
        `row.names<-`(NULL) %>% 
        .[2,]  
    })%>% 
      reduce(bind_rows) %>% 
      rename(p=`Pr(>|t|)`)
    
    #Run the equivalence test 
    
    equiv_res <- map(res, function(x){
      x %>% 
        parameters::equivalence_test( range= c(d_to_r(-equiv_d),d_to_r(equiv_d)), 
                                      rule="classic",
                                      ci=1-(0.05/nrow(nhst_res)),p_values=T)  %>% 
        as.data.frame() %>% 
        select(-p) %>% 
        `row.names<-`(NULL) %>% 
        .[2,]
    }) %>% 
      reduce(bind_rows)
    
    
    most_extreme <- max(abs(c(equiv_res$CI_low,equiv_res$CI_high)))
    
    equiv_res <- equiv_res %>% 
      filter(abs(CI_high)==most_extreme|
               abs(CI_low)==most_extreme)  
    # We only need take the result of the equivalence test for the most extreme parameter estimate

    #Combine and apply decision rules
    
    comb_res <- nhst_res %>% 
      bind_cols(equiv_res) %>% 
      mutate(iteration = iteration,
             corr = corr,
             sample_size=sample_size)
    
    simres_full<- comb_res

  row.names(simres_full)<- NULL
  return(simres_full)
}

#'
#'We apply this function with different parameters for the first power analysis
#'
# ----

# Set range of parameters
corrs <- c(seq(-0.01,-0.15,-0.01)) #start, end, increment
ns <- seq(12000,14000,1000)
equiv <- oddsratio_to_d(1.49) 


nsims=1000

## Load the checkpointing file, if it exists:
checkFile <- "checkpoint_aim2a.RData"
tempFile <- "tempCheckpoint_aim2a.RData"
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
    for(r in corrs){
      for(n in ns){
        
        scen_temp <- sim_all_sx(corr=r,
                                sample_size=n, 
                                iteration=iter,
                                equiv_d= equiv)
        fullsim <- rbind(fullsim, scen_temp)
        ## Save the results of the iteration.  (By first saving to a temporary file
        ## and then renaming it into the checkpointing file, we guard against being
        ## interrupted while saving.)
        save.image(tempFile)
        file.rename(tempFile, checkFile)
      }
    }
    dur <- proc.time()-ptm
    durs <-c(durs,dur[3])
    message(paste0("\nProjected finish in ", 
                   round( mean(durs)*(nsims-iter)/60,2), " mins" ))
  }
}

save(fullsim, file="./output/fullsim_power2a.RData")
## Clean up (_after_ saving the results)
if (file.exists(checkFile)) file.remove(checkFile)

