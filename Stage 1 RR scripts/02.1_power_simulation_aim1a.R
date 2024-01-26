# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Power analyses for hypothesis/equivalence tests in Aim 1 (`02.1_power_simulation_aim1a.R`) 
#'
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("{021_power_simulation_aim1a_codemap}.png", error = FALSE)
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
#'## Power analysis 1a
#'
#'### Simulation
#'
#'Here we create a function which simulates both the 14 yr depressive symptoms
#'variable and age at menarche (AAM) variable based on two inputs:
#'1. The correlation between the two
#'2. The 14-yr data availability (12k, 13k, 14k)
#'...and runs a linear model and accompanying equivalence test, saving 
#'standardised betas and inferences based on decision rules in a data.frame
# ----

set.seed(91734)


sim_dep_sx <- function(corr=0.1, 
                       sample_size=13000, 
                       iteration=1, 
                       equiv_d=0.1){
  
  # Mean and SD values are from  10.1192/bjp.bp.115.168617 Sequeira et al.
  # Supp table DS1
  
  simdat <- rnorm_multi(n = sample_size, 
                        mu = c(151.52, 5.71),
                        sd = c(14.11, 4.93),
                        r = corr, 
                        varnames = c("AaM", "DepSx14yr"),
                        empirical = FALSE) %>% 
    mutate(DepSx14yr = ifelse(DepSx14yr<0,0,round(DepSx14yr)),
           #Assume a floor effect; this attenuates the correlation slightly
           AaM = round(AaM/12) ) 
  
  # Note that the possibility of AaM values > 168 represents the situation
  # after we have imputed right censored values
  
  #Run the model 
  
  res <-lm(scale(DepSx14yr) ~ scale(AaM), data = simdat) 
  
  #Extract the results of the NHST
  
  nhst_res <- res %>% 
    summary() %>% 
    .$coefficients %>%
    as.data.frame()%>% 
    `row.names<-`(NULL) %>% 
    .[2,]
  
  #Run the equivalence test (using OR = 1.4 as SESOI)
  
  equiv_res <- parameters::equivalence_test(
    res, range= c(-effectsize::d_to_r(equiv_d),effectsize::d_to_r(equiv_d)), rule="classic")  %>% 
    as.data.frame() %>% 
    `row.names<-`(NULL) %>% 
    .[2,]
  
  
  #Combine and apply decision rules - inferiority (i.e., 
  #only use high equiv bound (negative in this case))
  
  comb_res <- nhst_res %>% 
    bind_cols(equiv_res) %>% 
    mutate(p = ifelse(Estimate<0,`Pr(>|t|)`/2,1), #One-tailed test p-value computation
           iteration = iteration,
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
corrs <- c(seq(0,-0.15,-0.01)) #start, end, increment
ns <- seq(12000,14000,1000)
equiv <- oddsratio_to_d(1.4) 


nsims=1000

## Load the checkpointing file, if it exists:
checkFile <- "checkpoint_aim1a.RData"
tempFile <- "tempCheckpoint_aim1a.RData"
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
        
        scen_temp <- sim_dep_sx(corr=r,
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

save(fullsim, file="./output/fullsim_power1a.RData")
## Clean up (_after_ saving the results)
if (file.exists(checkFile)) file.remove(checkFile)

