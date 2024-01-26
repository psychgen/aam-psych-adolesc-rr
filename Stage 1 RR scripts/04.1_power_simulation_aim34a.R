# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Power analyses for aim 3(&4) hypothesis/equivalence tests (`04.1_power_simulation_aim34a.R`) 
#'
# ---- echo=F, out.width='33%', out.extra = 'style="float:right;"'
knitr::include_graphics("041_power_simulation_aim34a_codemap.png", error = FALSE)
#'
#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to run the power analyses corresponding to 
#'Aim 3(&4) hypothesis/equivalence tests. 
#'Specifically, we will:
#'
#'  * Simulate a polygenic score of genome-wide significant SNPs explaining 
#'  variance in (simulated) aam at different levels, with causal associations 
#'  with depression symptoms simulated based on SNP-predicted values of AaM
#'  which associations with aam are randomly drawn from a normal distribution
#'  with an incrementally increasing mean, and in which availability of 14-year
#'  data (that is, genotyped sample size for the analyses) varies between 9.5k,
#'  10.5k, and 11.5k
#'  
#'  * Run one-sample MR (2 SLS linear models) on simulated data for 1000 iterations 
#'  per scenario, deriving power empirically for the NHST/Equivalence tests (alpha = 5%)
#'  NB we use two sets of equivalence bounds, one for depression based on the "small 
#'  telescopes" procedure described in the report, and one generic one for outcomes in
#'  any other domain. This is because if the criteria for triggering an Aim 4
#'  analysis are met by any of the domain-specific effects from Aim 2 analyses,
#'  the MR will be identical to the depression outcome MR
#'  
#'  * Simulate depression diagnosis outcome and run one-sample MR (structural
#'  mean models)  on simulated data for 1000 iterations 
#'  per scenario, deriving power empirically for the NHST/Equivalence tests (alpha = 5%)
#'  NB we use two sets of equivalence bounds, one for depression based on the "small 
#'  telescopes" procedure described in the report, and one generic one for outcomes in
#'  any other domain. This is because if the criteria for triggering an Aim 4
#'  analysis are met by any of the domain-specific effects from Aim 2 analyses,
#'  the MR will be identical to the depression outcome MR
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
library(AER)
library(modelr)
library(broom)
library(sandwich)


#'
#'## Power analysis 3/4a
#'
#'### Simulation
#'
#'Here we create a function which simulates the age at menarche (AAM) genetic
#'instrument and observed value, and depressive symptoms based on three inputs:
#'1. The explanatory power of the instrument in the exposure (aam ~ pgs)
#'2. 
#'2. The 14-yr data availability (as a proportion of 8-year availability)
#'...and runs a linear model and accompanying equivalence test, saving 
#'standardised betas and inferences based on decision rules in a data.frame
# ----

set.seed(91734)


sim_sx_mr <- function( causal_eff_d=0.25,
                       R2_inst=0.09, 
                       sample_size=10500, 
                       iteration=1, 
                       equiv_d_aim3=0.1,
                       equiv_d_aim4=0.1,
                       out_mean=5.71,
                       out_sd=4.93){
  
  # Mean and SD values are from  10.1192/bjp.bp.115.168617 Sequeira et al.
  # Supp table DS1
    
    #Calculate b_yz based on ACE and R2_inst:
    
    b_xz = sqrt(R2_inst)
    ace = rnorm(1, d_to_r(causal_eff_d), 0.02)
    b_yz = b_xz * ace
    
    confound= rnorm(1, 0, d_to_r(0.05))
    #if(confound<0){confound=0} #Allow masking?
    
    b_yx = ace + confound
    
    
    simdat <- rnorm_multi(n = sample_size, 
                          mu = c(151.52, 0, out_mean ),
                          sd = c(14.11, 1, out_sd ),
                          r =  c(b_xz, b_yx, b_yz) , 
                          varnames = c("AaM", "PGS_AaM","sx"),
                          empirical = TRUE) %>% 
      mutate(AaM = round(AaM/12),
             sx = ifelse(sx<0,0,round(sx)))
    #Assume a floor effect; this attenuates the correlation slightly 
    
     #Run the model and get 90%CIs for equivalence testing (will use bootstrapped
    #CIs in the analysis)
    
    mr_res <- ivreg(scale(sx) ~ scale(AaM) | scale(PGS_AaM), data = simdat) %>% 
      broom::tidy(conf.int=T, conf.level=0.9) %>% 
      mutate(ROPE_high_aim3 = d_to_r(equiv_d_aim3),
             ROPE_low_aim3 =- d_to_r(equiv_d_aim3),
             ROPE_high_aim4 = d_to_r(equiv_d_aim4),
             ROPE_low_aim4 =- d_to_r(equiv_d_aim4)) %>% 
      .[2,]
    
    #Combine and apply decision rules - inferiority (i.e., only use high equiv bound)
    
    comb_res <- mr_res %>% 
      mutate(p_one_tailed = ifelse(estimate<0,p.value/2,1), #One-tailed test p-value computation
             iteration = iteration,
             R2_inst = R2_inst,
             causal_eff_d=causal_eff_d,
             true_causal=ace,
             sample_size=sample_size)
    
    
  row.names(comb_res)<- NULL
  return(comb_res)
}
#'
#'We apply this function with different parameters for the second power analysis
#'
# ----

# Set range of parameters
r2s <- c(0.05,0.075,0.10) 
ns <- c(9500,10500)
ces <- seq(0.0,-0.30,-0.01)

equiv_d_aim3=0.25
equiv_d_aim4=0.20

nsims=1000

## Load the checkpointing file, if it exists:
checkFile <- "checkpoint_aim34a.RData"
tempFile <- "tempCheckpoint_aim34a.RData"
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
    
    for(r2 in r2s){
      for(n in ns){
        for(ce in ces){
          message(paste0(
            "\n R2 = ",r2,
            "\n N = ",n,
            "\n CE = ",ce,
            "\n (replicate ",iter,")"))
          scen_temp <- suppressMessages(sim_sx_mr(causal_eff_d=ce,
                                                  R2_inst=r2, 
                                                  sample_size=n, 
                                                  iteration=iter, 
                                                  equiv_d_aim3=equiv_d_aim3,
                                                  equiv_d_aim4=equiv_d_aim4,
                                                  out_mean=5.71,  
            # Mean and SD values are from  10.1192/bjp.bp.115.168617 Sequeira et al.
                                                  out_sd=4.93))  # Supp table DS1
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


save(fullsim, file="./output/fullsim_power34a.RData")
## Clean up (_after_ saving the results)
if (file.exists(checkFile)) file.remove(checkFile)

NA
NA
NA
NA
NA
NA
NA
NA
NA
