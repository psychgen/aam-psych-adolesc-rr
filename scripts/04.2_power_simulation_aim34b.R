# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
#'
#'# Power analyses for aim 3(&4) hypothesis/equivalence tests (`04.2_power_simulation_aim34b.R`) 
#'
# ---- echo=F, out.width='33%', out.extra = 'style="float:right;"'
knitr::include_graphics("042_power_simulation_aim34b_codemap.png", error = FALSE)
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
#'## Power analysis 3/4b
#'
#'### Simulation 
#'
#'Here we create a function which simulates both the 14 yr age at menarche (AAM)
#'and (depression cases), based on:
#'1. The effect size
#'2. 14-yr data availability (as a proportion of 8-year availability)
#'3. Prevalence of depression in sample
#'...and runs a linear model and accompanying equivalence test, saving 
#'standardised betas and inferences based on decision rules in a data.frame
# ----

set.seed(14516)


sim_dx_mr <- function(causal_eff_d=-0.2,
                      R2_inst=0.075, 
                      sample_size=10500, 
                      iteration=1, 
                      equiv_d_aim3=0.1, 
                      equiv_d_aim4=0.1,
                      case_rate=0.06){
  
  # Mean and SD values for AaM are from  10.1192/bjp.bp.115.168617 Sequeira et al.
  # Supp table DS1
  
  
  
  
    #Calculate OR_yz based on COR and R2_inst:
    
    b_xz = sqrt(R2_inst)
    Codds = d_to_oddsratio(causal_eff_d, log=T)
    odds_yz = b_xz * Codds
    
    confound= rnorm(1, 0, d_to_oddsratio(0.05,log=T))
    #if(confound<0){confound=0} #Allow masking?
    
    b_yx = Codds + confound
    
    
    simdat <- rnorm_multi(n = sample_size, 
                          mu = c(151.52, 0),
                          sd = c(14.11, 1),
                          r = b_xz,# rnorm(1, b_xz, 0.02), 
                          varnames = c("AaM", "PGS_AaM"),
                          empirical = TRUE) %>% 
      mutate(AaM = round(AaM/12)) 
    
    # Note that the possibility of AaM values > 168 represents the situation
    # after we have imputed right censored values
    
    # Get PGS predicted values of AaM
    
    res <-lm(scale(AaM) ~ PGS_AaM, data = simdat) 
    
    simdat <- simdat %>% 
      add_predictions(res) %>% 
      add_residuals(res)
    
    
    # Simulate the dx variable
    beta0 <- log((1/(1-case_rate))-1)

    
    beta1 <- b_yx 
    pi_x <- exp(beta0 + beta1 * scale(simdat$AaM) ) / (1 + exp(beta0 + beta1 * scale(simdat$AaM)))
    simdat$"dx" <- rbinom(n=length(simdat$AaM), size=1, prob=pi_x)
    obs_case_rate <- prop.table(table(simdat$dx))[[2]]
    
    #Run the model both ways, adjust SEs and get 90%CIs for equivalence testing (will use bootstrapped
    #CIs in the analysis)
    
    mr_res_logs2<- glm(dx ~ pred, family = binomial(link="logit"), data=simdat)
    
    mr_res_ivreg <- ivreg(dx ~ scale(AaM)|PGS_AaM, data=simdat)
     
    #And we also need to adjust the SEs to account for the uncertainty in the first stage:
 
    # Calculate robust standard errors (from https://stats.stackexchange.com/
    # questions/89999/how-to-replicate-statas-robust-binomial-glm-for-proportion-data-in-r)

    adjSEs <- function(x){
      cov.m1 <- vcovHC(x, type = "HC1") 
      ## This is used as it mirrors statas ivreg2 robust, used in Sequeira et al)
      std.err <- sqrt(diag(cov.m1))
      q.val <- qnorm(0.95)
      comb_res <- cbind(
        Estimate = coef(x)
        , "Robust SE" = std.err
        , z = (coef(x)/std.err)
        , "p"= 2 * pnorm(abs(coef(x)/std.err), lower.tail = FALSE)
        , CI_low = coef(x) - q.val  * std.err
        , CI_high = coef(x) + q.val  * std.err
      ) %>% 
        as.data.frame() %>% 
        .[2,] %>% 
        mutate(ROPE_high_aim3 = d_to_oddsratio(equiv_d_aim3, log=T),
               ROPE_low_aim3 =- d_to_oddsratio(equiv_d_aim3, log=T),
               ROPE_high_aim4 = d_to_oddsratio(equiv_d_aim4, log=T),
               ROPE_low_aim4 =- d_to_oddsratio(equiv_d_aim4, log=T))
      return(comb_res)
    }
    
    convert_ivreg_out <- function(x){
      y_x0 = coef(mr_res_ivreg)[[1]]
      y_x1 = coef(mr_res_ivreg)[[1]] + x
      odds <- log((y_x1 * (1-y_x0)) / ((1-y_x1)* y_x0))
      return(odds)
    }
    
    mr_log_res <- adjSEs(mr_res_logs2) %>% 
      mutate(method="logistic_2nd_stage")
    mr_ivr_res <- adjSEs(mr_res_ivreg) %>%
      mutate(across(c("Estimate","CI_low","CI_high"), ~convert_ivreg_out(.x)),
             method="standard_2SLS")
    
    
    #Combine and apply decision rules - inferiority, i.e., only use high bound
    
    comb_res <- mr_log_res %>% 
      bind_rows(mr_ivr_res) %>% 
      mutate(p_onetailed = ifelse(Estimate<0,p/2,1), #One-tailed test p-value computation
             iteration = iteration,
             r2_inst = R2_inst,
             effect_size = d_to_oddsratio(causal_eff_d),
             est_COR = exp(Estimate),
             est_COR_lci90 = exp(CI_low),
             est_COR_uci90 = exp(CI_high),
             sample_size=sample_size,
             case_rate=case_rate,
             obs_case_rate=obs_case_rate)
    
    
    simres_full<- comb_res
    
    
  row.names(simres_full)<- NULL
  return(simres_full)
}
#'
#'We apply this function with different parameters for the second power analysis
#'
# ----

# Set range of parameters
r2s <-c(0.05,0.075,0.10) # start, end, increment
ns <- c(9500,10500)
ces <- seq(0.0,-0.30,-0.02)
rates <- c(0.02,0.06,.10,.14)

equiv_d_aim3=0.25
equiv_d_aim4=0.20

nsims = 1000

## Load the checkpointing file, if it exists:
checkFile <- "checkpoint_aim34b.RData"
tempFile <- "tempCheckpoint_aim34b.RData"
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
    for(r2 in r2s){
      for(n in ns){
        for(ce in ces){
          for(rate in rates){
            message(paste0(
              "\n R2 = ",r2,
              "\n N = ",n,
              "\n CE = ",ce,
              "\n rate = ",rate,
              "\n (replicate ",iter,")"))
            scen_temp <- sim_dx_mr(causal_eff_d=ce,
                                                    R2_inst=r2, 
                                                    sample_size=n, 
                                                    iteration=iter, 
                                                    equiv_d_aim3=equiv_d_aim3,
                                                    equiv_d_aim4=equiv_d_aim4,
                                                    case_rate= rate )
            fullsim <- rbind(fullsim, scen_temp)
            
            
            ## Save the results of the iteration.  (By first saving to a temporary file
            ## and then renaming it into the checkpointing file, we guard against being
            ## interrupted while saving.)
            save.image(tempFile)
            file.rename(tempFile, checkFile)
          }
        }
      }
    }
    
    dur <- proc.time()-ptm
    durs <-c(durs,dur[3])
    message(paste0("\nProjected finish in ", 
                   round( mean(durs)*(nsims-iter)/60,2), " mins" ))
  }
  
} 


save(fullsim, file="./output/fullsim_power34b.RData")
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
