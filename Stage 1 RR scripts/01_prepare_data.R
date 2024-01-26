# ---- include = F
knitr::opts_knit$set(root.dir = rprojroot::find_rstudio_root_file())
knitr::opts_chunk$set(eval = FALSE)
#'
# ---- include = F
#'
#'# Data prepraration for the project (`01_prepare_data.R`) 
#'
# ---- echo=F, eval=T, out.width='33%',  fig.align='center'
knitr::include_graphics("01_prepare_data_codemap.png", error = FALSE)
#'

#'
#'
#'***
#'
#'&nbsp;
#'
#'The purpose of this script is to source and prepare variables for analysis. 
#'
#'Most of the preparation of raw variables in the MoBa dataset and linked
#'registry sources is done using the `phenotools` R package, code and documentation
#'for which is available at https://github.com/psychgen/phenotools.
#'
#'
#'&nbsp;
#'&nbsp;
#'
#'## MoBa/MBRN variables required:
#'
#'### Phenotypic:
#'
#'  * self-rep age at menarche and pubertal stage 14 year;
#'  * depressive sx (sMFQ) 8 & 14year; 
#'  * anxiety sx (SCARED) 8 & 14 year; 
#'  * conduct disorder sx (RS-DBD) 8 & 14 year;
#'  * oppositional defiant disorder (RS-DBD) sx 8 & 14 year;
#'  * ADHD sx (RS-DBD) 8 & 14 yr;
#'  * Covariates: BMI (derived) 8 & 14 yr; maternal/paternal age; child age at
#'  * Q completion; parental education; parental income; parental cohabitation;
#'  * number of children in family; financial problems; maternal prenatal 
#'  * depression (SCL); maternal postnatal depression(EPD)
#'
#'### Genetic:
#'
#'  * age at menarche PRS (GWSig hits only) 
#'  * childhood BMI PRS (GWSig hits only) 10.1093/hmg/ddv472
#'  * adult BMI PRS (GWSig hits only) 10.1093/hmg/ddy271
#'
#'
#'## Registry variables required:
#'
#'### Primary care diagnoses (KUHR; ICPC-2):
#'  * depressive disorders: P76 
#'  * anxiety disorders: P74, P79, P82
#'  * adhd: P81
#'  * conduct/behavioural disorders: P23 
#'  
#'### Secondary care diagnoses (NPR; ICD-10):
#'  * depressive disorders: F32-F33, F34.1 
#'  * anxiety disorders: F40-F44, F93.0-F93.2
#'  * adhd: F90
#'  * conduct/behavioural disorders: F91, F92   
#'  
#'***
#'
#'## Load required packages
# ---- 
library(phenotools)
library(tidyverse)
library(genotools)
#'
#'
#'
#'## Curate datatset with variables from MoBa phenotypic data
# ---- collapse = TRUE

covars <- c("KJONN",
            "age_at_q_ret_8yr",
            "parent_income_derived_q1",
            "parent_educ_derived_q1",
            "parent_cohab_18m",
            "parent_cohab_3yr",
            "p_age_at_birth",
            "m_age_at_birth",
            "financ_probs_18m",
            "scl_5item_m_q1",
            "scl_full_m_q3",
            "epds_full_m_6m",
            "n_children_8yr",
            "bmi_derived_c_8yr",
            "bmi_derived_c_14s")

q8yr_vars <- c("smfq_dep_c_8yr",
               "scared_anx_c_8yr",
               "rsdbd_adhd_c_8yr",
               "rsdbd_cd_c_8yr",
               "rsdbd_odd_c_8yr")

q14yr_vars <- c("aam_c_14s",
                "smfq_dep_c_14s",
                "scared_anx_c_14s",
                "rsdbd_adhd_c_14m",
                "rsdbd_cd_c_14s",
                "rsdbd_odd_c_14m")

if(!file.exists("data/pheno_data.RData")){
  pheno_data <- curate_dataset(
    variables_required=list(moba=c(covars,
                                   q8yr_vars,
                                   q14yr_vars),
                            npr=c("dep=F32,F33,F341",
                                  "anx=F40,F41,F42,F43,F44,F930,F931,F932",
                                  "adhd=F90",
                                  "beh=F91,F92"),
                            kuhr=c("dep=P76",
                                   "anx=P74, P79, P82",
                                   "adhd=P81",
                                   "beh=P23")),
    out_format = "merged_df",
    recursive=T,
    dx_owners="child")
  pheno_data <- pheno_data %>% 
    rename("sex"=KJONN_raw)
  filter(sex ==2)
  
  save(pheno_data, file="data/pheno_data.RData")
}else{
  load(file="data/pheno_data.RData")
}
#'
#'NB: all 14-year variables  are (by design) not available during 
#'preparation/submission of Stage 1 RR
#'
#'## Fetch and process polygenic scores
# ---- eval=FALSE
#Not run
pgs <- fetch_prs(c("aam","bmi","cbmi"),
                 thresholds = c( 1),
                 threshold_range = FALSE,
                 geno_data = "moba_full",
                 maf = "0.01",
                 clump = "10000_1_0.001")

pgs_procd <- process_prs(pgs,
                         geno_data = "moba_full",
                         regress_prs = TRUE,
                         return_covs = TRUE)
#'NB: polygenic scores generated using PRSice will be read in to R and processed 
#'in this section using functions from the `genotools` package, as outlined 
#'above . 
#'
#'
#'## Join genetic and phenotypic data for complete analytic dataset
# ---- 

fulldata <- pheno_data %>% 
  left_join(pgs_procd)


#save out
save(fulldata, file="data/analysis_dataset.RData")

# ---- eval=F,include=F
#This is the code to render and knit the script
ezknitr::ezspin("scripts/01_prepare_data.R",
                out_dir="scripts/reports",
                keep_rmd = T,
                verbose=T)

