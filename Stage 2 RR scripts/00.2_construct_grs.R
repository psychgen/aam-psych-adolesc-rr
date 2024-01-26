#00.2_construct_grs.R

##This script constructs the genetic instrument for age at menarche for
##one-sample MR analyses, regressing out covariates and standardising.
##Note that some of the script needs to be run on the cluster initially.

#load required packages
library(tidyverse)
library(data.table)
library(genotools)

#create age at menarche (aam) genetic risk score (GRS)

#quality control of summary statistics (remove INDELS, reformat Marker column and create CHR)

aam_sstats <- fread("./data/sstats/Menarche_1KG_NatGen2017_WebsiteUpload.txt",data.table=FALSE)

aam_sstats_subset <- aam_sstats %>%
  dplyr::filter(!str_detect(Markername,"INDEL")) %>%
  mutate(CHR = sapply(strsplit(Markername, ":"), '[', 1)) %>%
  mutate(CHR = str_remove(CHR, "chr")) %>%
  mutate(BP = sub(".*:","", Markername)) %>%
  select(-Markername)

bim <- fread("N:/durable/data/genetic/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc.bim",data.table=FALSE)

bim <- bim %>%
  select("CHR" = V1, "SNP" = V2, "BP" = V4) %>%
  mutate(CHR = as.character(CHR),
         BP = as.character(BP))

aam_sstats_subset <- aam_sstats_subset %>%
  left_join(bim) %>%
  filter(!is.na(SNP))

write.table(aam_sstats_subset, file = "./data/sstats/Menarche_1KG_NatGen2017_WebsiteUpload_subset.txt", sep = "\t", quote=FALSE,row.names=FALSE,col.names=TRUE)

#create submission script based on publicly available summary statistics:
###https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip###
#information about columns in the sumstats from README file:

#Markername - chr:position (build 37)
#Allele 1 - effect allele
#Allele 2 - other allele
#Effect - beta measured in years (from the final meta-analysis excluding 23andMe, max N ~ 252k).
#Pvalue - p-value for association with menarche (from the final meta-analysis excluding 23andMe, max N ~ 252k)
#minor_allele - minor allele based on study meta-analysis.

make_prsice(user= "Adrian",
            jobname = "aam",
            outcome_type="continuous",
            cpu_time="08:00:00",
            memory="32G",
            inputs_dir="/cluster/p/p471/cluster",
            outputs_dir="/cluster/p/p471/cluster/common/raw_prsice_output/",
            sumstats_filename = "Menarche_1KG_NatGen2017_WebsiteUpload_subset.txt",
            genotype_dir="data/genetic_data/MoBaPsychGen_v1",
            genotype_data="MoBaPsychGen_v1-ec-eur-batch-basic-qc",
            prsice_dir="/cluster/p/p471/cluster/common/prsice/",
            A1 ="Allele1",
            A2 = "Allele2",
            stat = "Effect",
            pvalue = "Pvalue",
            snp = "SNP",
            bp = "BP",
            chr = "CHR",
            thresholds = "5e-08",
            clump_kb = "10000",
            clump_p = "1.000000",
            clump_r2 = "0.001",
            lower = "5e-08",
            maf = "0.01",
            mhc= "exclude",
            add_metadata = TRUE,
            meta_dir = "//ess01/p471/data/durable/common/pgs_directory/add_to_dir/",
            source = "https://www.reprogen.org/Menarche_1KG_NatGen2017_WebsiteUpload.zip",
            pmid ="10.1038/ng.3841")


pgs <- fetch_pgs(c("aam"),
                 thresholds = c(5e-08),
                 threshold_range = FALSE,
                 geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                 pgs_software = "prsice2",
                 pgs_meta = genotools::pgs_metadata,
                 maf = "0.01",
                 clump = "10000_1_0")

pgs_procd <- process_pgs(pgs,
                         geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                         covs_dir = "//ess01/P471/data/durable/data/genetic/MoBaPsychGen_v1",
                         regress_pgs = FALSE,
                         return_covs = TRUE)

pgs_nsnps <- get_pgs_nsnps(c("aam"),
                           thresholds=c(5e-08, 1),
                           geno_data = "MoBaPsychGen_v1-ec-eur-batch-basic-qc",
                           pgs_directory="//ess01/P471/data/durable/common/pgs_directory/pgs",
                           pgs_software="prsice2",
                           maf="0.01",
                           clump="10000_1_0")
pgs_nsnps

#remove parents
pgs_procd <- pgs_procd %>%
  filter(Role == "Child")

#regress on batch and 20 PCs manually
m <- lm(pgs_procd$`aam_p<5e-08` ~ 
          pgs_procd$PC1 + pgs_procd$PC2 + pgs_procd$PC3 + pgs_procd$PC4 + pgs_procd$PC5 + 
          pgs_procd$PC6 + pgs_procd$PC7 + pgs_procd$PC8 + pgs_procd$PC9 + pgs_procd$PC10 +
          pgs_procd$PC11 + pgs_procd$PC12 + pgs_procd$PC13 + pgs_procd$PC14 + pgs_procd$PC15 + 
          pgs_procd$PC16 + pgs_procd$PC17 + pgs_procd$PC18 + pgs_procd$PC19 + pgs_procd$PC20 +
          pgs_procd$genotyping_batch, na.action = na.exclude)

pgs_procd$aam_res <- rstandard(m)

pgs_procd <- pgs_procd %>%
  mutate(preg_id = as.numeric(preg_id)) %>%
  mutate(BARN_NR = as.numeric(BARN_NR)) %>%
  select(preg_id,BARN_NR,"PGS_AaM"=aam_res)

#save genetic data
save(pgs_procd, file="data/pgs_procd.RData")

