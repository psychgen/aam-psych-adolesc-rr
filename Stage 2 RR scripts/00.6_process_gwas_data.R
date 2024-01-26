#00.6_process_gwas_data.R

##This script does processing of the MoBa GWAS summary data, saving
##them out for 2-sample MR analyses, and generates Manhattan/QQ-plots. 

#load required packages
library(qqman)
library(dplyr)
library(data.table)
library(ggplot2)

#check top hits from GWAS
dep <- fread("N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/Transformed sx GWAS/smfq_dep_c_14c.txt") 
anx <- fread("N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/Transformed sx GWAS/scared_anx_c_14c.txt") 
cd <- fread("N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/Transformed sx GWAS/rsdbd_cd_c_14c.txt") 
odd <- fread("N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/Transformed sx GWAS/rsdbd_odd_c_14m.txt") 
adhd <- fread("N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/Transformed sx GWAS/rsdbd_adhd_c_14m.txt") 

#exponentiate p-values and drop columns
dep <- dep %>% mutate(P = 10^-LOG10P) %>%
  select(-c(TEST,LOG10P,EXTRA))
anx <- anx %>% mutate(P = 10^-LOG10P) %>%
  select(-c(TEST,LOG10P,EXTRA))
cd <- cd %>% mutate(P = 10^-LOG10P) %>%
  select(-c(TEST,LOG10P,EXTRA))
odd <- odd %>% mutate(P = 10^-LOG10P) %>%
  select(-c(TEST,LOG10P,EXTRA))
adhd <- adhd %>% mutate(P = 10^-LOG10P) %>%
  select(-c(TEST,LOG10P,EXTRA))

#check hits
sig_dep <- dep %>% filter(P < 5e-08)
sig_anx <- anx %>% filter(P < 5e-08)
sig_cd <- cd %>% filter(P < 5e-08)
sig_odd <- odd %>% filter(P < 5e-08)
sig_adhd <- adhd %>% filter(P < 5e-08)

#write out
write_delim(dep,file="N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/smfq_dep_c_14c.txt",col_names=TRUE)
write_delim(anx,file="N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/scared_anx_c_14c.txt",col_names=TRUE)
write_delim(cd,file="N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/rsdbd_cd_c_14c.txt",col_names=TRUE)
write_delim(odd,file="N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/rsdbd_odd_c_14m.txt",col_names=TRUE)
write_delim(adhd,file="N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/rsdbd_adhd_c_14m.txt",col_names=TRUE)

#make Manhattan and QQ plots for the Regenie GWAS
setwd("N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/")

Sys.getenv(c("DISPLAY"))
options(bitmapType='cairo') 

files <- list.files(path = "N:/durable/projects/age_at_menarche_RR/scripts/GWAS/output_files/", pattern = ".txt")

GWAS_plots <- function(file){
  name <- stringr::str_remove(file, ".txt")

# manhattan plot
results_lin <- fread(file, head=TRUE)
jpeg(paste0("Linear_manhattan_",name,".jpeg"), 
     width = 1024, height = 768, pointsize = 14, quality = 100)
manhattan(results_lin, chr="CHROM",bp="GENPOS",p="P",snp="ID", main = paste0("Manhattan plot: ", name), ylim = c(0,10),  cex = 0.6, cex.axis = 0.9, col = c("navy", "royalblue"), suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), chrlabs = c(1:22))
dev.off()

# QQ plot:
jpeg(paste0("QQPlot_", name, ".jpeg") , 
     width = 600, height = 400, pointsize = 20, quality = 100)
qq(results_lin$P)
dev.off()

}

purrr::map(files, GWAS_plots)
