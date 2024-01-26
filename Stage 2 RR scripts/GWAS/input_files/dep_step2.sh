#!/bin/bash
#SBATCH --job-name=dep_step2
#SBATCH --account=p471_tsd
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=25G
#SBATCH --cpus-per-task=18



source /cluster/bin/jobsetup
module purge
module load singularity/3.7.3

set -o errexit

OUTM="/ess/p471/cluster/projects/aam_RR/GWAS/output_files"

export SINGULARITYENV_LC_ALL=C

singularity exec -B /ess/p471/cluster /ess/p471/cluster/regenie_v3.1.1.gz.sif regenie \
--step 2 \
--bed /ess/p471/cluster/data/genetic_data/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc \
--covarFile /ess/p471/cluster/projects/aam_RR/GWAS/input_files/covar_gwas.txt \
--phenoFile /ess/p471/cluster/projects/aam_RR/GWAS/input_files/pheno_gwas_dep_tr.txt \
--bsize 200 \
--strict \
--pred /ess/p471/cluster/projects/aam_RR/GWAS/output_files/dep_pred.list \
--out ${OUTM}/dep_step2