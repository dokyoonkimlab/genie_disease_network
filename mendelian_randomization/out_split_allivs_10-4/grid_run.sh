#!/bin/bash
#$ -N mrs
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=15:00:00
#$ -t 1-135

module unload gcc/4.9.2
R_LIBS=/cbica/home/manu148/R/x86_64-redhat-linux-gnu-library/3.6
Rscript ../run_MR_split_allIVs.R "/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/continous,/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/discrete" "/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/continous_split1,/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/discrete_split1" "/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/continous_split2,/cbica/projects/kimlab_hpeace/projects/genie_phewas/QC/imputed/association1/plato/discrete_split2" ../out_split_allivs_10-4 0.0001 F F $SGE_TASK_ID F

