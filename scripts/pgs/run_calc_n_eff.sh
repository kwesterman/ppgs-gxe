#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=30G
#$ -l h_rt=2:00:00

#$ -j y
#$ -cwd


pwy_subset=$1


source /broad/software/scripts/useuse
use R-4.1

Rscript pgs/calc_n_eff.R ${pwy_subset}
