#!/bin/sh


#$ -l os=RedHat7
#$ -l h_vmem=5G
#$ -l h_rt=6:00:00


#$ -pe smp 8
#$ -binding linear:8
#$ -R y

#$ -j y
#$ -cwd


source /broad/software/scripts/useuse
use R-4.1

Rscript pgs/test_ppgs_by_e.R
