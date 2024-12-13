#!/bin/bash


#$ -l h_vmem=15G
#$ -l h_rt=01:00:00

#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use R-4.1


Rscript pgs/optimize_ppgs.R
