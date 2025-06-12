#!/bin/bash


#$ -l h_vmem=10G
#$ -l h_rt=6:00:00

#$ -cwd
#$ -j y


source /broad/software/scripts/useuse
use GCC-5.2
use R-4.0


tag=$1
#sumstats_file=$1
#prefix=$2


magma_dir=../opt/magma_v1.10
ldref_dir=../../ipgs/data/processed/ld_ref/
working_dir=../data/processed/magma
prsice_datadir=../data/processed/prsice

sumstats_file=../data/processed/gwas/${tag}_merged_clean


# Gene analysis (based on p-values)

${magma_dir}/magma \
	--bfile ${ldref_dir}/ukb_20k_hg19 \
	--pval ${sumstats_file} use=SNP,robust_P_int N=360000 \
	--gene-model snp-wise=mean \
	--gene-annot ${working_dir}/ukb_20k_hg19_2.1.genes.annot \
	--out ${working_dir}/${tag}

# Pathway analysis (based on gene-level results)
${magma_dir}/magma \
    --gene-results ${working_dir}/${tag}.genes.raw \
    --set-annot ${prsice_datadir}/c2.cp.kegg_legacy.v2024.1.Hs.entrez.gmt \
    --out ${working_dir}/${tag}
