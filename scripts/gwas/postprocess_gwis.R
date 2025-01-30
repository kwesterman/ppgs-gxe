library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
filepath <- args[1]
e <- args[2]
y <- args[3]

gwis_dir <- dirname(filepath)

print(paste0("Postprocessing: ", filepath))


### Read in summary stats and subset to columns of interest

ss_cols <- c(
  CHR = "CHR", SNP = "RSID", POS = "POS",
  EA = "Effect_Allele", NEA = "Non_Effect_Allele", AF = "AF",
  N = "N_Samples", P_int = "P_Value_Interaction", P_joint = "P_Value_Joint", 
  robust_P_int = "robust_P_Value_Interaction", robust_P_joint = "robust_P_Value_Joint"
)

high_qual_variants <- read_tsv("../data/processed/ukb_rsIDs_maf0.005_info0.5.txt", 
                               col_names=F, col_types="c")[[1]]

ss_df <- fread(filepath, stringsAsFactors=F, data.table=F) %>%
  select(all_of(ss_cols), matches("^Beta"), matches("^Var_Beta")) %>%
  mutate(across(contains("P_"), ~ as.numeric(.))) %>%
  filter(SNP %in% high_qual_variants)

ss_df %>%
  write_tsv(gsub("_merged", "_merged_clean", filepath))
ss_df %>%
  filter(P_int < 0.05) %>%
  write_tsv(gsub("_merged", "_merged_clean_nom", filepath))
