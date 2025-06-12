library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
pwy_subset <- args[1]


liver_bms <- c("alt_log", "ast_log", "ggt_log")
collections <- c("kegg_legacy", "hallmark", "kegg_medicus")
pwy_prefixes <- c(kegg_legacy = "KEGG", hallmark = "HALLMARK", kegg_medicus = "KEGG_MEDICUS")


calc_n_eff <- function(y, collection) {
  print(paste0("Running: ", y, " / ", collection, "..."))
  ppgs_fn <- paste0("../data/processed/pgs/", y, "_prset_", collection, ".all_score")
  ppgs_df <- read_delim(ppgs_fn, delim = " ", show_col_types = FALSE) %>%
    select(id = IID, contains("_0.001")) %>%
    rename_with(~ gsub("_0.001", "", .x))
  
  if (pwy_subset == "all") {
    pwy_vec <- setdiff(names(ppgs_df), c("id", "Base"))
  } else if (pwy_subset == "sig") {
    pwy_vec <- readLines(paste0("../data/processed/pgs/", y, "_", collection, "_sig_pwy_names.txt"))
    pwy_vec <- paste0(pwy_prefixes[collection], "_", pwy_vec)
  }
  
  pca <- ppgs_df %>%
    select(all_of(pwy_vec)) %>%
    #mutate(across(everything(), ~ ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))) %>%
    prcomp(scale. = TRUE)
  eigenvals <- pca$sdev ** 2
  n_eff <- sum(eigenvals) ** 2 / sum(eigenvals ** 2)
  n_eff
}


n_eff_tbl <- expand_grid(y = liver_bms, pathway_group = collections) %>%
  rowwise() %>%
  mutate(n_eff = calc_n_eff(y, pathway_group)) %>%
  ungroup()

write_csv(n_eff_tbl, paste0("../data/processed/pgs/", pwy_subset, "_n_eff_tbl.csv"))
