library(tidyverse)


collection <- "kegg_legacy"


ukb_phenos <- read_csv("../data/processed/ukb_phenos_unrelated.csv", show_col_types = FALSE) %>%
  mutate(across(matches("^gPC"), ~ . * bmi, .names = "bmiBy{.col}"))

gPCs <- paste0("gPC", 1:10)
primary_covars <- c("sex", "age", "age_squared", "ageBySex", gPCs)


run_lrt <- function(y) {

  print(paste0("Running ", y, "..."))

  ppgs_fn <- paste0("../data/processed/pgs/", y, "_prset_", collection, ".all_score")
  ppgs_df <- read_delim(ppgs_fn, delim = " ", show_col_types = FALSE) %>%
    select(id = IID, contains("_0.001")) %>%
    rename_with(~ gsub("_0.001", "", .x))

  pwy_vec <- setdiff(colnames(ppgs_df), c("id", "Base"))

  regression_df <- inner_join(ukb_phenos, ppgs_df, by = "id")

  null_fmla <- paste0(
    y, " ~ ",
    paste(primary_covars, collapse = " + "), " + ",
    paste0("bmiBy", gPCs, collapse = " + "), " + ",
    "Base * bmi + ",
    paste(pwy_vec, collapse = " + ")
  )
  alt_fmla <- paste0(
    null_fmla, " + ", 
    paste0(pwy_vec, " * bmi", collapse = " + ")
  )

  null_fit <- lm(null_fmla, data = regression_df)
  alt_fit <- lm(alt_fmla, data = regression_df)

  lmtest::lrtest(null_fit, alt_fit)
}

liver_bms <- c("alt_log", "ast_log", "ggt_log")

for (y in liver_bms) {
  lrt <- run_lrt(y)
  saveRDS(lrt, paste0("../data/processed/pgs/", y, "_", collection, "_pwy_vs_base_lrt.rds"))
}
