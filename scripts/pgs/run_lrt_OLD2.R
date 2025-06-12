library(tidyverse)


collection <- "kegg_legacy"


ukb_phenos <- read_csv("../data/processed/ukb_phenos_unrelated.csv", show_col_types = FALSE) %>%
  mutate(across(matches("^gPC"), ~ . * bmi, .names = "bmiBy{.col}"))

gPCs <- paste0("gPC", 1:10)
primary_covars <- c("sex", "age", "age_squared", "ageBySex", gPCs)


run_lrt <- function(y, type) {

  print(paste0("Running ", y, " - ", type, "..."))

  ppgs_fn <- paste0("../data/processed/pgs/", y, "_prset_", collection, ".all_score")
  ppgs_df <- read_delim(ppgs_fn, delim = " ", show_col_types = FALSE) %>%
    select(id = IID, contains("_0.001")) %>%
    rename_with(~ gsub("_0.001", "", .x))

  stopifnot(type %in% c("gw", "pwy"))
  if (type == "gw") {
    pwy_vec <- "Base"
  } else if (type == "pwy") {
    pwy_vec <- setdiff(colnames(ppgs_df), c("id", "Base"))
  }

  regression_df <- inner_join(ukb_phenos, ppgs_df, by = "id")

  null_fmla <- paste0(
    y, " ~ ",
    paste(primary_covars, collapse = " + "), " + ",
    paste0("bmiBy", gPCs, collapse = " + "), " + ",
    "bmi + ",
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
  lrt_gw <- run_lrt(y, type = "gw")
  saveRDS(lrt_gw, paste0("../data/processed/pgs/", y, "_", collection, "_gw_lrt.rds"))
  lrt_pwy <- run_lrt(y, type = "pwy")
  saveRDS(lrt_pwy, paste0("../data/processed/pgs/", y, "_", collection, "_pwy_lrt.rds"))
}
