library(tidyverse)
library(parallel)


ukb_test_df <- read_csv(paste0("../data/processed/ukb_testing_set.csv"), show_col_types = FALSE) %>%
  mutate(across(matches("^gPC"), ~ . * bmi, .names = "bmiBy{.col}"))

exposures <- c("bmi")
exposures_clean <- c("BMI")

outcomes <- c("alt_log")#, "ast_log", "ggt_log")
outcomes_clean <- c("log(ALT)")#, "log(AST)", "log(GGT)")

pathway_groups <- c("kegg_legacy")
pathway_groups_clean <- c("KEGG Legacy")

thresholds <- c("0.001")
threshold_names <- paste0("Pt_", c("0.001"))
threshold_names_clean <- c("P<0.001")

sources <- c("full", "held_out")

ppgs_config_df <- expand_grid(
  e = exposures,
  y = outcomes,
  pathway_group = pathway_groups,
  threshold = thresholds,
  source = sources
)

gPCs <- paste0("gPC", 1:10)
primary_covars <- c("sex", "age", "age_squared", "ageBySex", gPCs)
liver_bms <- c("alt_log")

read_ppgs_full <- function(tag, pg, threshold) {
  ppgs_df <- read_delim(paste0("../data/processed/pgs/", tag, "_prset_", pg, ".all_score"),
                        delim = " ", show_col_types = FALSE) %>%
    select(id = IID, contains(paste0("_", threshold))) %>%
    rename_with(~ gsub("_[0-9\\.e]+$", "", .x))
  ppgs_df
}

read_ppgs_trainset <- function(tag, pg) {
  ppgs_df <- read_delim(paste0("../data/processed/pgs/", tag, "_prset_", pg, ".best"),
                        delim = " ", show_col_types = FALSE) %>%
    select(id = IID, everything(), -FID, -In_Regression)
  ppgs_df
}

test_single_ppgs_gxe <- function(e, y, regression_df, covars, standardize) {
  if (standardize) {
    regression_df <- mutate(regression_df,
                            across(all_of(c(e, y, "ppgs")), ~ as.vector(scale(.x))))
  }
  lm_form_str <- paste0(y, " ~ ppgs * ", e)
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  lm_fit <- lm(as.formula(lm_form_str), data = regression_df)
  lm_fit %>%
    broom::tidy() %>%
    filter(term %in% c(paste0("ppgs:", e), "ppgs", e))
}

test_all_ppgs_gxe <- function(e, y, pg, threshold, source, covars, 
                              standardize = TRUE) {
  print(paste0("Testing ", e, " - ", y, " in the held-out testing set using ", 
               pg, " pathways and ", threshold, " threshold and ", source, " pPGS source."))
  if (source == "full") {
    ppgs_df <- read_ppgs_full(y, pg, threshold)
  } else if (source == "held_out") {
    ppgs_df <- read_ppgs_trainset(paste0(y, "_trainset"), pg)
  } else {
    stop("Source must be 'full' or 'held_out'")
  }
  pathways <- setdiff(names(ppgs_df), "id")
  regression_df <- inner_join(ukb_test_df, ppgs_df, by = "id")
  mclapply(pathways, function(p) {
    regression_df$ppgs <- regression_df[[p]]
    test_single_ppgs_gxe(e, y, regression_df, covars, standardize)
  }, mc.cores = 8) %>%
    setNames(pathways) %>%
    bind_rows(.id = "pathway")
}

ppgs_res_df <- ppgs_config_df %>%
  select(e, y, pathway_group, threshold, source) %>%
  rowwise() %>%
  mutate(pathway_fits = list(
    test_all_ppgs_gxe(e, y, pathway_group, threshold, source,
                      c(primary_covars, paste0(e, "By", gPCs)))
  )) %>%
  unnest(pathway_fits)

saveRDS(ppgs_res_df, "../data/processed/pgs/ppgs_res_df_trainset.rds")
