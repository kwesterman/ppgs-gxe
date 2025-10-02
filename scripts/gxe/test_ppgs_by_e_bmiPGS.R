library(tidyverse)
library(parallel)


ukb_phenos <- read_csv("../data/processed/ukb_phenos_unrelated.csv", show_col_types = FALSE) %>%
  mutate(across(matches("^gPC"), ~ . * bmi, .names = "bmiBy{.col}"))

bmi_pgs_df <- read_table("../data/processed/pgs/bmi_prset_hallmark.all_score", show_col_types = FALSE) %>%
  select(id = IID, bmi_pgs = `Base_0.001`)

exposures <- c("bmi")
exposures_clean <- c("BMI")

outcomes <- c("alt_log")
outcomes_clean <- c("log(ALT)")

pathway_groups <- c("kegg_legacy")
pathway_groups_clean <- c("KEGG Legacy")

thresholds <- c("0.001")
threshold_names <- paste0("Pt_", c("0.001"))
threshold_names_clean <- c("P<0.001")

adjustments <- c("bmi_pgs_main", "bmi_pgs_by_bmi")

ppgs_config_df <- expand_grid(
  e = exposures,
  y = outcomes,
  pathway_group = pathway_groups,
  threshold = thresholds,
  adjustment = adjustments
)

gPCs <- paste0("gPC", 1:10)
primary_covars <- c("sex", "age", "age_squared", "ageBySex", gPCs)
liver_bms <- c("alt_log")

read_ppgs <- function(tag, pg, threshold) {
  ppgs_df <- read_delim(paste0("../data/processed/pgs/", tag, "_prset_", pg, ".all_score"),
                        delim = " ", show_col_types = FALSE) %>%
    select(id = IID, contains(paste0("_", threshold))) %>%
    rename_with(~ gsub("_[0-9\\.e]+$", "", .x))
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

test_all_ppgs_gxe <- function(e, y, pg, threshold, covars, adjustment,
                              standardize = TRUE) {
  print(paste0("Testing ", e, " - ", y, " in the held-out testing set using ", 
               pg, " pathways and ", threshold, " threshold, with ", adjustment, " adjustment."))
  ppgs_df <- read_ppgs(y, pg, threshold)
  pathways <- setdiff(names(ppgs_df), "id")
  if (adjustment == "bmi_pgs_main") {
    covars <- c(covars, "bmi_pgs")
  } else if (adjustment == "bmi_pgs_bmi_int") {
    covars <- c(covars, "bmi_pgs * bmi")
  }
  regression_df <- ukb_phenos %>%
    inner_join(ppgs_df, by = "id") %>%
    inner_join(bmi_pgs_df, by = "id")
  mclapply(pathways, function(p) {
    regression_df$ppgs <- regression_df[[p]]
    test_single_ppgs_gxe(e, y, regression_df, covars, standardize)
  }, mc.cores = 8) %>%
    setNames(pathways) %>%
    bind_rows(.id = "pathway")
}

ppgs_res_df <- ppgs_config_df %>%
  select(e, y, pathway_group, threshold, adjustment) %>%
  rowwise() %>%
  mutate(pathway_fits = list(
    test_all_ppgs_gxe(e, y, pathway_group, threshold,
                      c(primary_covars, paste0(e, "By", gPCs)),
                      adjustment)
  )) %>%
  unnest(pathway_fits)

saveRDS(ppgs_res_df, "../data/processed/pgs/ppgs_res_df_bmiPGS.rds")
