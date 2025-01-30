library(tidyverse)
library(parallel)


ukb_subset_names <- c("training", "testing")
ukb_subset_names_clean <- c("Training", "Testing")
ukb_subsets <- map(set_names(ukb_subset_names), function(set) {
  read_csv(paste0("../data/processed/ukb_", set, "_set.csv"), show_col_types = FALSE) %>%
    mutate(across(matches("^gPC"), ~ . * bmi, .names = "bmiBy{.col}"))
})

exposures <- c("bmi")
exposures_clean <- c("BMI")

outcomes <- c("alt_log", "ast_log", "ggt_log")
outcomes_clean <- c("log(ALT)", "log(AST)", "log(GGT)")

pathway_groups <- c("hallmark", "kegg")
pathway_groups_clean <- c("Hallmark", "KEGG")

thresholds <- c("5e-08", "0.001")
threshold_names <- paste0("Pt_", c("5e-08", "0.001"))
threshold_names_clean <- c("P<5e-8", "P<0.001")

ppgs_config_df <- expand_grid(
  e = exposures,
  y = outcomes,
  pathway_group = pathway_groups,
  threshold = thresholds
)

gPCs <- paste0("gPC", 1:10)
primary_covars <- c("sex", "age", "age_squared", "ageBySex", gPCs)
liver_bms <- c("alt_log", "ast_log", "ggt_log")

read_pgs_prset <- function(tag, pg) {
  ppgs_df <- read_delim(paste0("../data/processed/pgs/", tag, "_prset_", pg, ".all_score"),
                        delim = " ", show_col_types = FALSE) %>%
    select(id = IID, everything(), -FID)
  # rename_with(~ str_replace_all(., "_[\\d\\.]+$", ""), -id)
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

test_all_ppgs_gxe <- function(e, y, pg, threshold, set, covars, 
                              ancestry_test = "all", standardize = TRUE) {
  print(paste0("Testing ", e, " - ", y, " using ", pg, " pathways and ", threshold, " threshold."))
  ppgs_df <- read_pgs_prset(y, pg) %>%
    select(id, contains(threshold))
  pathways <- gsub(paste0("_", threshold), "", setdiff(names(ppgs_df), "id"))
  regression_df <- inner_join(ukb_subsets[[set]], ppgs_df, by = "id")
  if (ancestry_test != "all") {
    regression_df <- inner_join(regression_df, ancestry_df, by = "id") %>%
      filter(ancestry == ancestry_test)
  }
  mclapply(pathways, function(p) {
    regression_df$ppgs <- regression_df[[paste0(p, "_", threshold)]]
    test_single_ppgs_gxe(e, y, regression_df, covars, standardize)
  }, mc.cores = 8) %>%
    setNames(pathways) %>%
    bind_rows(.id = "pathway")
}

ppgs_res_df <- ppgs_config_df %>%
  select(e, y, pathway_group, threshold) %>%
  rowwise() %>%
  mutate(pathway_fits = list(
    test_all_ppgs_gxe(e, y, pathway_group, threshold, "testing",
                      c(primary_covars, paste0(e, "By", gPCs)))
  )) %>%
  unnest(pathway_fits)

saveRDS(ppgs_res_df, "../data/processed/pgs/ppgs_res_df.rds")
