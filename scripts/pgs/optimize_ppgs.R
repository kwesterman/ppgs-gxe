library(tidyverse)


# args <- commandArgs(trailingOnly = TRUE)
# bm <- args[1]


# ukb_subset_names <- c("training", "validation", "testing")
# ukb_subset_names_clean <- c("Training", "Validation", "Testing")
# ukb_subsets <- map(set_names(ukb_subset_names), function(set) {
#   read_csv(paste0("../data/processed/ukb_", set, "_set.csv"), show_col_types = FALSE) %>%
#     mutate(across(matches("^gPC"), ~ . * bmi, .names = "bmiBy{.col}"))
# })
ukb_validation_df <- read_csv(paste0("../data/processed/ukb_validation_set.csv"), 
                              show_col_types = FALSE) %>%
  mutate(across(matches("^gPC"), ~ . * bmi, .names = "bmiBy{.col}"))


exposures <- c("bmi")
exposures_clean <- c("BMI")

# outcomes <- c("alb", "alt_log")
outcomes <- c("hscrp_log", "hba1c", "ldl_statinadj", "tg_log", "alt_log",
              "ggt_log", "hdl", "alb", "lipA_log", "vitD")

# threshold_vec <- 5 * 10 ** seq(-8, -2)
thresholds <- c("5e-08", "0.005")
# threshold_names <- paste0("Pt_", c("5e-08", "0.005"))
# threshold_names_clean <- c("P<5e-8", "P<0.005")

gPCs <- paste0("gPC", 1:10)
primary_covars <- c("sex", "age", "age_squared", "ageBySex", gPCs)

pgs_config_df <- expand_grid(
  e = exposures,
  y = outcomes,
  threshold = c(thresholds, "all")
) %>%
  mutate(tag = paste0(y, "_prset"))


read_pgs_prset <- function(tag, thresh) {
  ppgs_df <- read_delim(paste0("../data/processed/pgs/", tag, ".all_score"), 
                        delim = " ", show_col_types = FALSE) %>%
    select(id = IID, everything(), -FID, -contains("Base"))
  ppgs_df
}

get_ppgs_weights <- function(e, y, tag, thresh, covars) {
  pgs_df <- read_pgs_prset(tag, thresh)
  if (!(thresh == "all")) {
    pgs_df <- select(pgs_df, id, contains(thresh))
  }
  pathways <- setdiff(names(pgs_df), "id")
  regression_df <- inner_join(ukb_validation_df, pgs_df, by = "id")
  lm_form_str <- paste0(y, " ~ pgs * ", e)
  if (!is.null(covars)) {
    lm_form_str <- paste0(lm_form_str, " + ", paste(covars, collapse = " + "))
  }
  univariate_ppgs_gxe_df <- parallel::mclapply(
    setNames(pathways, pathways), 
    function(p) {
      regression_df$pgs <- regression_df[[p]]
      lm_fit <- lm(as.formula(lm_form_str), data = regression_df)
      lm_fit %>%
        broom::tidy() %>%
        filter(term == paste0("pgs:", e))
    }, 
    mc.cores = 1) %>%
    bind_rows(.id = "pathway")
  univariate_ppgs_gxe_df
}

calc_ppgs <- function(tag, weights_df) {
  pgs_df <- read_pgs_prset(tag)
  pgs_df$ppgs <- as.vector(as.matrix(pgs_df[weights_df$pathway]) %*% weights_df$estimate)
  select(pgs_df, id, ppgs)
}

ppgs_weights_df <- pgs_config_df %>%
  select(e, y, tag, threshold) %>%
  distinct() %>% 
  rowwise() %>%
  mutate(pathway_fits = list(
    get_ppgs_weights(e, y, paste0(y, "_prset"), threshold,
                     c(primary_covars, paste0(e, "By", gPCs)))
  ))

ppgs_vec_df <- ppgs_weights_df %>%
  rowwise() %>%
  mutate(ppgs_df = list(calc_ppgs(tag, pathway_fits)))

saveRDS(ppgs_vec_df, "../data/processed/pgs/ppgs_vec_df.rds")
