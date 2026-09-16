# =============================================================================
# External validation of the Paediatric Chiang Mai University Extubation
# Failure Predictive Score (Ped-CMU ExFPS)
#
# Companion analysis script for the manuscript:
#   "External Validation of the Ped-CMU Extubation Failure Predictive Score and
#    Clinician Adherence in Pediatric Cardiac Patients"
#
# This script reproduces, top to bottom:
#   1. Data import and harmonisation of the development and external-validation
#      datasets into a single pooled ("membership") dataset.
#   2. Table 1 - baseline characteristics, development vs. validation.
#   3. Table 2 - distribution of the score and its linear predictor.
#   4. Table 3 - predictor-outcome associations in each dataset.
#   5. Case-mix dissimilarity (membership model).
#   6. Discrimination - AUROC, with single-cohort and dev-vs-val ROC curves.
#   7. Calibration - calibration curve, slope and calibration-in-the-large.
#   8. Table 4 - diagnostic indices at the pre-defined cut-off (score >= 5).
#   9. Decision curve analysis.
#  10. Exploratory analysis of clinician adherence to score-based management.
#
# Score definition (linear predictor on the logit scale):
#   score = 10 * history_of_reintubation
#         +  4 * pneumonia
#         +  1 * acyanosis
#         +  6 * cyanosis (SpO2 > 85%)
#   linear predictor = -3.286 + 0.231 * score
#
# Reproducibility
#   Tested with R 4.5.2. Set the two paths in the "Paths" block below, then run
#   the script top to bottom. All packages are available on CRAN:
#     install.packages(c(
#       "readxl", "dplyr", "tidyr", "tibble", "forcats", "stringr",
#       "janitor", "haven", "ggplot2", "gtsummary", "gt", "broom",
#       "pROC", "epiR", "binom"))
# =============================================================================


# ---- 0. Paths and packages --------------------------------------------------

# Edit these two lines to point at the folder holding the (de-identified) data
# and the desired output location, then run the script top to bottom.
# The de-identified datasets are available from the corresponding author on
# reasonable request.
data_dir   <- "data"      # folder containing the .xlsx datasets
output_dir <- "output"    # folder for tables, figures and exported data

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Data manipulation
library(readxl)    # read_xlsx
library(dplyr)     # data wrangling verbs
library(tidyr)     # pivot_*, complete
library(tibble)    # tibble
library(forcats)   # fct_relevel
library(stringr)   # str_detect
library(janitor)   # clean_names
library(haven)     # write_dta
# Visualisation
library(ggplot2)
# Tables
library(gtsummary) # tbl_summary, tbl_regression, add_p, add_n
library(gt)        # gtsave
# Modelling and validation
library(broom)     # tidy
library(pROC)      # roc, ggroc, coords
library(epiR)      # epi.tests
library(binom)     # binom.confint

# Helper: format p-values consistently (< 0.001, otherwise 3 decimals).
format_pvalue <- function(p) {
  ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
}


# ---- 1. Import and harmonise the two datasets -------------------------------

df     <- readxl::read_xlsx(file.path(data_dir, "external_validation_dataset.xlsx"))
df_dev <- readxl::read_xlsx(file.path(data_dir, "development_dataset.xlsx"), sheet = "Sheet1")

# One validation record (record_id 71, an infant aged 5.7 months) has body
# weight and height entered in each other's fields (58 kg, 3.8 cm). Swap them
# back to 3.8 kg and 58 cm and recompute body-mass index with height in metres.
df <- df |>
  dplyr::mutate(
    bw_raw = bw,
    bw     = dplyr::if_else(record_id == 71, ht, bw),
    ht     = dplyr::if_else(record_id == 71, bw_raw, ht),
    bmi    = dplyr::if_else(record_id == 71, bw / (ht / 100)^2, as.numeric(bmi)),
    bw_raw = NULL
  )

# Direct identifiers are dropped. These are no-ops on the de-identified data and
# simply document the de-identification step.
df$record_id <- NULL
df$hn        <- NULL
df$name      <- NULL

# External-validation dataset: type coercion and resolution of duplicated
# cyanosis columns produced on export (...22 / ...34).
df <- df |>
  dplyr::mutate(
    bmi        = as.numeric(bmi),
    age_extumo = as.numeric(age_extumo),
    cyanosis_icu_y...22 = NULL
  ) |>
  dplyr::rename(cyanosis_icu = cyanosis_icu_y...34)

# Development dataset: derive the two cyanosis indicators and align column names
# with the external-validation dataset.
df_dev <- df_dev |>
  dplyr::mutate(
    `Acyanosis score` = dplyr::if_else(acycyless85cy85_2recode == 2, 1, 0),
    `Cyanosis score`  = dplyr::if_else(acycyless85cy85_2recode == 3, 1, 0)
  ) |>
  dplyr::rename(
    "Event of extubation"   = "N",
    "age_extumo"            = "Age_extumo",
    "sex"                   = "Male",
    "bw"                    = "BW",
    "ht"                    = "Ht",
    "bmi"                   = "BMI",
    "genetic_y"             = "Genetic_Y",
    "genetic_c"             = "Genetic_C",
    "cardiac_c"             = "Cardiac_C",
    "postop_y"              = "Postop_Y",
    "pah"                   = "PAH",
    "dur_intu"              = "Dur_intu",
    "Pneumonia score"       = "pneumonia_score",
    "Reintubation score"    = "reintubation_score",
    "Mx 2 steriod"          = "steriodpre",
    "e_sat"                 = "Sat",
    "icu_stay"              = "ICUstay",
    "death"                 = "Death",
    "mode_v"                = "Mode_V",
    "cyanotic_native_y"     = "Cyatotic_native_YN",
    "cyanosis_icu"          = "Cyanosis_ICU_YN",
    "outcome_of_extubation" = "extubationfailure",
    "hx_reintu"             = "reintubation_Y",
    "pneumonia"             = "Pneumonia"
  )

# A small number of validation records carry a cyanosis score of 6; recode to 1.
df <- df |>
  dplyr::mutate(
    `Cyanosis score` = dplyr::if_else(`Cyanosis score` == 6, 1, `Cyanosis score`)
  )


# ---- 2. Pooled ("membership") dataset and score -----------------------------

# Stack the two datasets on their shared columns; dataset = 0 (development),
# 1 (external validation).
common_name <- intersect(names(df), names(df_dev))
df_member <- dplyr::bind_rows(
  list(`0` = df_dev[, common_name], `1` = df[, common_name]),
  .id = "dataset"
) |>
  dplyr::mutate(dataset = as.numeric(dataset))

# Single cyanosis indicator: 0 = none, 1 = acyanosis, 2 = cyanosis (SpO2 > 85%).
# History of reintubation is collapsed to a 0/1 indicator.
df_member <- df_member |>
  dplyr::mutate(
    cyanosis = dplyr::case_when(
      `Acyanosis score` == 1 ~ 1,
      `Cyanosis score`  == 1 ~ 2,
      TRUE ~ 0
    ),
    hx_reintu = dplyr::if_else(hx_reintu == 2, 1, hx_reintu)
  )

# Ped-CMU ExFPS score and its linear predictor. cyanosis is coded 0/1/2, so
# acyanosis is cyanosis == 1 (+1 point) and cyanosis with SpO2 > 85% is
# cyanosis == 2 (+6 points).
df_member <- df_member |>
  dplyr::mutate(
    score    = (10 * hx_reintu) + (4 * pneumonia) +
               (1 * (cyanosis == 1)) + (6 * (cyanosis == 2)),
    lp_score = -3.286 + (0.231 * score)
  )


# ---- 3. Table 1: baseline characteristics -----------------------------------

table1 <- df_member |>
  dplyr::select(
    dataset, age_extumo, sex, bw, ht, bmi, genetic_y, genetic_c, cardiac_c,
    cyanotic_native_y, cyanosis_icu, e_sat, pah, dur_intu, postop_y,
    pneumonia, hx_reintu, `Mx 2 steriod`, balance_kg, outcome_of_extubation, death
  ) |>
  dplyr::mutate(
    genetic_c = factor(
      genetic_c,
      levels = 1:5,
      labels = c("Down's syndrome", "Heterotaxy syndrome", "Digeorge's syndrome",
                 "VacTERL syndrome", "Other")
    ),
    cardiac_c = factor(
      cardiac_c,
      levels = 1:18,
      labels = c("Single ventricle", "VSD", "PDA", "COA", "ASD",
                 "Hemitruncus / truncus", "IAA", "TOF", "DORV", "TGA", "ccTGA",
                 "AVSD", "Abnormal coro / ALCAPA", "PS", "Ebstein",
                 "Hypertropic", "TAPVR", "Other")
    )
  ) |>
  dplyr::rename(
    "Age (month)"                                = "age_extumo",
    "Male"                                       = "sex",
    "Body weight (kg)"                           = "bw",
    "Height (cm)"                                = "ht",
    "Body mass index (kg/m2)"                    = "bmi",
    "Having genetic disease"                     = "genetic_y",
    "Genetic disease"                            = "genetic_c",
    "Cardiac disease"                            = "cardiac_c",
    "Native anatomy cyanosis"                    = "cyanotic_native_y",
    "Status ICU lesion cyanosis"                 = "cyanosis_icu",
    "Oxygen saturation"                          = "e_sat",
    "Paediatric pulmonary arterial hypertension" = "pah",
    "Duration of intubation"                     = "dur_intu",
    "Post-operation"                             = "postop_y",
    "Pneumonia"                                  = "pneumonia",
    "History of intubation"                      = "hx_reintu",
    "Treated with steroid"                       = "Mx 2 steriod",
    "Balance per kg"                             = "balance_kg",
    "Extubation failure"                         = "outcome_of_extubation"
  ) |>
  dplyr::mutate(
    `Extubation failure` = factor(`Extubation failure`, levels = c(1, 0),
                                  labels = c("Fail", "Success")),
    dataset = factor(dataset, levels = c(1, 0),
                     labels = c("External validation dataset", "Development dataset"))
  ) |>
  gtsummary::tbl_summary(
    by = dataset,
    missing = "no",
    statistic = list(all_continuous() ~ "{mean} (± {sd})")
  ) |>
  gtsummary::add_n(
    statistic = "{N_miss} ({p_miss}%)",
    col_label = "**Missing**",
    last = FALSE
  ) |>
  gtsummary::add_p(
    test = list(
      all_continuous()  ~ "t.test",
      all_categorical() ~ "fisher.test"
    ),
    pvalue_fun = ~ gtsummary::style_pvalue(.x, digits = 3)
  )

table1
gtsave(gtsummary::as_gt(table1), filename = file.path(output_dir, "table1_baseline.docx"))


# ---- 4. Table 2: score and linear-predictor distribution --------------------

table2 <- df_member |>
  dplyr::select(dataset, score, lp_score) |>
  dplyr::mutate(
    dataset = factor(dataset, levels = c(1, 0),
                     labels = c("External validation dataset", "Development dataset"))
  ) |>
  dplyr::rename(
    "Linear predictor of Ped-CMU ExFPS" = "lp_score",
    "Ped-CMU ExFPS score"               = "score"
  ) |>
  gtsummary::tbl_summary(
    by = dataset,
    missing = "no",
    type = list(c("Linear predictor of Ped-CMU ExFPS",
                  "Ped-CMU ExFPS score") ~ "continuous"),
    statistic = list(c("Linear predictor of Ped-CMU ExFPS",
                       "Ped-CMU ExFPS score") ~ "{mean} (± {sd})")
  ) |>
  gtsummary::add_p(
    test = list(c("Linear predictor of Ped-CMU ExFPS",
                  "Ped-CMU ExFPS score") ~ "t.test"),
    pvalue_fun = ~ gtsummary::style_pvalue(.x, digits = 3)
  )

table2
gtsave(gtsummary::as_gt(table2), filename = file.path(output_dir, "table2_score_distribution.docx"))


# ---- 5. Table 3: predictor-outcome associations, dev vs. val ----------------

dev_model <- glm(
  outcome_of_extubation ~ hx_reintu + pneumonia + factor(cyanosis),
  family = binomial(link = "logit"),
  data = dplyr::filter(df_member, dataset == 0)
)
val_model <- glm(
  outcome_of_extubation ~ hx_reintu + pneumonia + factor(cyanosis),
  family = binomial(link = "logit"),
  data = dplyr::filter(df_member, dataset == 1)
)

broom::tidy(dev_model, exponentiate = TRUE, conf.int = TRUE)
broom::tidy(val_model, exponentiate = TRUE, conf.int = TRUE)

# Development model: add a beta (log-OR) column alongside the odds ratio.
table3_dev <- dev_model |>
  gtsummary::tbl_regression(
    exponentiate = TRUE,
    pvalue_fun = gtsummary::label_style_pvalue(digits = 3)
  ) |>
  gtsummary::modify_column_unhide(estimate) |>
  gtsummary::modify_table_body(
    ~ .x |> dplyr::mutate(coef = round(log(estimate), 3), .before = estimate)
  ) |>
  gtsummary::modify_header(coef = "**Beta**")

# Validation model: the cyanosis estimates are not identifiable (zero events in
# one category), so they are blanked out.
table3_val <- val_model |>
  gtsummary::tbl_regression(
    exponentiate = TRUE,
    pvalue_fun = gtsummary::label_style_pvalue(digits = 3)
  ) |>
  gtsummary::modify_column_unhide(estimate) |>
  gtsummary::modify_table_body(
    ~ .x |>
      dplyr::mutate(coef = round(log(estimate), 3), .before = estimate) |>
      dplyr::mutate(
        across(
          c(coef, estimate, p.value),
          ~ dplyr::if_else(stringr::str_detect(variable, "cyanosis"), NA_real_, .x)
        ),
        ci = dplyr::if_else(stringr::str_detect(variable, "cyanosis"), NA_character_, ci)
      )
  ) |>
  gtsummary::modify_header(coef = "**Beta**")

table3 <- gtsummary::tbl_merge(
  tbls = list(table3_dev, table3_val),
  tab_spanner = c("Development model", "Validation model")
) |>
  gtsummary::modify_table_styling(
    columns = label,
    rows = stringr::str_detect(label, "cyanosis"),
    footnote = paste(
      "Cyanosis estimates could not be calculated in the validation model",
      "due to zero events in one category"
    )
  )

table3
gtsave(gtsummary::as_gt(table3), filename = file.path(output_dir, "table3_predictor_outcome.docx"))


# ---- 6. Case-mix dissimilarity: membership model ----------------------------

# Linear-predictor distribution in each dataset (overlaid histograms + fitted
# normal densities).
histo_stat <- df_member |>
  dplyr::group_by(dataset) |>
  dplyr::summarise(
    mean_lp = mean(lp_score, na.rm = TRUE),
    sd_lp   = sd(lp_score, na.rm = TRUE),
    .groups = "drop"
  )
histo_dev <- dplyr::filter(histo_stat, dataset == 0)
histo_val <- dplyr::filter(histo_stat, dataset == 1)

ggplot(df_member, aes(x = lp_score, fill = factor(dataset), group = factor(dataset))) +
  geom_histogram(aes(y = after_stat(density)), alpha = 0.5, position = "identity") +
  geom_function(fun = dnorm, args = list(mean = histo_dev$mean_lp, sd = histo_dev$sd_lp),
                colour = "red") +
  geom_function(fun = dnorm, args = list(mean = histo_val$mean_lp, sd = histo_val$sd_lp),
                colour = "darkgreen") +
  labs(x = "Linear predictor", y = "Density", fill = "Dataset") +
  theme_classic(base_size = 12)

# Membership model: how well the predictors and the outcome distinguish the two
# datasets (Debray et al. 2015). The response is dataset membership, so the ROC
# is scored against dataset, not against extubation failure.
member_model <- glm(
  dataset ~ hx_reintu + pneumonia + factor(cyanosis) + outcome_of_extubation,
  family = binomial(link = "logit"),
  data = df_member
)
broom::tidy(member_model, conf.int = TRUE)

case_mix_roc <- pROC::roc(
  response  = df_member$dataset,
  predictor = predict(member_model, type = "response"),
  ci = TRUE
)


# ---- 7. Discrimination: AUROC -----------------------------------------------

# Validation cohort only.
dfm <- dplyr::filter(df_member, dataset == 1)

val_lp_model <- glm(
  outcome_of_extubation ~ lp_score,
  family = binomial(link = "logit"),
  data = dfm
)

roc <- pROC::roc(
  response  = dfm$outcome_of_extubation,
  predictor = predict(val_lp_model, type = "response"),
  ci = TRUE
)

roc_coords <- pROC::coords(
  roc, x = "all",
  ret = c("threshold", "sensitivity", "specificity", "ppv", "npv",
          "1-specificity", "1-sensitivity", "lr_pos", "lr_neg"),
  transpose = FALSE
)

roc_plot <- pROC::ggroc(roc, alpha = 1, colour = "#7c0000", legacy.axes = TRUE) +
  annotate("text", x = 0.6, y = 0.1,
           label = sprintf("AUROC = %.2f, 95%%CI %.2f-%.2f",
                           roc$ci[2], roc$ci[1], roc$ci[3])) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 1,
           colour = "grey", linetype = "dashed") +
  geom_point(aes(x = roc_coords$`1-specificity`, y = roc_coords$sensitivity),
             colour = "#7c0000", size = 2) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  theme_classic(base_size = 12) +
  theme(
    axis.line   = element_line(colour = "black", linewidth = 0.8),
    axis.ticks  = element_line(colour = "black", linewidth = 0.6),
    axis.text   = element_text(colour = "black", size = 11),
    axis.title  = element_text(size = 12, face = "bold"),
    plot.title  = element_text(hjust = 0.5, size = 14, face = "bold")
  )
roc_plot
ggsave(file.path(output_dir, "roc.pdf"), roc_plot)

cat(sprintf("AUROC: %.2f (%.2f - %.2f)\n", roc$ci[2], roc$ci[1], roc$ci[3]))

# Development vs. validation ROC curves on the same axes.
roc_dev <- pROC::roc(response = df_member$outcome_of_extubation[df_member$dataset == 0],
                     predictor = df_member$lp_score[df_member$dataset == 0], ci = TRUE)
roc_val <- pROC::roc(response = df_member$outcome_of_extubation[df_member$dataset == 1],
                     predictor = df_member$lp_score[df_member$dataset == 1], ci = TRUE)

coords_dev <- pROC::coords(roc_dev, x = "all",
  ret = c("sensitivity", "1-specificity"), transpose = FALSE)
coords_val <- pROC::coords(roc_val, x = "all",
  ret = c("sensitivity", "1-specificity"), transpose = FALSE)

label_dev <- sprintf("Development %.2f (%.2f-%.2f)",
                     roc_dev$auc, roc_dev$ci[1], roc_dev$ci[3])
label_val <- sprintf("Validation %.2f (%.2f-%.2f)",
                     roc_val$auc, roc_val$ci[1], roc_val$ci[3])

roc_compare_plot <- pROC::ggroc(list(Development = roc_dev, Validation = roc_val),
                                legacy.axes = TRUE) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour = "grey", linetype = "dashed") +
  geom_point(data = coords_dev, aes(x = `1-specificity`, y = sensitivity),
             colour = "red", size = 2, inherit.aes = FALSE) +
  geom_point(data = coords_val, aes(x = `1-specificity`, y = sensitivity),
             colour = "darkgreen", size = 2, inherit.aes = FALSE) +
  scale_color_manual(
    name = "Area under ROC curve",
    values = c(Development = "red", Validation = "darkgreen"),
    labels = c(Development = label_dev, Validation = label_val)
  ) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  theme_classic(base_size = 12) +
  theme(
    axis.line     = element_line(colour = "black", linewidth = 0.8),
    axis.ticks    = element_line(colour = "black", linewidth = 0.6),
    axis.text     = element_text(colour = "black", size = 11),
    axis.title    = element_text(size = 12, face = "bold"),
    plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = c(0.7, 0.15),
    legend.direction = "vertical",
    legend.text   = element_text(size = 10),
    legend.title  = element_text(size = 12, face = "bold")
  ) +
  guides(color = guide_legend(ncol = 1))
roc_compare_plot
ggsave(file.path(output_dir, "roc_eval.pdf"), roc_compare_plot)
# Figure 2 for the manuscript: 300 dpi LZW TIFF (ragg writes the dpi header).
ggsave(file.path(output_dir, "figure2_roc.tiff"), roc_compare_plot,
       device = ragg::agg_tiff, width = 939 / 220, height = 939 / 220,
       units = "in", dpi = 300, compression = "lzw")


# ---- 8. Calibration ---------------------------------------------------------

# Predicted probability for the validation cohort.
cal <- tibble::tibble(
  pred = plogis(dfm$lp_score),
  obs  = dfm$outcome_of_extubation,
  lp   = dfm$lp_score
) |>
  tidyr::drop_na()

# Group predictions into deciles for the grouped observed-vs-expected points.
cal$bin <- cut(
  cal$pred,
  breaks = unique(quantile(cal$pred, probs = seq(0, 1, by = 0.1))),
  include.lowest = TRUE, labels = FALSE
)

bin_stats <- cal |>
  dplyr::group_by(bin) |>
  dplyr::summarise(
    exp = mean(pred),
    obs = mean(obs),
    n   = dplyr::n(),
    se  = sqrt((obs * (1 - obs)) / n),
    lci = pmax(0, obs - 1.96 * se),
    uci = pmin(1, obs + 1.96 * se),
    .groups = "drop"
  )

cal <- dplyr::arrange(cal, pred)

# Calibration slope: outcome ~ linear predictor.
cal_model_slope <- glm(obs ~ lp, family = binomial(link = "logit"), data = cal)
cal_slope <- tibble::tibble(
  slope = coef(cal_model_slope)[2],
  lci   = confint(cal_model_slope)[2, 1],
  uci   = confint(cal_model_slope)[2, 2]
)

# Calibration-in-the-large: intercept-only model with the linear predictor as
# an offset.
cal_model_citl <- glm(obs ~ 1, offset = lp, family = binomial(link = "logit"), data = cal)
cal_citl <- tibble::tibble(
  citl = coef(cal_model_citl)[1],
  lci  = confint(cal_model_citl)[1],
  uci  = confint(cal_model_citl)[2]
)

# Logistic calibration curve with a confidence band.
# The band is computed on the logit scale and back-transformed (Jensen's
# inequality: averaging on the probability scale would be biased).
cal_curve <- tibble::tibble(
  x = cal$pred,
  y = predict(cal_model_slope, newdata = cal, type = "response")
)
link_pred <- predict(cal_model_slope, newdata = cal, type = "link", se.fit = TRUE)
cal_curve$lci <- plogis(link_pred$fit - 1.96 * link_pred$se.fit)
cal_curve$uci <- plogis(link_pred$fit + 1.96 * link_pred$se.fit)

calibration_plot <- ggplot() +
  # Perfect-calibration reference line.
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1, linetype = "Perfect calibration"),
               colour = "gray40", linewidth = 0.5) +
  # Calibration curve and its 95% CI band.
  geom_ribbon(data = cal_curve, aes(x = x, ymin = lci, ymax = uci,
              fill = "Calibration curve 95% CI"), alpha = 0.2) +
  geom_line(data = cal_curve, aes(x = x, y = y, colour = "Calibration curve"),
            linewidth = 1) +
  # Grouped (decile) observed-vs-expected points and CIs.
  geom_point(data = bin_stats, aes(x = exp, y = obs, colour = "Grouped observations"),
             size = 2.5) +
  geom_errorbar(data = bin_stats, aes(x = exp, ymin = lci, ymax = uci,
                colour = "Grouped observations"), width = 0.02) +
  # Slope and CITL annotation, right-aligned so it stays inside the panel.
  annotate("text", x = 1, y = 0.12, hjust = 1, size = 3.2,
           label = sprintf("Slope: %.2f (%.2f, %.2f)\nCITL: %.2f (%.2f, %.2f)",
                           cal_slope$slope, cal_slope$lci, cal_slope$uci,
                           cal_citl$citl, cal_citl$lci, cal_citl$uci)) +
  # Rug plots: events (red) above, non-events (blue) below.
  geom_segment(data = dplyr::filter(cal, obs == 0), aes(x = pred, xend = pred),
               y = -0.08, yend = -0.05, colour = "blue", linewidth = 0.3) +
  geom_segment(data = dplyr::filter(cal, obs == 1), aes(x = pred, xend = pred),
               y = -0.05, yend = -0.02, colour = "red", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = 1, y = -0.05, yend = -0.05), colour = "black") +
  scale_colour_manual(name = NULL,
    values = c("Calibration curve" = "#5bb0ff", "Grouped observations" = "forestgreen")) +
  scale_fill_manual(name = NULL, values = c("Calibration curve 95% CI" = "#5bb0ff")) +
  scale_linetype_manual(name = NULL, values = c("Perfect calibration" = "dashed")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
  labs(x = "Predicted probability", y = "Observed proportion", title = "Calibration plot") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(colour = "black", linewidth = 0.5),
    axis.ticks        = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text         = element_text(colour = "black"),
    plot.title        = element_text(hjust = 0.5, face = "bold"),
    aspect.ratio      = 1,
    # Legend to the right: the figure is wider than it is tall.
    legend.position   = "right"
  )
calibration_plot
ggsave(file.path(output_dir, "calibration_plot.pdf"), calibration_plot)
# Figure 3 for the manuscript: 300 dpi LZW TIFF.
ggsave(file.path(output_dir, "figure3_calibration.tiff"), calibration_plot,
       device = ragg::agg_tiff, width = 1430 / 220, height = 884 / 220,
       units = "in", dpi = 300, compression = "lzw")


# ---- 9. Table 4: diagnostic indices at the cut-off (score >= 5) -------------

# Use the score from section 2 and dichotomise at the pre-specified high-risk
# cut-off of 5.
df_diag <- df_member |>
  dplyr::mutate(
    disease = factor(outcome_of_extubation, levels = c(1, 0)),
    test    = factor(dplyr::if_else(score >= 5, 1, 0), levels = c(1, 0))
  ) |>
  dplyr::select(disease, test, dataset)

tab_dev <- with(dplyr::filter(df_diag, dataset == 0), table(test, disease))
tab_val <- with(dplyr::filter(df_diag, dataset == 1), table(test, disease))

diag_dev <- epiR::epi.tests(tab_dev, method = "exact", digits = 2, conf.level = 0.95)
diag_val <- epiR::epi.tests(tab_val, method = "exact", digits = 2, conf.level = 0.95)

keep_stats <- c("tp", "se", "sp", "diag.ac", "pv.pos", "pv.neg", "lr.pos", "lr.neg")
stat_labels <- c(
  tp      = "True prevalence",
  se      = "Sensitivity",
  sp      = "Specificity",
  diag.ac = "Diagnosis accuracy",
  pv.pos  = "Positive predictive value",
  pv.neg  = "Negative predictive value",
  lr.pos  = "Likelihood ratio of positive",
  lr.neg  = "Likelihood ratio of negative"
)

table4 <- dplyr::bind_rows(
  dplyr::mutate(dplyr::filter(diag_dev$detail, statistic %in% keep_stats), Dataset = "Development"),
  dplyr::mutate(dplyr::filter(diag_val$detail, statistic %in% keep_stats), Dataset = "Validation")
) |>
  dplyr::mutate(
    statistic = dplyr::recode(statistic, !!!stat_labels),
    estimate  = sprintf("%.2f%% (%.2f, %.2f)", est * 100, lower * 100, upper * 100)
  ) |>
  dplyr::select(Dataset, statistic, estimate) |>
  tidyr::pivot_wider(names_from = Dataset, values_from = estimate) |>
  gt::gt()

table4
gtsave(table4, filename = file.path(output_dir, "table4_diagnostic_indices.docx"))

# Export the pooled, cleaned dataset (Stata format) for downstream use.
haven::write_dta(janitor::clean_names(df_member), file.path(output_dir, "data.dta"))


# ---- 10. Decision curve analysis --------------------------------------------

dca_data <- tibble::tibble(
  y    = as.numeric(dfm$outcome_of_extubation),
  pred = plogis(dfm$lp_score)
) |>
  tidyr::drop_na()

n_total    <- nrow(dca_data)
prevalence <- mean(dca_data$y)
thresholds <- seq(0.01, 1, by = 0.001)

# Net benefit at each threshold for treat-all, treat-none and the model.
dca_result <- do.call(dplyr::bind_rows, lapply(thresholds, function(t) {
  w        <- t / (1 - t)
  flagged  <- dca_data$pred >= t
  tp_rate  <- sum(flagged & dca_data$y == 1) / n_total
  fp_rate  <- sum(flagged & dca_data$y == 0) / n_total
  tibble::tibble(
    threshold = t,
    treat_all  = prevalence - (1 - prevalence) * w,
    treat_none = 0,
    model      = tp_rate - fp_rate * w
  )
})) |>
  tidyr::pivot_longer(
    cols = c(treat_all, treat_none, model),
    names_to = "strategy", values_to = "net_benefit"
  ) |>
  dplyr::mutate(strategy = factor(strategy, levels = c("treat_none", "treat_all", "model")))

dca_plot <- ggplot(dca_result, aes(x = threshold, y = net_benefit,
                                    group = strategy, colour = strategy)) +
  geom_line(linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 0.5)) +
  scale_y_continuous(limits = c(-0.06, 0.075)) +
  scale_colour_manual(
    name = NULL,
    values = c(treat_all = "red", treat_none = "blue", model = "darkgreen"),
    labels = c(treat_all = "Treat all", treat_none = "Treat none", model = "Model")
  ) +
  labs(x = "Threshold probability", y = "Net benefit", title = "Decision curve analysis") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    aspect.ratio       = 1,
    axis.line          = element_line(colour = "black", linewidth = 0.5),
    axis.title.x       = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y       = element_text(size = 16, margin = margin(r = 10)),
    axis.text          = element_text(size = 12),
    plot.title         = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position    = "bottom",
    legend.text        = element_text(size = 14)
  )
dca_plot


# ---- 11. Exploratory: adherence to score-based management -------------------
# These analyses are exploratory and use the external-validation cohort only.

# Per-indication management columns and the derived score / risk group.
df_mgmt <- df |>
  dplyr::select(
    hx_reintu, pneumonia, `Acyanosis score`, `Cyanosis score`,
    outcome_of_extubation, dplyr::starts_with("Mx"), `M3 negtiveIO`, balance_kg
  ) |>
  dplyr::mutate(
    cyanosis       = dplyr::case_when(
      `Acyanosis score` == 1 ~ 1,
      `Cyanosis score`  == 1 ~ 2,
      TRUE ~ 0
    ),
    `M3 negtiveIO` = dplyr::if_else(balance_kg < 0, 1, 0),
    score          = (10 * hx_reintu) + (4 * pneumonia) +
                     (1 * (cyanosis == 1)) + (6 * (cyanosis == 2)),
    risk_group     = factor(dplyr::if_else(score >= 5, 1, 0),
                            levels = c(0, 1), labels = c("Low risk", "High risk"))
  )

# 11a. Proportion receiving each management, by risk group.
mgmt_cols <- c("Mx1 Positive pressure", "Mx 2 steriod", "Mx3 Positive pressure",
               "Mx 3 decrease Fio2", "M3 negtiveIO")

mgmt_proportion <- df_mgmt |>
  dplyr::group_by(risk_group) |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(mgmt_cols),
      list(
        n       = ~ sum(!is.na(.x)),
        success = ~ sum(.x, na.rm = TRUE),
        mean    = ~ mean(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = -risk_group,
    names_to = c("category", ".value"),
    names_pattern = "(.+)_(n|success|mean)"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    ci  = list(binom::binom.confint(success, n, methods = "wilson")),
    lci = ci$lower,
    uci = ci$upper
  ) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    category = dplyr::recode(category,
      "Mx1 Positive pressure" = "NIPPV assistance to prevent pneumonia",
      "Mx 2 steriod"          = "Systemic steroid",
      "Mx3 Positive pressure" = "NIPPV assistance to prevent cyanosis",
      "Mx 3 decrease Fio2"    = "Decreased FiO2 as possible",
      "M3 negtiveIO"          = "Controlled negative I/O to prevent imbalance"
    ),
    category = factor(category, levels = c(
      "Decreased FiO2 as possible",
      "Controlled negative I/O to prevent imbalance",
      "Systemic steroid",
      "NIPPV assistance to prevent cyanosis",
      "NIPPV assistance to prevent pneumonia"
    ))
  )

mgmt_proportion_plot <- ggplot(mgmt_proportion,
    aes(x = category, y = mean, fill = risk_group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(
    name = "Risk group",
    values = c("Low risk" = "#92cced", "High risk" = "#ea801c")
  ) +
  labs(x = "Management", y = "Proportion") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line   = element_line(),
    axis.ticks  = element_line(),
    legend.position = "bottom",
    axis.title.x = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(vjust = 0.5, face = "bold", size = 12, margin = margin(r = 10)),
    axis.text    = element_text(size = 8.5)
  )
mgmt_proportion_plot

# Protocol-adherence indicators (overall and per indication). A patient is
# "adherent" when treated if indicated and untreated if not.
df_protocol <- df_mgmt |>
  dplyr::mutate(
    `Mx 3 decrease Fio2` = dplyr::if_else(is.na(`Mx 3 decrease Fio2`), 1, `Mx 3 decrease Fio2`),
    protocol_reintu = dplyr::if_else(
      (hx_reintu == 1 & `Mx 2 steriod` == 1) | (hx_reintu != 1 & `Mx 2 steriod` != 1), 1, 0
    ),
    protocol_pneumonia = dplyr::if_else(
      (pneumonia == 1 & `Mx1 Positive pressure` == 1) |
        (pneumonia != 1 & `Mx1 Positive pressure` != 1), 1, 0
    ),
    protocol_cyanosis = dplyr::if_else(
      (cyanosis == 2 & `Mx3 Positive pressure` == 1 & `M3 negtiveIO` == 1) |
        (cyanosis != 2 & `Mx3 Positive pressure` != 1 & `M3 negtiveIO` != 1), 1, 0
    ),
    protocol = dplyr::if_else(
      protocol_reintu == 1 & protocol_pneumonia == 1 & protocol_cyanosis == 1, 1, 0
    ),
    get_mx = dplyr::if_else(
      `Mx1 Positive pressure` == 1 | `Mx 2 steriod` == 1 | `Mx 3 decrease Fio2` == 1 |
        `Mx3 Positive pressure` == 1 | `M3 negtiveIO` == 1, 1, 0
    )
  )

# 11b. Extubation-failure proportion by protocol adherence, among high-risk
# patients.
df_highrisk <- dplyr::filter(df_protocol, risk_group == "High risk")
fisher_protocol <- fisher.test(table(df_highrisk$protocol, df_highrisk$outcome_of_extubation))

protocol_summary <- df_highrisk |>
  dplyr::group_by(protocol) |>
  dplyr::summarise(
    n          = dplyr::n(),
    n_failure  = sum(outcome_of_extubation == 1),
    proportion = mean(outcome_of_extubation == 1),
    fraction   = sprintf("%d / %d", n_failure, n),
    .groups = "drop"
  ) |>
  tidyr::complete(protocol = c(0, 1),
                  fill = list(n = 0, n_failure = 0, proportion = 0)) |>
  dplyr::mutate(protocol = factor(protocol, levels = c(0, 1),
                labels = c("Not adherent", "Adherent")))

protocol_plot <- ggplot(protocol_summary, aes(x = protocol, y = proportion, fill = protocol)) +
  geom_col(width = 0.7) +
  geom_text(aes(y = 0.01, label = sprintf("n = %s", fraction)), size = 4) +
  annotate("text", x = 1.5, y = 0.36,
           label = sprintf("P-value = %s", format_pvalue(fisher_protocol$p.value)),
           size = 4) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(name = NULL,
    values = c("Not adherent" = "#92cced", "Adherent" = "#ea801c")) +
  labs(x = "Protocol adherence", y = "Proportion of extubation failure") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line   = element_line(),
    axis.ticks  = element_line(),
    legend.position = "bottom",
    axis.title.x = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(vjust = 0.5, face = "bold", size = 12, margin = margin(r = 10)),
    axis.text    = element_text(size = 10)
  )
protocol_plot

# 11c. Extubation-failure proportion by whether any management was received,
# across the whole validation cohort.
fisher_mx <- fisher.test(table(df_protocol$get_mx, df_protocol$outcome_of_extubation))

mgmt_summary <- df_protocol |>
  dplyr::group_by(get_mx) |>
  dplyr::summarise(
    n          = dplyr::n(),
    n_failure  = sum(outcome_of_extubation, na.rm = TRUE),
    proportion = mean(outcome_of_extubation, na.rm = TRUE),
    fraction   = sprintf("%d / %d", n_failure, n),
    .groups = "drop"
  ) |>
  dplyr::mutate(get_mx = factor(get_mx, levels = c(0, 1),
                labels = c("No management", "Received management")))

mgmt_overall_plot <- ggplot(mgmt_summary, aes(x = get_mx, y = proportion, fill = get_mx)) +
  geom_col(width = 0.7) +
  geom_text(aes(y = 0.01, label = sprintf("n = %s", fraction)), size = 4) +
  annotate("text", x = 1.5, y = 0.36,
           label = sprintf("P-value = %s", format_pvalue(fisher_mx$p.value)),
           size = 4) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(
    values = c("No management" = "#92cced", "Received management" = "#ea801c")
  ) +
  labs(x = "Management", y = "Proportion of extubation failure") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line   = element_line(),
    axis.ticks  = element_line(),
    legend.position = "none",
    axis.title.x = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(t = 10)),
    axis.title.y = element_text(vjust = 0.5, face = "bold", size = 12, margin = margin(r = 10)),
    axis.text    = element_text(size = 10)
  )
mgmt_overall_plot


# ---- 12. Checks and saved results -------------------------------------------

# Cohort sizes, event counts and Table 4 are not affected by the score or
# membership-model corrections, so a rerun must reproduce them exactly.
stopifnot(
  sum(df_member$dataset == 0) == 352,
  sum(df_member$dataset == 1) == 142,
  sum(df_member$outcome_of_extubation[df_member$dataset == 0]) == 40,
  sum(df_member$outcome_of_extubation[df_member$dataset == 1]) == 11,
  round(diag_val$detail$est[diag_val$detail$statistic == "se"], 4) == 0.8182,
  round(diag_val$detail$est[diag_val$detail$statistic == "pv.neg"], 4) == 0.9778
)

# Every number quoted in the manuscript, in one place.
results <- list(
  membership_auroc = as.numeric(case_mix_roc$ci),
  score_by_dataset = df_member |>
    dplyr::group_by(dataset) |>
    dplyr::summarise(
      score_mean = mean(score), score_sd = sd(score),
      lp_mean = mean(lp_score), lp_sd = sd(lp_score),
      .groups = "drop"
    ),
  score_ttest_p  = t.test(score ~ dataset, data = df_member)$p.value,
  auroc_val      = as.numeric(roc_val$ci),
  auroc_dev      = as.numeric(roc_dev$ci),
  cal_slope      = cal_slope,
  cal_citl       = cal_citl,
  mean_predicted = mean(cal$pred),
  observed_rate  = mean(cal$obs),
  table4         = table4
)
saveRDS(results, file.path(output_dir, "results.rds"))

# =============================================================================
# End of script
# =============================================================================
