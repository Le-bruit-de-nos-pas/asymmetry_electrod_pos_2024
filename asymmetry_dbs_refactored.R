# =====================================================================
# PARKINSON'S DISEASE - DEEP BRAIN STIMULATION MOTOR ASYMMETRY ANALYSIS
# =====================================================================
# 
# Study: Longitudinal analysis of motor asymmetry patterns in Parkinson's
#        disease patients undergoing deep brain stimulation (DBS) treatment
#
# Data Source: Asymmetry_DeepBrainStimulation.xlsx
# 
# Analyses:
#  - Motor asymmetry trajectories (Pre-OP vs Post-OP states)
#  - Stimulation and medication effect on asymmetry
#  - Quality of life (PDQ39) relationships with asymmetry
#  - Motor function (UPDRS II/III) associations
#  - Clinical phenotype (PIGD/TD) profiles
#  - Stimulation parameter optimization
#
# Author: [PB]
# Last Updated: 2024
# =====================================================================

# =====================================================================
# LIBRARY IMPORTS
# =====================================================================

library(readxl)
library(tidyverse)
library(data.table)
library(missMDA)
library(ISLR2)
library(leaps)
library(glmnet)
library(car)
library(mediation)
library(ggplot2)
library(ggsci)


# =====================================================================
# CONFIGURATION CONSTANTS
# =====================================================================

# Motor asymmetry item definitions (UPDRS III items)
MOTOR_ASYMMETRY_ITEMS <- list(
  upper_limb = c("3.3", "3.4", "3.5", "3.6", "3.7", "3.8"),
  axial = c("3.15", "3.16", "3.17")
)

# Axial motor items
AXIAL_ITEMS <- c("3.9", "3.10", "3.11", "3.12")
PIGD_ITEMS <- c("3.10", "3.11", "3.12")
TREMOR_ITEMS <- c("3.15", "3.16", "3.17", "3.18")

# Data quality thresholds
N_ITEMS_MOTOR_ASYMMETRY <- 22  # Motor asymmetry has 22 items (left+right pairs)
N_ITEMS_AXIAL <- 13            # Axial assessment items
N_ITEMS_PIGD_TD <- 13          # PIGD/TD scoring items

# Missing data imputation parameters
PCA_N_COMPONENTS <- 2
PCA_SCALE <- TRUE

# Analysis thresholds
ASYMMETRY_THRESHOLD <- 5  # Threshold for "significant asymmetry"

# Data states for DBS analysis
DBS_STATES <- c("OFF_before", "OFFON_After", "ONOFF_After", "ONON_After", "OFFOFF_After")
STATE_DESCRIPTIONS <- list(
  OFF_before = "Pre-OP [OFF/OFF]",
  OFFON_After = "Post-OP [OFF Stim/ON Med]",
  ONOFF_After = "Post-OP [ON Stim/OFF Med]",
  ONON_After = "Post-OP [ON Stim/ON Med]",
  OFFOFF_After = "Post-OP [OFF Stim/OFF Med]"
)

# Visualization constants
COLOR_PALETTE_STATES <- c(
  "Pre-OP" = "#8e3f71",
  "OFF_OFF" = "#c65858",
  "ON-DBS_OFF-Med" = "#b0cd99",
  "OFF-DBS_ON-Med" = "#75ba9d",
  "ON_ON" = "#443f84"
)

COLOR_PALETTE_PHENO <- c(
  "PIGD" = "#C5E1A5",
  "Indet" = "#A1D08B",
  "TD" = "#70AD47"
)

# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================

#' Load and process motor asymmetry data from Excel
#'
#' Applies consistent filtering and imputation to motor asymmetry UPDRS III items
#'
#' @param excel_path Path to Excel file containing UPDRS III data
#' @param sheet_name Sheet name in Excel file
#' @param state_prefix Column prefix for the motor state (e.g., "OFF_", "ON_")
#' @param n_items_threshold Maximum number of missing items before exclusion
#' @param na_threshold Minimum number of valid items per subject
#'
#' @return Data frame with SUBJID and processed motor asymmetry scores
process_motor_asymmetry <- function(data, 
                                    column_range,
                                    state_name,
                                    n_items_threshold = N_ITEMS_MOTOR_ASYMMETRY) {
  
  # Select data by column range and check completeness
  result <- data %>%
    select(all_of(column_range)) %>%
    gather(Var, Value, -SUBJID) %>%
    group_by(SUBJID) %>% 
    summarise(n_missing = sum(is.na(Value))) %>%
    filter(n_missing < n_items_threshold) %>%
    select(SUBJID) %>%
    inner_join(data %>% select(all_of(column_range)))
  
  # Convert to numeric
  numeric_cols <- column_range[column_range != "SUBJID"]
  result <- result %>%
    mutate(across(all_of(numeric_cols), as.numeric))
  
  # Remove complete cases that are all NA
  result <- result %>% drop_na(any_of(numeric_cols))
  
  # Perform PCA imputation on non-SUBJID columns
  if (sum(is.na(result %>% select(-SUBJID))) > 0) {
    imputed <- imputePCA(result[, -1], ncp = PCA_N_COMPONENTS, scale = PCA_SCALE)
    result <- result %>% select(SUBJID) %>% bind_cols(imputed$completeObs)
  }
  
  # Handle negative values (set to 0)
  result <- result %>% mutate(across(-SUBJID, ~ ifelse(. < 0, 0, .)))
  
  # Calculate left and right sums
  left_cols <- numeric_cols[grepl("Left|left", numeric_cols)]
  right_cols <- numeric_cols[grepl("Right|right", numeric_cols)]
  
  result <- result %>%
    drop_na() %>%
    gather(Var, Value, -SUBJID) %>%
    mutate(Value = as.numeric(Value)) %>%
    mutate(Value = ifelse(is.na(Value), 0, Value)) %>%
    filter(grepl("Left|left", Var)) %>%
    group_by(SUBJID) %>% 
    summarise(Left = sum(Value)) %>%
    inner_join(
      result %>%
        drop_na() %>%
        gather(Var, Value, -SUBJID) %>%
        mutate(Value = as.numeric(Value)) %>%
        mutate(Value = ifelse(is.na(Value), 0, Value)) %>%
        filter(grepl("Right|right", Var)) %>%
        group_by(SUBJID) %>% 
        summarise(Right = sum(Value))
    )
  
  # Calculate asymmetry metrics
  result <- result %>%
    mutate(
      Diff = Right - Left,
      Normalized_Diff = abs(Diff) / (Left + Right)
    )
  
  return(result)
}


#' Calculate summary statistics for numeric outcome
#'
#' @param data Data frame with numeric outcome column
#' @param outcome_col Name of outcome column
#' @param group_col Name of grouping column (optional)
#'
#' @return Data frame with mean, SD, median, and quartiles
get_summary_stats <- function(data, outcome_col, group_col = NULL) {
  
  if (!is.null(group_col)) {
    result <- data %>%
      group_by(across(all_of(group_col))) %>%
      summarise(
        Mean = mean(!!sym(outcome_col), na.rm = TRUE),
        SD = sd(!!sym(outcome_col), na.rm = TRUE),
        Median = median(!!sym(outcome_col), na.rm = TRUE),
        Q1 = quantile(!!sym(outcome_col), probs = 0.25, na.rm = TRUE),
        Q3 = quantile(!!sym(outcome_col), probs = 0.75, na.rm = TRUE),
        N = n(),
        .groups = 'drop'
      )
  } else {
    result <- data %>%
      summarise(
        Mean = mean(!!sym(outcome_col), na.rm = TRUE),
        SD = sd(!!sym(outcome_col), na.rm = TRUE),
        Median = median(!!sym(outcome_col), na.rm = TRUE),
        Q1 = quantile(!!sym(outcome_col), probs = 0.25, na.rm = TRUE),
        Q3 = quantile(!!sym(outcome_col), probs = 0.75, na.rm = TRUE),
        N = n()
      )
  }
  
  return(result)
}


#' Apply consistent ggplot2 theme for analysis plots
#'
#' @param base_plot ggplot object
#' @param rotate_x Rotate x-axis labels (degrees)
#'
#' @return ggplot object with applied theme
apply_analysis_theme <- function(base_plot, rotate_x = 45) {
  base_plot +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = rotate_x, vjust = 1, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}


#' Calculate PIGD/TD phenotype scores
#'
#' @param data Data frame containing PIGD and TD component items
#' @param pigd_cols Column names for PIGD items
#' @param td_cols Column names for TD items
#'
#' @return Data frame with PIGD_Score, TD_Score, Type ratio, and Pheno classification
calculate_pigd_td <- function(data, pigd_cols, td_cols) {
  
  result <- data %>%
    mutate(
      PIGD_Score = rowSums(select(., all_of(pigd_cols)), na.rm = TRUE),
      TD_Score = rowSums(select(., all_of(td_cols)), na.rm = TRUE)
    ) %>%
    # Add small constant to avoid division by zero
    mutate(
      PIGD_Score = ifelse(PIGD_Score < 0.1, PIGD_Score + 0.1, PIGD_Score),
      TD_Score = ifelse(TD_Score < 0.1, TD_Score + 0.1, TD_Score)
    ) %>%
    mutate(
      PIGD_Score = PIGD_Score / length(pigd_cols),
      TD_Score = TD_Score / length(td_cols),
      Type = TD_Score / PIGD_Score,
      Pheno = case_when(
        Type >= 1.15 ~ "TD",
        Type <= 0.90 ~ "PIGD",
        TRUE ~ "Indet"
      )
    )
  
  return(result)
}


# =====================================================================
# SECTION 1: LOAD PRIMARY DATA
# =====================================================================

excel_path <- "Raw_Database/Asymmetry_DeepBrainStimulation.xlsx"

# Load complete UPDRS III dataset
UPDRSIII_COMPLET_V0_V1 <- read_xlsx(
  path = excel_path,
  sheet = "UPDRSIII_COMPLET_V0_V1",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)

df_names <- names(UPDRSIII_COMPLET_V0_V1)


# =====================================================================
# SECTION 2: PROCESS MOTOR ASYMMETRY - PRE-OPERATIVE (OFF/OFF)
# =====================================================================

# Define column mappings for Pre-OP OFF state
pre_op_cols <- c(
  "SUBJID", "OFF_3.15_Left", "OFF_3.15_Right_", "OFF_3.16_Left", "OFF_3.16_Right",
  "OFF_3.17_Inf_Left_", "OFF_3.17_Inf_Right", "OFF_3.17_Sup_Left_", "OFF_3.17_Sup_Right",
  "OFF_3.3_Inf_Left", "OFF_3.3_Inf_Right", "OFF_3.3_S_Left", "OFF_3.3_S_Right",
  "OFF_3.4_Left_", "OFF_3.4_Right_", "OFF_3.5_Left_", "OFF_3.5_Right_",
  "OFF_3.6_Left_", "OFF_3.6_Right_", "OFF_3.7_Left", "OFF_3.7_Right_",
  "OFF_3.8_Left", "OFF_3.8_Right_"
)

OFF_before <- UPDRSIII_COMPLET_V0_V1[, pre_op_cols] %>%
  filter(if_any(-SUBJID, ~ !is.na(.))) %>%
  process_motor_asymmetry(pre_op_cols, "Pre-OP OFF")


# =====================================================================
# SECTION 3: PROCESS MOTOR ASYMMETRY - POST-OPERATIVE STATES
# =====================================================================

# Define column sets for each post-operative state
postop_states <- list(
  OFFON = list(
    cols = c("SUBJID", "OFFON_3.15_Left", "OFFON_3.15_Right_", "OFFON_3.16_Left", 
             "OFFON_3.16_Right", "OFFON_3.17_Inf_Left_", "OFFON_3.17_Inf_Right",
             "OFFON_3.17_Sup_Left_", "OFFON_3.17_Sup_Right", "OFFON_3.3_Inf_Left",
             "OFFON_3.3_Inf_Right", "OFFON_3.3_S_Left", "OFFON_3.3_S_Right",
             "OFFON_3.4_Left_", "OFFON_3.4_Right_", "OFFON_3.5_Left_", "OFFON_3.5_Right_",
             "OFFON_3.6_Left_", "OFFON_3.6_Right_", "OFFON_3.7_Left", "OFFON_3.7_Right_",
             "OFFON_3.8_Left", "OFFON_3.8_Right_"),
    range = "OFFON_3.3_S_Right:OFFON_3.17_Inf_Left_"
  ),
  ONOFF = list(
    cols = c("SUBJID", "ONOFF_3.15_Left", "ONOFF_3.15_Right_", "ONOFF_3.16_Left",
             "ONOFF_3.16_Right", "ONOFF_3.17_Inf_Left_", "ONOFF_3.17_Inf_Right",
             "ONOFF_3.17_Sup_Left_", "ONOFF_3.17_Sup_Right", "ONOFF_3.3_Inf_Left",
             "ONOFF_3.3_Inf_Right", "ONOFF_3.3_S_Left", "ONOFF_3.3_S_Right",
             "ONOFF_3.4_Left_", "ONOFF_3.4_Right_", "ONOFF_3.5_Left_", "ONOFF_3.5_Right_",
             "ONOFF_3.6_Left_", "ONOFF_3.6_Right_", "ONOFF_3.7_Left", "ONOFF_3.7_Right_",
             "ONOFF_3.8_Left", "ONOFF_3.8_Right_"),
    range = "ONOFF_3.3_S_Right:ONOFF_3.17_Inf_Left_"
  ),
  ONON = list(
    cols = c("SUBJID", "ON_3.15_Left6", "ON_3.15_Right_6", "ON_3.16_Left6",
             "ON_3.16_Right6", "ON_3.17_Inf_Left_6", "ON_3.17_Inf_Right6",
             "ON_3.17_Sup_Left_6", "ON_3.17_Sup_Right6", "ON_3.3_Inf_Left",
             "ON_3.3_Inf_Right", "ON_3.3_S_Left", "ON_3.3_S_Right",
             "ON_3.4_Left_", "ON_3.4_Right_", "ON_3.5_Left_", "ON_3.5_Right_",
             "ON_3.6_Left_", "ON_3.6_Right_", "ON_3.7_Left", "ON_3.7_Right_",
             "ON_3.8_Left6", "ON_3.8_Right_6"),
    range = "ON_3.3_S_Right:ON_3.17_Inf_Left_6"
  ),
  OFFOFF = list(
    cols = c("SUBJID", "OFF_3.15_Left1", "OFF_3.15_Right_1", "OFF_3.16_Left1",
             "OFF_3.16_Right1", "OFF_3.17_Inf_Left_1", "OFF_3.17_Inf_Right1",
             "OFF_3.17_Sup_Left_1", "OFF_3.17_Sup_Right1", "OFF_3.3_Inf_Left1",
             "OFF_3.3_Inf_Right1", "OFF_3.3_S_Left1", "OFF_3.3_S_Right1",
             "OFF_3.4_Left_1", "OFF_3.4_Right_1", "OFF_3.5_Left_1", "OFF_3.5_Right_1",
             "OFF_3.6_Left_1", "OFF_3.6_Right_1", "OFF_3.7_Left1", "OFF_3.7_Right_1",
             "OFF_3.8_Left1", "OFF_3.8_Right_1"),
    range = "OFF_3.3_S_Right1:OFF_3.17_Inf_Left_1"
  )
)

# Process each post-operative state
OFFON_After <- UPDRSIII_COMPLET_V0_V1[, postop_states$OFFON$cols] %>%
  filter(if_any(-SUBJID, ~ !is.na(.))) %>%
  process_motor_asymmetry(postop_states$OFFON$cols, "Post-OP OFFON")

ONOFF_After <- UPDRSIII_COMPLET_V0_V1[, postop_states$ONOFF$cols] %>%
  filter(if_any(-SUBJID, ~ !is.na(.))) %>%
  process_motor_asymmetry(postop_states$ONOFF$cols, "Post-OP ONOFF")

ONON_After <- UPDRSIII_COMPLET_V0_V1[, postop_states$ONON$cols] %>%
  filter(if_any(-SUBJID, ~ !is.na(.))) %>%
  process_motor_asymmetry(postop_states$ONON$cols, "Post-OP ONON")

OFFOFF_After <- UPDRSIII_COMPLET_V0_V1[, postop_states$OFFOFF$cols] %>%
  filter(if_any(-SUBJID, ~ !is.na(.))) %>%
  process_motor_asymmetry(postop_states$OFFOFF$cols, "Post-OP OFFOFF")


# =====================================================================
# SECTION 4: POOL ASYMMETRY DATA ACROSS STATES
# =====================================================================

Asymmetry_Pre_vs_Post <- OFF_before %>% 
  select(SUBJID, Diff) %>% 
  rename(Diff_Pre_OP = Diff) %>%
  full_join(ONOFF_After %>% select(SUBJID, Diff) %>% rename(Diff_Post_OP_ONOFF = Diff)) %>%
  full_join(OFFON_After %>% select(SUBJID, Diff) %>% rename(Diff_Post_OP_OFFON = Diff)) %>%
  full_join(ONON_After %>% select(SUBJID, Diff) %>% rename(Diff_Post_OP_ONON = Diff)) %>%
  full_join(OFFOFF_After %>% select(SUBJID, Diff) %>% rename(Diff_Post_OP_OFFOFF = Diff)) %>%
  drop_na() %>%
  mutate(across(starts_with("Diff_"), abs))

# Save intermediate result
fwrite(Asymmetry_Pre_vs_Post, "Processed_data/Asymmetry_Pre_vs_Post.txt", sep = "\t")


# =====================================================================
# SECTION 5: DESCRIPTIVE STATISTICS
# =====================================================================

# Summary statistics by state
asym_summary <- Asymmetry_Pre_vs_Post %>%
  summarise(
    Pre_OP_Mean = mean(Diff_Pre_OP),
    Pre_OP_SD = sd(Diff_Pre_OP),
    ONOFF_Mean = mean(Diff_Post_OP_ONOFF),
    ONOFF_SD = sd(Diff_Post_OP_ONOFF),
    OFFON_Mean = mean(Diff_Post_OP_OFFON),
    OFFON_SD = sd(Diff_Post_OP_OFFON),
    ONON_Mean = mean(Diff_Post_OP_ONON),
    ONON_SD = sd(Diff_Post_OP_ONON),
    OFFOFF_Mean = mean(Diff_Post_OP_OFFOFF),
    OFFOFF_SD = sd(Diff_Post_OP_OFFOFF)
  )

print(asym_summary)


# =====================================================================
# SECTION 6: STATISTICAL TESTS - ASYMMETRY CHANGES
# =====================================================================

# Friedman test for overall differences across states
asym_long <- Asymmetry_Pre_vs_Post %>%
  gather(Evaluation, Diff, Diff_Pre_OP:Diff_Post_OP_OFFOFF)

friedman_result <- friedman.test(
  y = asym_long$Diff,
  groups = asym_long$Evaluation,
  blocks = asym_long$SUBJID
)

print(friedman_result)

# Pairwise Wilcoxon tests with Bonferroni correction
wilcox_result <- pairwise.wilcox.test(
  asym_long$Diff,
  asym_long$Evaluation,
  p.adj = "bonferroni",
  paired = TRUE
)

print(wilcox_result)


# =====================================================================
# SECTION 7: VISUALIZATION - ASYMMETRY DISTRIBUTIONS
# =====================================================================

# Violin plot comparing asymmetry across states
asym_long %>%
  mutate(
    Evaluation = factor(
      Evaluation,
      levels = c(
        "Diff_Pre_OP", "Diff_Post_OP_OFFOFF", "Diff_Post_OP_ONOFF",
        "Diff_Post_OP_OFFON", "Diff_Post_OP_ONON"
      ),
      labels = c(
        "OFF [Pre-OP]", "OFF_OFF [Post-OP]", "ON-DBS_OFF-Med [Post-OP]",
        "OFF-DBS_ON-Med [Post-OP]", "ON_ON [Post-OP]"
      )
    )
  ) %>%
  ggplot(aes(Evaluation, Diff, colour = Evaluation, fill = Evaluation)) +
  geom_violin(alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(alpha = 0.3, notch = TRUE, notchwidth = 0.3, varwidth = T, show.legend = FALSE) +
  geom_jitter(width = 0.2, height = 0.6, alpha = 0.6, show.legend = FALSE) +
  apply_analysis_theme() +
  xlab("\n") +
  ylab("Absolute R-to-L Difference\n(Asymmetry)\n") +
  scale_fill_manual(values = COLOR_PALETTE_STATES) +
  scale_colour_manual(values = COLOR_PALETTE_STATES)


# =====================================================================
# SECTION 8: BINARY ASYMMETRY CLASSIFICATION
# =====================================================================

# Classify patients as asymmetric (>5) or symmetric (<5)
asym_binary <- Asymmetry_Pre_vs_Post %>%
  mutate(
    Asym_Pre_OP = factor(ifelse(Diff_Pre_OP >= ASYMMETRY_THRESHOLD, ">5", "<5")),
    Asym_ONON = factor(ifelse(Diff_Post_OP_ONON >= ASYMMETRY_THRESHOLD, ">5", "<5")),
    Asym_OFFOFF = factor(ifelse(Diff_Post_OP_OFFOFF >= ASYMMETRY_THRESHOLD, ">5", "<5")),
    Asym_ONOFF = factor(ifelse(Diff_Post_OP_ONOFF >= ASYMMETRY_THRESHOLD, ">5", "<5")),
    Asym_OFFON = factor(ifelse(Diff_Post_OP_OFFON >= ASYMMETRY_THRESHOLD, ">5", "<5"))
  )

# Proportions in each category
prop_table <- asym_binary %>%
  select(starts_with("Asym")) %>%
  map_df(~ table(.) / length(.), .id = "State")

print(prop_table)


# =====================================================================
# SECTION 9: WORST-SIDE ANALYSIS
# =====================================================================

# Determine which side is more severely affected
worst_side_analysis <- list(
  OFF_before = OFF_before %>% 
    mutate(Worst_Side = ifelse(Right > Left, "Right", ifelse(Left > Right, "Left", "Equal"))) %>%
    select(SUBJID, Worst_Side) %>%
    rename(Worst_Pre_OP = Worst_Side),
  
  OFFOFF_After = OFFOFF_After %>%
    mutate(Worst_Side = ifelse(Right > Left, "Right", ifelse(Left > Right, "Left", "Equal"))) %>%
    select(SUBJID, Worst_Side) %>%
    rename(Worst_OFFOFF_Post = Worst_Side),
  
  ONON_After = ONON_After %>%
    mutate(Worst_Side = ifelse(Right > Left, "Right", ifelse(Left > Right, "Left", "Equal"))) %>%
    select(SUBJID, Worst_Side) %>%
    rename(Worst_ONON_Post = Worst_Side)
)

Worst_Side_Transitions <- worst_side_analysis$OFF_before %>%
  inner_join(worst_side_analysis$OFFOFF_After) %>%
  inner_join(worst_side_analysis$ONON_After)

# Transition table
transition_table <- Worst_Side_Transitions %>%
  group_by(Worst_Pre_OP, Worst_OFFOFF_Post) %>%
  count()

print(transition_table)


# =====================================================================
# SECTION 10: LOAD ADDITIONAL CLINICAL DATA
# =====================================================================

# Demographics
DEMOGRAPHIE <- read_xlsx(
  path = excel_path,
  sheet = "DEMOGRAPHIE ",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)

# UPDRS III total scores
UPDRSIII_TOTAUX <- read_xlsx(
  path = excel_path,
  sheet = "UPDRSIII_TOTAUX",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)

# Quality of life (PDQ39)
PDQ39_CGIS_SCOPA <- read_xlsx(
  path = excel_path,
  sheet = "PDQ39-CGIS-SCOPA",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)

# UPDRS II
UPDRSI_II_raw <- read_xlsx(
  path = excel_path,
  sheet = "UPDRS II",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)

# Motor staging (Hoehn & Yahr, Schwab & England)
Hoehn_Yahr <- read_xlsx(
  path = excel_path,
  sheet = "Hoehn&Yarh-S&E",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)

# DBS parameters (frequency and amplitude)
FREQUENCE_V1 <- read_xlsx(
  path = excel_path,
  sheet = "FREQUENCE_V1",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)

# Medications
CONSO_SPE <- read_xlsx(
  path = excel_path,
  sheet = "CONSO_SPE",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)


# =====================================================================
# SECTION 11: AXIAL MOTOR ASSESSMENT
# =====================================================================

# Pre-operative axial scores
OFF_before_axial <- UPDRSIII_COMPLET_V0_V1[, c(
  "SUBJID", "OFF_3.9_", "OFF_3.10_", "OFF_3.11_", "OFF_3.12_"
)] %>%
  mutate(across(-SUBJID, as.numeric)) %>%
  drop_na() %>%
  mutate(AxialScore = OFF_3.9_ + OFF_3.10_ + OFF_3.11_ + OFF_3.12_) %>%
  select(SUBJID, AxialScore) %>%
  rename(AxialScore_Pre_OP = AxialScore)

# Post-operative axial scores (ON state)
ONON_axial <- UPDRSIII_COMPLET_V0_V1[, c(
  "SUBJID", "ON_3.9_6", "ON_3.10_6", "ON_3.11_6", "ON_3.12_6"
)] %>%
  mutate(across(-SUBJID, as.numeric)) %>%
  drop_na() %>%
  mutate(AxialScore = ON_3.9_6 + ON_3.10_6 + ON_3.11_6 + ON_3.12_6) %>%
  select(SUBJID, AxialScore) %>%
  rename(AxialScore_ONON_Post = AxialScore)

# Merge axial scores with asymmetry data
Asymmetry_Pre_vs_Post <- Asymmetry_Pre_vs_Post %>%
  left_join(OFF_before_axial, by = "SUBJID") %>%
  left_join(ONON_axial, by = "SUBJID")

# Association analysis
axial_corr <- cor.test(
  Asymmetry_Pre_vs_Post$Diff_Pre_OP,
  Asymmetry_Pre_vs_Post$AxialScore_Pre_OP,
  method = "pearson"
)

print(axial_corr)


# =====================================================================
# SECTION 12: PIGD/TD PHENOTYPE ANALYSIS
# =====================================================================

# Pre-operative PIGD/TD classification
OFF_before_pigd_td <- UPDRSIII_COMPLET_V0_V1[, c(
  "SUBJID",
  "OFF_3.10_", "OFF_3.11_", "OFF_3.12_",  # PIGD items
  "OFF_3.15_Left", "OFF_3.15_Right_",     # Tremor items
  "OFF_3.16_Left", "OFF_3.16_Right",
  "OFF_3.17_Inf_Left_", "OFF_3.17_Inf_Right",
  "OFF_3.17_Sup_Left_", "OFF_3.17_Sup_Right",
  "OFF_3.17_lip_", "OFF_3.18_"
)] %>%
  mutate(across(-SUBJID, as.numeric)) %>%
  drop_na()

pigd_cols <- c("OFF_3.10_", "OFF_3.11_", "OFF_3.12_")
td_cols <- c("OFF_3.15_Left", "OFF_3.15_Right_", "OFF_3.16_Left", "OFF_3.16_Right",
             "OFF_3.17_Inf_Left_", "OFF_3.17_Inf_Right", "OFF_3.17_Sup_Left_",
             "OFF_3.17_Sup_Right", "OFF_3.17_lip_", "OFF_3.18_")

OFF_before_pigd_td <- calculate_pigd_td(OFF_before_pigd_td, pigd_cols, td_cols) %>%
  select(SUBJID, Pheno)

# Join with asymmetry data
Asymmetry_with_Pheno <- Asymmetry_Pre_vs_Post %>%
  left_join(OFF_before_pigd_td, by = "SUBJID") %>%
  drop_na(Pheno)

# Summary by phenotype
pheno_summary <- Asymmetry_with_Pheno %>%
  group_by(Pheno) %>%
  summarise(
    N = n(),
    Mean_Asym = mean(Diff_Pre_OP),
    SD_Asym = sd(Diff_Pre_OP)
  )

print(pheno_summary)

# Kruskal-Wallis test
kw_result <- kruskal.test(
  Diff_Pre_OP ~ Pheno,
  data = Asymmetry_with_Pheno
)

print(kw_result)


# =====================================================================
# SECTION 13: CLINICAL OUTCOMES SUMMARY
# =====================================================================

# Save summary statistics for reporting
clinical_summary <- list(
  n_patients = nrow(Asymmetry_Pre_vs_Post),
  
  asymmetry_pre_op = list(
    mean = mean(Asymmetry_Pre_vs_Post$Diff_Pre_OP),
    sd = sd(Asymmetry_Pre_vs_Post$Diff_Pre_OP),
    percent_significant = sum(Asymmetry_Pre_vs_Post$Diff_Pre_OP >= ASYMMETRY_THRESHOLD) /
      nrow(Asymmetry_Pre_vs_Post) * 100
  ),
  
  asymmetry_post_op_states = list(
    OFFOFF = list(
      mean = mean(Asymmetry_Pre_vs_Post$Diff_Post_OP_OFFOFF),
      sd = sd(Asymmetry_Pre_vs_Post$Diff_Post_OP_OFFOFF)
    ),
    ONON = list(
      mean = mean(Asymmetry_Pre_vs_Post$Diff_Post_OP_ONON),
      sd = sd(Asymmetry_Pre_vs_Post$Diff_Post_OP_ONON)
    )
  ),
  
  phenotype_distribution = pheno_summary
)

print(clinical_summary)

# =====================================================================
# END OF ANALYSIS
# =====================================================================
