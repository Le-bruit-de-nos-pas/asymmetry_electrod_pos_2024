# ============================================================================
# Analysis of Parkinson's Disease Motor Asymmetry and Related Outcomes
# ============================================================================
# Analysis of motor asymmetry patterns across visits (Year 1, 3, and 5) in
# Parkinson's disease patients undergoing deep brain stimulation.
# Examines relationships with axial scores, quality of life (PDQ39), and
# baseline demographic/clinical predictors.
# ============================================================================

# ============================================================================
# SETUP AND CONFIGURATION
# ============================================================================

library(readxl)
library(tidyverse)
library(data.table)
library(missMDA)
library(ISLR2)
library(leaps)
library(glmnet)
library(car)
library(mediation)

# Configuration constants
DATA_DIR <- "Raw_Database"
VISITED_ITEMS_THRESHOLD <- 22  # Min items to exclude subject
THRESHOLD_NEGATIVE_VALUES <- 0  # Minimum allowed value for motor scores
N_PCA_COMPONENTS <- 2
ASYMMETRY_OUTLIER_THRESHOLD <- 25
MOTOR_ASYMMETRY_COLS_V0V1 <- c(
  "ON_3.15_Left6", "ON_3.15_Right_6", "ON_3.16_Left6", "ON_3.16_Right6",
  "ON_3.17_Inf_Left_6", "ON_3.17_Inf_Right6", "ON_3.17_Sup_Left_6",
  "ON_3.17_Sup_Right6", "ON_3.3_Inf_Left", "ON_3.3_Inf_Right",
  "ON_3.3_S_Left", "ON_3.3_S_Right", "ON_3.4_Left_", "ON_3.4_Right_",
  "ON_3.5_Left_", "ON_3.5_Right_", "ON_3.6_Left_", "ON_3.6_Right_",
  "ON_3.7_Left", "ON_3.7_Right_", "ON_3.8_Left6", "ON_3.8_Right_6"
)
AXIAL_ITEMS_INDICES <- c("3.9", "3.10", "3.11", "3.12")
COLOR_PALETTE <- c("Year 1" = "#bbbbbb", "Year 3" = "#4689cc", "Year 5" = "#2841b0")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Load and preprocess motor asymmetry data from V0/V1
#'
#' @param file_path Path to Excel file
#' @param sheet_name Name of sheet containing data
#'
#' @return Processed data frame with Left, Right, Diff columns
load_and_process_v0v1_asymmetry <- function(file_path, sheet_name) {
  raw_data <- read_xlsx(
    path = file_path,
    sheet = sheet_name,
    skip = 0,
    col_types = "text",
    trim_ws = TRUE
  )

  # Select relevant columns
  cols_to_keep <- c("SUBJID", MOTOR_ASYMMETRY_COLS_V0V1)
  selected_cols <- intersect(cols_to_keep, names(raw_data))
  data <- raw_data[which(names(raw_data) %in% selected_cols)]
  data <- data[-1, ] # Remove header row

  # Filter subjects with adequate data
  data <- data %>%
    gather(variable, value, -SUBJID) %>%
    group_by(SUBJID) %>%
    summarise(n_missing = sum(is.na(value))) %>%
    filter(n_missing < VISITED_ITEMS_THRESHOLD) %>%
    select(SUBJID) %>%
    inner_join(data)

  # Convert to numeric and impute missing values
  data <- data %>%
    mutate(across(-SUBJID, as.numeric))

  # Impute missing values using PCA
  imputed <- imputePCA(data[, -1], ncp = N_PCA_COMPONENTS, scale = TRUE)
  data <- data %>%
    select(SUBJID) %>%
    bind_cols(imputed$completeObs)

  # Ensure all values are non-negative
  data[data < THRESHOLD_NEGATIVE_VALUES] <- THRESHOLD_NEGATIVE_VALUES

  # Calculate laterality scores
  left_cols <- names(data)[grepl("Left", names(data))]
  right_cols <- names(data)[grepl("Right", names(data))]

  data <- data %>%
    drop_na() %>%
    gather(variable, value, -SUBJID) %>%
    mutate(value = as.numeric(value), value = ifelse(is.na(value), 0, value)) %>%
    group_by(SUBJID, variable) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    drop_na()

  left_score <- data %>%
    select(SUBJID, all_of(left_cols)) %>%
    rowwise() %>%
    mutate(Left = sum(c_across(-SUBJID), na.rm = TRUE)) %>%
    select(SUBJID, Left)

  right_score <- data %>%
    select(SUBJID, all_of(right_cols)) %>%
    rowwise() %>%
    mutate(Right = sum(c_across(-SUBJID), na.rm = TRUE)) %>%
    select(SUBJID, Right)

  result <- left_score %>%
    inner_join(right_score) %>%
    mutate(Diff = Right - Left)

  return(result)
}

#' Calculate and format visit-based summary statistics
#'
#' @param data Data frame with VISIT column and numeric outcome
#' @param outcome_col Name of outcome column
#'
#' @return Data frame with summary statistics by visit
get_visit_summary_stats <- function(data, outcome_col) {
  data %>%
    group_by(VISIT) %>%
    summarise(
      mean = mean(get(outcome_col), na.rm = TRUE),
      sd = sd(get(outcome_col), na.rm = TRUE),
      median = median(get(outcome_col), na.rm = TRUE),
      Q1 = quantile(get(outcome_col), 0.25, na.rm = TRUE),
      Q3 = quantile(get(outcome_col), 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Apply common ggplot theme for consistent styling
#'
#' @return ggplot theme object
apply_plot_theme <- function() {
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.line = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, vjust = -0.5),
    axis.title.y = element_text(size = 12, vjust = -0.5),
    plot.margin = margin(5, 5, 5, 5, "pt")
  )
}

#' Convert visit numeric code to label
#'
#' @param visit Numeric visit code (1, 3, or 5)
#'
#' @return Character visit label
convert_visit_to_label <- function(visit) {
  dplyr::case_when(
    visit == 1 ~ "Year 1",
    visit == 3 ~ "Year 3",
    visit == 5 ~ "Year 5",
    TRUE ~ as.character(visit)
  )
}

# ============================================================================
# DATA LOADING AND PREPROCESSING
# ============================================================================

# Load Year 0/1 asymmetry data
v0v1_asymmetry_path <- file.path(DATA_DIR, "Asymmetry_DeepBrainStimulation.xlsx")
asymmetry_v0v1 <- load_and_process_v0v1_asymmetry(
  v0v1_asymmetry_path,
  "UPDRSIII_COMPLET_V0_V1"
) %>%
  mutate(VISIT = 1) %>%
  select(SUBJID, VISIT, Left, Right, Diff)

# Load Year 3/5 asymmetry data
v3v5_asymmetry_path <- file.path(DATA_DIR, "Raquel_Margherita_Juil 24.xlsx")
raw_v3v5 <- read_xlsx(
  path = v3v5_asymmetry_path,
  sheet = "UPDRSIII_COMPLET_V3_V5 ",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
)

motor_asymmetry_cols_v3v5 <- c(
  "ON_MSDROIT_RIG", "ON_MSGCHE_RIG", "ON_MIDROIT_RIG", "ON_MIGCHE_RIG",
  "ON_MS_DROIT_DOIGT", "ON_MSGCHE_DOIGT", "ON_MSDROIT_MAINS", "ON_MSGCHE_MAINS",
  "ON_MSDROIT_MA", "ON_MSGCHE_MA",
  "ON_MIDROIT_PIED", "ON_MIGCHE_PIED", "ON_MIDROIT_JAMBE", "ON_MIGCHE_JAMBE",
  "ON_TREMBLPOST_MSDROIT", "ON_TREMBLPOST_MSGCHE",
  "ON_TREMBLMAIN_MSGCHE", "ON_TREMBLMAIN_MSDROIT",
  "ON_MSDROIT_AMPLI_TREMBL", "ON_MSGCHE_AMPLI_TREMBL",
  "ON_MIDROIT_AMPLI_TREMBL", "ON_MIGCHE_AMPLI_TREMBL"
)

left_symptom_cols_v3v5 <- c(
  "ON_MSGCHE_RIG", "ON_MIGCHE_RIG", "ON_MSGCHE_DOIGT", "ON_MSGCHE_MAINS",
  "ON_MSGCHE_MA", "ON_MIGCHE_PIED", "ON_MIGCHE_JAMBE", "ON_TREMBLPOST_MSGCHE",
  "ON_TREMBLMAIN_MSGCHE", "ON_MSGCHE_AMPLI_TREMBL", "ON_MIGCHE_AMPLI_TREMBL"
)

right_symptom_cols_v3v5 <- c(
  "ON_MSDROIT_RIG", "ON_MIDROIT_RIG", "ON_MS_DROIT_DOIGT", "ON_MSDROIT_MAINS",
  "ON_MSDROIT_MA", "ON_MIDROIT_PIED", "ON_MIDROIT_JAMBE", "ON_TREMBLPOST_MSDROIT",
  "ON_TREMBLMAIN_MSDROIT", "ON_MSDROIT_AMPLI_TREMBL", "ON_MIDROIT_AMPLI_TREMBL"
)

# Process V3/V5 data
asymmetry_v3v5 <- raw_v3v5 %>%
  select(SUBJID, VISIT, all_of(motor_asymmetry_cols_v3v5)) %>%
  slice(-1) %>%
  mutate(across(-c(SUBJID, VISIT), as.numeric))

# Filter subjects with adequate data
asymmetry_v3v5 <- asymmetry_v3v5 %>%
  gather(variable, value, -SUBJID, -VISIT) %>%
  group_by(SUBJID, VISIT) %>%
  summarise(n_missing = sum(is.na(value)), .groups = "drop") %>%
  filter(n_missing < VISITED_ITEMS_THRESHOLD) %>%
  select(SUBJID) %>%
  inner_join(asymmetry_v3v5)

# Impute missing values
imputed_v3v5 <- imputePCA(
  asymmetry_v3v5[, -c(1, 2)],
  ncp = N_PCA_COMPONENTS,
  scale = TRUE
)
asymmetry_v3v5 <- asymmetry_v3v5 %>%
  select(SUBJID, VISIT) %>%
  bind_cols(imputed_v3v5$completeObs)

# Ensure non-negative values
asymmetry_v3v5[asymmetry_v3v5 < THRESHOLD_NEGATIVE_VALUES] <- THRESHOLD_NEGATIVE_VALUES

# Clean visit codes and calculate laterality
asymmetry_v3v5 <- asymmetry_v3v5 %>%
  mutate(VISIT = ifelse(grepl("V3", VISIT), 3, 5))

# Calculate left-side scores
left_scores_v3v5 <- asymmetry_v3v5 %>%
  drop_na() %>%
  gather(variable, value, -SUBJID, -VISIT) %>%
  mutate(value = as.numeric(value), value = ifelse(is.na(value), 0, value)) %>%
  filter(variable %in% left_symptom_cols_v3v5) %>%
  group_by(SUBJID, VISIT) %>%
  summarise(Left = sum(value), .groups = "drop")

# Calculate right-side scores
right_scores_v3v5 <- asymmetry_v3v5 %>%
  drop_na() %>%
  gather(variable, value, -SUBJID, -VISIT) %>%
  mutate(value = as.numeric(value), value = ifelse(is.na(value), 0, value)) %>%
  filter(variable %in% right_symptom_cols_v3v5) %>%
  group_by(SUBJID, VISIT) %>%
  summarise(Right = sum(value), .groups = "drop")

# Combine and calculate asymmetry
asymmetry_v3v5 <- left_scores_v3v5 %>%
  inner_join(right_scores_v3v5) %>%
  mutate(Diff = Right - Left)

# Merge V0/1 and V3/5 data
asymmetry_combined <- asymmetry_v0v1 %>%
  select(SUBJID) %>%
  distinct() %>%
  inner_join(asymmetry_v0v1) %>%
  bind_rows(asymmetry_v3v5)

# Filter subjects with complete data across all three visits
asymmetry_combined <- asymmetry_combined %>%
  group_by(SUBJID) %>%
  filter(n() == 3) %>%
  ungroup()

n_subjects <- n_distinct(asymmetry_combined$SUBJID)

# Save subject list for downstream analysis
subject_ids <- asymmetry_combined %>%
  select(SUBJID) %>%
  distinct()
fwrite(subject_ids, "Asym_415_allpats.txt")

# ============================================================================
# DESCRIPTIVE STATISTICS AND TESTING
# ============================================================================

# Left motor scores
left_summary <- get_visit_summary_stats(asymmetry_combined, "Left")
friedman.test(
  y = asymmetry_combined$Left,
  groups = asymmetry_combined$VISIT,
  blocks = asymmetry_combined$SUBJID
)
pairwise.wilcox.test(
  asymmetry_combined$Left,
  asymmetry_combined$VISIT,
  p.adj = "bonferroni",
  paired = TRUE
)

# Right motor scores
right_summary <- get_visit_summary_stats(asymmetry_combined, "Right")
friedman.test(
  y = asymmetry_combined$Right,
  groups = asymmetry_combined$VISIT,
  blocks = asymmetry_combined$SUBJID
)
pairwise.wilcox.test(
  asymmetry_combined$Right,
  asymmetry_combined$VISIT,
  p.adj = "bonferroni",
  paired = TRUE
)

# Asymmetry (absolute difference)
asymmetry_combined <- asymmetry_combined %>%
  mutate(Diff = abs(Diff))

asymmetry_summary <- get_visit_summary_stats(asymmetry_combined, "Diff")
friedman.test(
  y = asymmetry_combined$Diff,
  groups = asymmetry_combined$VISIT,
  blocks = asymmetry_combined$SUBJID
)
pairwise.wilcox.test(
  asymmetry_combined$Diff,
  asymmetry_combined$VISIT,
  p.adj = "bonferroni",
  paired = TRUE
)

# ============================================================================
# VISUALIZATION - LEFT SCORES
# ============================================================================

plot_data_left <- asymmetry_combined %>%
  mutate(VISIT_label = convert_visit_to_label(VISIT))

# Density plot
plot_left_density <- plot_data_left %>%
  ggplot(aes(Left, colour = VISIT_label, fill = VISIT_label)) +
  geom_density(alpha = 0.4) +
  labs(
    x = "\nLeft MDS UPDRS III",
    y = "Patient density\n"
  ) +
  scale_fill_manual(values = COLOR_PALETTE) +
  scale_colour_manual(values = COLOR_PALETTE) +
  xlim(0, 40) +
  apply_plot_theme()

ggsave(file = "left_density.svg", plot = plot_left_density, width = 4, height = 4)

# Boxplot with data points
plot_left_box <- plot_data_left %>%
  ggplot(aes(VISIT_label, Left, colour = VISIT_label, fill = VISIT_label)) +
  geom_boxplot(
    alpha = 0.6,
    notch = TRUE,
    width = 0.5,
    outlier.colour = NULL,
    outlier.fill = NULL,
    outlier.alpha = 0.00001
  ) +
  geom_jitter(alpha = 0.8, height = 0.5, shape = 1) +
  labs(
    x = "\nYearly Visit",
    y = "Left MDS UPDRS III\n"
  ) +
  scale_fill_manual(values = COLOR_PALETTE) +
  scale_colour_manual(values = COLOR_PALETTE) +
  ylim(0, 40) +
  apply_plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "left_boxplot.svg", plot = plot_left_box, width = 6, height = 6)

# ============================================================================
# VISUALIZATION - RIGHT SCORES
# ============================================================================

plot_data_right <- asymmetry_combined %>%
  mutate(VISIT_label = convert_visit_to_label(VISIT))

# Boxplot with data points
plot_right_box <- plot_data_right %>%
  ggplot(aes(VISIT_label, Right, colour = VISIT_label, fill = VISIT_label)) +
  geom_boxplot(
    alpha = 0.6,
    notch = TRUE,
    width = 0.5,
    outlier.colour = NULL,
    outlier.fill = NULL,
    outlier.alpha = 0.00001
  ) +
  geom_jitter(alpha = 0.8, height = 0.5, shape = 1) +
  labs(
    x = "\nYearly Visit",
    y = "Right MDS UPDRS III\n"
  ) +
  scale_fill_manual(values = COLOR_PALETTE) +
  scale_colour_manual(values = COLOR_PALETTE) +
  ylim(0, 40) +
  apply_plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "right_boxplot.svg", plot = plot_right_box, width = 6, height = 6)

# ============================================================================
# VISUALIZATION - ASYMMETRY SCORES
# ============================================================================

plot_data_asym <- asymmetry_combined %>%
  mutate(VISIT_label = convert_visit_to_label(VISIT))

# Boxplot with data points
plot_asym_box <- plot_data_asym %>%
  ggplot(aes(VISIT_label, Diff, colour = VISIT_label, fill = VISIT_label)) +
  geom_boxplot(
    alpha = 0.6,
    notch = TRUE,
    width = 0.5,
    outlier.colour = NULL,
    outlier.fill = NULL,
    outlier.alpha = 0.00001
  ) +
  geom_jitter(alpha = 0.8, height = 0.5, shape = 1) +
  labs(
    x = "\nYearly Visit",
    y = "Right-to-Left Asymmetry\nMDS UPDRS III\n"
  ) +
  scale_fill_manual(values = COLOR_PALETTE) +
  scale_colour_manual(values = COLOR_PALETTE) +
  ylim(0, 25) +
  apply_plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "asymmetry_boxplot.svg", plot = plot_asym_box, width = 6, height = 6)

# ============================================================================
# NORMALIZED ASYMMETRY
# ============================================================================

asymmetry_combined <- asymmetry_combined %>%
  mutate(
    Total = Left + Right,
    Normalized_Asymmetry = abs(Diff) / Total
  ) %>%
  drop_na(Normalized_Asymmetry)

# Filter subjects with complete data after normalization
asymmetry_combined <- asymmetry_combined %>%
  group_by(SUBJID) %>%
  filter(n() == 3) %>%
  ungroup()

# Normalized summary statistics
norm_asym_summary <- get_visit_summary_stats(
  asymmetry_combined,
  "Normalized_Asymmetry"
)

friedman.test(
  y = asymmetry_combined$Normalized_Asymmetry,
  groups = asymmetry_combined$VISIT,
  blocks = asymmetry_combined$SUBJID
)
pairwise.wilcox.test(
  asymmetry_combined$Normalized_Asymmetry,
  asymmetry_combined$VISIT,
  p.adj = "bonferroni",
  paired = TRUE
)

# Visualization
plot_norm_asym <- asymmetry_combined %>%
  mutate(VISIT_label = convert_visit_to_label(VISIT)) %>%
  ggplot(aes(VISIT_label, Normalized_Asymmetry, colour = VISIT_label, fill = VISIT_label)) +
  geom_boxplot(
    alpha = 0.6,
    notch = TRUE,
    width = 0.5,
    outlier.colour = NULL,
    outlier.fill = NULL,
    outlier.alpha = 0.00001
  ) +
  geom_jitter(alpha = 0.8, height = 0.01, shape = 1) +
  labs(
    x = "\nYearly Visit",
    y = "Normalized Right-to-Left Asymmetry\nMDS UPDRS III\n"
  ) +
  scale_fill_manual(values = COLOR_PALETTE) +
  scale_colour_manual(values = COLOR_PALETTE) +
  apply_plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "normalized_asymmetry.svg", plot = plot_norm_asym, width = 6, height = 6)

# ============================================================================
# AXIAL SCORE ANALYSIS
# ============================================================================

# Load axial score data
axial_v0v1 <- read_xlsx(
  path = v0v1_asymmetry_path,
  sheet = "UPDRSIII_COMPLET_V0_V1",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
) %>%
  select(SUBJID, ON_3.9_6, ON_3.10_6, ON_3.11_6, ON_3.12_6) %>%
  slice(-1) %>%
  mutate(
    across(-SUBJID, as.numeric),
    VISIT = 1,
    Axial_Score = ON_3.9_6 + ON_3.10_6 + ON_3.11_6 + ON_3.12_6
  ) %>%
  select(SUBJID, VISIT, Axial_Score)

axial_v3v5 <- read_xlsx(
  path = v3v5_asymmetry_path,
  sheet = "UPDRSIII_COMPLET_V3_V5 ",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
) %>%
  select(SUBJID, VISIT, ON_LEVER, ON_MARCHE, ON_FREEZING, ON_STAB_POST) %>%
  slice(-1) %>%
  mutate(
    across(-SUBJID, as.numeric),
    VISIT = ifelse(grepl("V3", VISIT), 3, 5),
    Axial_Score = ON_LEVER + ON_MARCHE + ON_FREEZING + ON_STAB_POST
  ) %>%
  select(SUBJID, VISIT, Axial_Score)

axial_combined <- bind_rows(axial_v0v1, axial_v3v5) %>%
  drop_na() %>%
  group_by(SUBJID) %>%
  filter(n() == 3) %>%
  ungroup()

# Merge with asymmetry data
axial_data <- axial_combined %>%
  inner_join(
    asymmetry_combined %>%
      select(SUBJID, VISIT, Diff, Total, Normalized_Asymmetry)
  )

# Statistical testing
axial_summary <- get_visit_summary_stats(axial_data, "Axial_Score")

friedman.test(
  y = axial_data$Axial_Score,
  groups = axial_data$VISIT,
  blocks = axial_data$SUBJID
)
pairwise.wilcox.test(
  axial_data$Axial_Score,
  axial_data$VISIT,
  p.adj = "bonferroni",
  paired = TRUE
)

# Correlation analysis
for (visit in c(1, 3, 5)) {
  data_subset <- axial_data %>% filter(VISIT == visit)
  cat("\nVisit:", visit, "\n")
  print(cor.test(
    abs(data_subset$Diff),
    data_subset$Axial_Score,
    method = "spearman"
  ))
}

# Linear regression
model_axial <- lm(Axial_Score ~ abs(Diff) + Total, data = axial_data)
summary(model_axial)

# Mediation analysis
med_model_axial <- lm(Total ~ Diff, data = axial_data)
outcome_model_axial <- lm(Axial_Score ~ Diff + Total, data = axial_data)
med_results_axial <- mediation::mediate(
  med_model_axial,
  outcome_model_axial,
  treat = "Diff",
  mediator = "Total",
  boot = TRUE,
  sims = 1000
)
summary(med_results_axial)

# Visualization
plot_axial <- axial_data %>%
  mutate(VISIT_label = convert_visit_to_label(VISIT)) %>%
  ggplot(aes(VISIT_label, Axial_Score, colour = VISIT_label, fill = VISIT_label)) +
  geom_boxplot(
    alpha = 0.6,
    notch = TRUE,
    width = 0.5,
    outlier.colour = NULL,
    outlier.fill = NULL,
    outlier.alpha = 0.00001
  ) +
  geom_jitter(alpha = 0.8, height = 0.5, shape = 1) +
  labs(
    x = "\nYearly Visit",
    y = "Axial Score\n"
  ) +
  scale_fill_manual(values = COLOR_PALETTE) +
  scale_colour_manual(values = COLOR_PALETTE) +
  apply_plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "axial_boxplot.svg", plot = plot_axial, width = 6, height = 6)

# Scatterplot: Asymmetry vs Axial Score
plot_axial_scatter <- axial_data %>%
  filter(abs(Diff) < ASYMMETRY_OUTLIER_THRESHOLD) %>%
  mutate(VISIT_label = convert_visit_to_label(VISIT)) %>%
  ggplot(aes(abs(Diff), Axial_Score)) +
  geom_jitter(width = 0.5, height = 0.5, alpha = 0.5, size = 2, shape = 1) +
  geom_smooth(method = "lm", colour = "#2841b0", fill = "#2841b0") +
  facet_wrap(~VISIT_label) +
  labs(
    x = "\nAbsolute Right-to-left Asymmetry",
    y = "Axial Score\n"
  ) +
  coord_cartesian(ylim = c(0, 5)) +
  apply_plot_theme()

ggsave(file = "axial_asymmetry_scatter.svg", plot = plot_axial_scatter, width = 12, height = 6)

# ============================================================================
# GAIT AND FREEZING ANALYSIS
# ============================================================================

axial_data <- axial_data %>%
  mutate(
    Gait_Impairment = as.factor(ifelse(ON_3.10_6 >= 2, 1, 0)),
    Freezing = as.factor(ifelse(ON_3.11_6 >= 1, 1, 0))
  )

# Gait analysis
for (visit in c(1, 3, 5)) {
  data_subset <- axial_data %>% filter(VISIT == visit)
  cat("\nGait - Visit", visit, "\n")
  print(wilcox.test(Diff ~ Gait_Impairment, data = data_subset))
}

# Freezing analysis
for (visit in c(1, 3, 5)) {
  data_subset <- axial_data %>% filter(VISIT == visit)
  cat("\nFreezing - Visit", visit, "\n")
  print(wilcox.test(Diff ~ Freezing, data = data_subset))
}

# ============================================================================
# PDQ39 QUALITY OF LIFE ANALYSIS
# ============================================================================

# Load PDQ39 data
pdq39_v0v1 <- read_xlsx(
  path = v0v1_asymmetry_path,
  sheet = "PDQ39-CGIS-SCOPA",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
) %>%
  select(SUBJID, starts_with("PDQ39")) %>%
  slice(-1) %>%
  group_by(SUBJID) %>%
  mutate(VISIT = 1) %>%
  ungroup() %>%
  select(SUBJID, VISIT, PDQ39_SCORE) %>%
  mutate(PDQ39_SCORE = as.numeric(PDQ39_SCORE)) %>%
  filter(PDQ39_SCORE != 0) %>%
  drop_na()

pdq39_v3v5 <- read_xlsx(
  path = v3v5_asymmetry_path,
  sheet = "PDQ39-CGIS-SCOPA _V3_V5",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
) %>%
  select(SUBJID, VISIT, PDQ39_SCORE) %>%
  slice(-1) %>%
  filter(VISIT != "Visit") %>%
  mutate(
    VISIT = ifelse(grepl("V3", VISIT), 3, 5),
    PDQ39_SCORE = as.numeric(PDQ39_SCORE)
  ) %>%
  drop_na() %>%
  filter(PDQ39_SCORE != 0)

pdq39_combined <- bind_rows(pdq39_v0v1, pdq39_v3v5)

# Merge with asymmetry and axial data
pdq39_data <- asymmetry_combined %>%
  inner_join(pdq39_combined) %>%
  inner_join(axial_combined %>% select(SUBJID, VISIT, Axial_Score)) %>%
  mutate(
    Total = Left + Right,
    Normalized_Asymmetry = abs(Diff) / Total,
    Diff = abs(Diff)
  )

# Filter subjects with complete data across all visits
pdq39_data <- pdq39_data %>%
  group_by(SUBJID) %>%
  filter(n() == 3) %>%
  ungroup()

# Statistical testing
pdq39_summary <- get_visit_summary_stats(pdq39_data, "PDQ39_SCORE")

friedman.test(
  y = pdq39_data$PDQ39_SCORE,
  groups = pdq39_data$VISIT,
  blocks = pdq39_data$SUBJID
)
pairwise.wilcox.test(
  pdq39_data$PDQ39_SCORE,
  pdq39_data$VISIT,
  p.adj = "bonferroni",
  paired = TRUE
)

# Correlation analysis
for (visit in c(1, 3, 5)) {
  data_subset <- pdq39_data %>% filter(VISIT == visit)
  cat("\nPDQ39 vs Asymmetry - Visit", visit, "\n")
  print(cor.test(
    abs(data_subset$Diff),
    data_subset$PDQ39_SCORE,
    method = "spearman"
  ))
}

# Linear regression
model_pdq39 <- lm(PDQ39_SCORE ~ Diff + Total, data = pdq39_data)
summary(model_pdq39)

# Mediation analysis
med_model_pdq39 <- lm(Total ~ Diff, data = pdq39_data)
outcome_model_pdq39 <- lm(PDQ39_SCORE ~ Diff + Total, data = pdq39_data)
med_results_pdq39 <- mediation::mediate(
  med_model_pdq39,
  outcome_model_pdq39,
  treat = "Diff",
  mediator = "Total",
  boot = TRUE,
  sims = 1000
)
summary(med_results_pdq39)

# Visualization
plot_pdq39 <- pdq39_data %>%
  mutate(VISIT_label = convert_visit_to_label(VISIT)) %>%
  ggplot(aes(VISIT_label, PDQ39_SCORE, colour = VISIT_label, fill = VISIT_label)) +
  geom_boxplot(
    alpha = 0.6,
    notch = TRUE,
    width = 0.5,
    outlier.colour = NULL,
    outlier.fill = NULL,
    outlier.alpha = 0.00001
  ) +
  geom_jitter(alpha = 0.8, height = 0.5, shape = 1) +
  labs(
    x = "\nYearly Visit",
    y = "PDQ39 Score\n"
  ) +
  scale_fill_manual(values = COLOR_PALETTE) +
  scale_colour_manual(values = COLOR_PALETTE) +
  apply_plot_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "pdq39_boxplot.svg", plot = plot_pdq39, width = 6, height = 6)

# ============================================================================
# PREDICTION MODEL: BASELINE PREDICTORS OF YEAR 5 ASYMMETRY
# ============================================================================

# Load demographic and baseline clinical data
demographics <- read_xlsx(
  path = v0v1_asymmetry_path,
  sheet = "DEMOGRAPHIE ",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
) %>%
  filter(SUBJID %in% subject_ids$SUBJID) %>%
  mutate(
    Age = as.numeric(AGE),
    Sex = ifelse(SEXE == "Homme", 1, 0),
    Disease_Duration = as.numeric(D_SCREEN) - as.numeric(D_1ER_SYMPT)
  ) %>%
  select(SUBJID, Age, Sex, Disease_Duration)

# Load baseline UPDRS III totals
updrs_totals <- read_xlsx(
  path = v0v1_asymmetry_path,
  sheet = "UPDRSIII_TOTAUX",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
) %>%
  filter(SUBJID %in% subject_ids$SUBJID) %>%
  mutate(UPDRS_Total_OFF = as.numeric(TOT_OFF_DRUG_V0)) %>%
  select(SUBJID, UPDRS_Total_OFF)

# Load baseline axial scores (OFF medication)
axial_baseline <- read_xlsx(
  path = v0v1_asymmetry_path,
  sheet = "UPDRSIII_COMPLET_V0_V1",
  skip = 0,
  col_types = "text",
  trim_ws = TRUE
) %>%
  select(SUBJID, OFF_3.9_, OFF_3.10_, OFF_3.11_, OFF_3.12_) %>%
  slice(-1) %>%
  mutate(across(-SUBJID, as.numeric)) %>%
  drop_na() %>%
  mutate(
    Axial_Score_Baseline = OFF_3.9_ + OFF_3.10_ + OFF_3.11_ + OFF_3.12_
  ) %>%
  select(SUBJID, Axial_Score_Baseline)

# Prepare prediction dataset
prediction_data <- asymmetry_combined %>%
  filter(VISIT == 1) %>%
  select(SUBJID) %>%
  inner_join(demographics) %>%
  left_join(updrs_totals) %>%
  left_join(axial_baseline) %>%
  drop_na()

# Merge with Year 5 asymmetry outcomes
prediction_data <- prediction_data %>%
  left_join(
    asymmetry_combined %>%
      filter(VISIT == 5) %>%
      select(SUBJID, Diff) %>%
      distinct()
  ) %>%
  drop_na()

# Fit model with standardized coefficients
model_predictors <- lm(
  Diff ~ Age + Sex + Disease_Duration + Axial_Score_Baseline + UPDRS_Total_OFF,
  data = prediction_data
)
summary(model_predictors)

# Standardize coefficients
prediction_data_std <- prediction_data %>%
  mutate(across(
    c(Age, Disease_Duration, Axial_Score_Baseline, UPDRS_Total_OFF),
    scale
  ))

model_predictors_std <- lm(
  Diff ~ Age + Sex + Disease_Duration + Axial_Score_Baseline + UPDRS_Total_OFF,
  data = prediction_data_std
)
summary(model_predictors_std)

# ============================================================================
# END OF ANALYSIS
# ============================================================================
