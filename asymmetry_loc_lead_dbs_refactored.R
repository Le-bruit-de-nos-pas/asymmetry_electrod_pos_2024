# =====================================================================
# PARKINSON'S DBS - ELECTRODE LOCATION AND VTA COVERAGE ANALYSIS
# =====================================================================
#
# Study: Analysis of deep brain stimulation electrode locations and
#        volume of tissue activated (VTA) coverage in relation to
#        motor asymmetry outcomes in Parkinson's disease patients
#
# Key Analyses:
#  - Electrode placement sites (STN, SM anatomical regions)
#  - VTA coverage for substantia nigra functional zones
#  - Correlation between VTA coverage and axial symptom improvement
#  - Clinical outcome groups (Asym->Sym vs Asym->Asym)
#  - Electrode coordinate analysis (X, Y, Z positions)
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


# =====================================================================
# CONFIGURATION CONSTANTS
# =====================================================================

# File paths
EXCEL_FILE <- "Asymmetry_DeepBrainStimulation.xlsx"
VTA_EXCEL_FILE <- "VTA_asymmetry.xlsx"
VTA_COVERAGE_FILE <- "VTA_patients_MB.xlsx"
ACTIVEPLOTS_FILE <- "activeplots.xlsx"
ASYMMETRY_DATA_FILE <- "Asymmetry_Pre_vs_Post.txt"
AXIAL_DATA_FILE <- "UPDRSI_II.txt"
PDQ39_DATA_FILE <- "PDQ39.txt"

# Analysis parameters
GROUP_THRESHOLD <- 22  # Minimum items for data inclusion
MISSING_DATA_THRESHOLD <- 13  # Missing items threshold
ASYMMETRY_THRESHOLD <- 5  # Clinical asymmetry cutoff

# Contact definitions
AXIAL_ITEMS <- c("3.9", "3.10", "3.11", "3.12")
SM_THRESHOLD <- 1  # Subthalamic margin indicator

# Group definitions
GROUP_LABELS <- list(
  Asym_to_Sym = "Asym_to_Sym(A)",
  Asym_to_Asym = "Aym_to_Asym(C)"
)

# Visualization constants
COLOR_PALETTE <- c(
  "Asym -> Sym" = "#0099E0",
  "Asym = Asym" = "#D45769"
)

SNC_ZONES <- list(
  somatomotor = "lateral SNc",
  limbic = "medial SNc",
  associative = "ventral SN"
)

STN_ZONES <- list(
  motor = "STN motor",
  associative = "STN associative",
  limbic = "STN limbic"
)

COORDINATE_AXES <- c("X", "Y", "Z")


# =====================================================================
# UTILITY FUNCTIONS
# =====================================================================

#' Load and prepare VTA asymmetry data
#'
#' @param vta_file Path to VTA asymmetry Excel file
#' @param active_plots_file Path to active plots file
#' @return Processed VTA data frame with standardized site coding
load_vta_data <- function(vta_file = VTA_EXCEL_FILE,
                          active_plots_file = ACTIVEPLOTS_FILE) {
  
  # Load electrode contact information
  activeplots <- read_xlsx(active_plots_file, col_types = "text", trim_ws = TRUE) %>%
    mutate(PLOT = ifelse(grepl("G", PLOT), "L", "R")) %>%
    mutate(CONTACT = paste0(CONTACT, PLOT)) %>%
    select(-PLOT)
  
  # Load and process VTA asymmetry data
  vta_data <- read_xlsx(vta_file, col_types = "text", trim_ws = TRUE) %>%
    rename(SUBJID = ID) %>%
    select(-c(group, `...3`, electrode)) %>%
    gather(feature, SITE, -SUBJID) %>%
    drop_na()
  
  # Standardize column naming
  vta_data <- vta_data %>%
    mutate(
      feature = str_replace_all(feature, "LH", "L"),
      feature = str_replace_all(feature, "RH", "R"),
      feature = str_replace_all(feature, "C0_", "0"),
      feature = str_replace_all(feature, "C1_", "1"),
      feature = str_replace_all(feature, "C2_", "2"),
      feature = str_replace_all(feature, "C3_", "3")
    )
  
  # Join with active electrode contact information
  vta_data <- vta_data %>%
    left_join(activeplots) %>%
    arrange(SUBJID) %>%
    filter(str_detect(feature, CONTACT))
  
  return(vta_data)
}


#' Extract and classify electrode sites
#'
#' Identifies STN vs SM sites and categorizes by functional/anatomical regions
#'
#' @param vta_data Processed VTA data
#' @return Data frame with site classifications
extract_electrode_sites <- function(vta_data) {
  
  sites <- vta_data %>%
    filter(feature %in% c("0L", "1L", "2L", "3L", "0R", "1R", "2R", "3R")) %>%
    select(SUBJID, feature, SITE) %>%
    distinct() %>%
    mutate(SITE = str_replace_all(SITE, ",", "/")) %>%
    separate_rows(SITE, sep = "/", convert = TRUE) %>%
    mutate(SITE = str_replace_all(SITE, " ", "")) %>%
    distinct()
  
  # Classify sites as STN_SM or Other
  sites <- sites %>%
    mutate(
      SITE = ifelse(SITE == "STN_SM", "STN_SM", "Other"),
      Contact_Side = str_sub(feature, 2L, 2L)
    ) %>%
    select(SUBJID, Contact_Side, SITE)
  
  return(sites)
}


#' Calculate axial motor scores from UPDRS III items
#'
#' Sums items 3.9-3.12 across specified motor states
#'
#' @param excel_file Path to UPDRSIII_COMPLET_V0_V1 sheet
#' @param states List of motor states to process (Pre-OP, OFFON, ONOFF, ONON, OFFOFF)
#' @return Data frame with axial scores by state
calculate_axial_scores <- function(excel_file = EXCEL_FILE, states = NULL) {
  
  UPDRSIII <- read_xlsx(
    excel_file,
    sheet = "UPDRSIII_COMPLET_V0_V1",
    col_types = "text",
    trim_ws = TRUE
  )
  
  df_names <- names(UPDRSIII)
  
  # Define column mappings for each state
  col_mappings <- list(
    PreOP = c("SUBJID", "OFF_3.9_", "OFF_3.10_", "OFF_3.11_", "OFF_3.12_"),
    ONOFF = c("SUBJID", "ONOFF_3.9_", "ONOFF_3.10_", "ONOFF_3.11_", "ONOFF_3.12_"),
    OFFON = c("SUBJID", "OFFON_3.9_", "OFFON_3.10_", "OFFON_3.11_", "OFFON_3.12_")
  )
  
  # Process each state
  axial_scores <- NULL
  
  for (state in names(col_mappings)) {
    cols <- col_mappings[[state]]
    state_data <- UPDRSIII[, cols] %>%
      filter(if_any(-SUBJID, ~ !is.na(.))) %>%
      mutate(across(-SUBJID, as.numeric)) %>%
      drop_na()
    
    # Calculate axial score (sum of 4 items)
    score_name <- paste0("AxialScore_", state)
    state_data[[score_name]] <- rowSums(state_data[, -1], na.rm = TRUE)
    
    state_data <- state_data %>% select(SUBJID, all_of(score_name))
    
    if (is.null(axial_scores)) {
      axial_scores <- state_data
    } else {
      axial_scores <- axial_scores %>% inner_join(state_data, by = "SUBJID")
    }
  }
  
  return(axial_scores)
}


#' Apply standardized visualization theme
#'
#' @param base_plot ggplot object
#' @return ggplot with applied theme
apply_density_theme <- function(base_plot) {
  base_plot +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "top",
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      axis.line = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 12, vjust = -0.5),
      axis.title.y = element_text(size = 12, vjust = -0.5),
      plot.margin = margin(5, 5, 5, 5, "pt")
    )
}


#' Compare outcome groups on quantitative variable
#'
#' Performs Wilcoxon test and returns summary statistics
#'
#' @param data Data frame
#' @param outcome Column name for outcome variable
#' @param group Column name for grouping variable
#' @param group_a First group identifier
#' @param group_b Second group identifier
#'
#' @return List with test results and summary stats
compare_groups_wilcox <- function(data, outcome, group,
                                  group_a, group_b) {
  
  values_a <- as.numeric(data[[outcome]][data[[group]] == group_a])
  values_b <- as.numeric(data[[outcome]][data[[group]] == group_b])
  
  test <- wilcox.test(values_a, values_b)
  
  result <- list(
    Group_A_Mean = mean(values_a, na.rm = TRUE),
    Group_A_SD = sd(values_a, na.rm = TRUE),
    Group_B_Mean = mean(values_b, na.rm = TRUE),
    Group_B_SD = sd(values_b, na.rm = TRUE),
    W = test$statistic,
    P_Value = test$p.value
  )
  
  return(result)
}


#' Create comparison visualization for two outcome groups
#'
#' @param data Data frame with outcome and grouping variables
#' @param x_var Variable for x-axis (quantitative)
#' @param group_var Grouping variable (qualitative)
#' @param title Plot title
#' @param x_label X-axis label
#'
#' @return ggplot object
create_group_comparison_plot <- function(data, x_var, group_var,
                                         title, x_label) {
  
  plot <- data %>%
    mutate(!!sym(group_var) := ifelse(
      !!sym(group_var) == GROUP_LABELS$Asym_to_Sym,
      "Asym -> Sym", "Asym = Asym"
    )) %>%
    ggplot(aes(x = !!sym(x_var), 
               colour = !!sym(group_var), 
               fill = !!sym(group_var))) +
    geom_density(linewidth = 2, alpha = 0.25) +
    facet_wrap(~ as.factor(!!sym(group_var)), ncol = 2) +
    labs(
      title = title,
      x = x_label,
      y = "Patient Density\nGaussian Kernel Estimate"
    ) +
    apply_density_theme() +
    scale_fill_manual(values = COLOR_PALETTE) +
    scale_colour_manual(values = COLOR_PALETTE)
  
  return(plot)
}


#' Perform correlation analysis with multiple variables
#'
#' @param data Data frame
#' @param outcome Numeric outcome variable
#' @param predictor_cols Vector of predictor column names
#' @param method Correlation method ("pearson", "spearman")
#'
#' @return Data frame with correlation results
correlation_batch <- function(data, outcome, predictor_cols, method = "spearman") {
  
  results <- data.frame(
    predictor = character(),
    correlation = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  outcome_vals <- data[[outcome]]
  
  for (pred in predictor_cols) {
    pred_vals <- data[[pred]]
    test <- cor.test(outcome_vals, pred_vals, method = method)
    
    results <- rbind(results, data.frame(
      predictor = pred,
      correlation = test$estimate,
      p_value = test$p.value,
      stringsAsFactors = FALSE
    ))
  }
  
  return(results)
}


# =====================================================================
# SECTION 1: LOAD AND PREPARE DATA
# =====================================================================

# Load VTA and electrode data
vta_data <- load_vta_data()
groups <- read_xlsx(VTA_EXCEL_FILE, col_types = "text", trim_ws = TRUE) %>%
  select(ID, group) %>%
  rename(SUBJID = ID)

# Extract electrode sites
sites <- extract_electrode_sites(vta_data)

# Load asymmetry outcomes
asymmetry_data <- fread(ASYMMETRY_DATA_FILE)

# Prepare analysis cohort
analysis_cohort <- sites %>%
  select(SUBJID) %>%
  distinct() %>%
  inner_join(groups)


# =====================================================================
# SECTION 2: ELECTRODE SITE CLASSIFICATION ANALYSIS
# =====================================================================

# Classify patients by site location
site_classification <- sites %>%
  mutate(Site_Indicator = 1) %>%
  spread(key = SITE, value = Site_Indicator, fill = 0) %>%
  mutate(
    Has_SM = ifelse(STN_SM == 1, 1, 0),
    Site_Type_Left = Contact_Side == "L" & Has_SM,
    Site_Type_Right = Contact_Side == "R" & Has_SM
  ) %>%
  select(-c(STN_SM, Other, Contact_Side)) %>%
  inner_join(groups)

# Summary: STN_SM location prevalence
site_summary <- site_classification %>%
  group_by(group) %>%
  summarise(
    N_With_SM = sum(Has_SM),
    N_Total = n(),
    Percent_SM = round(100 * sum(Has_SM) / n(), 1)
  )

print("STN_SM Site Distribution by Group:")
print(site_summary)


# =====================================================================
# SECTION 3: MOTOR ASYMMETRY BY ELECTRODE LOCATION
# =====================================================================

# Prepare asymmetry data with site information
asym_by_location <- sites %>%
  mutate(
    Site_Binary = 1,
    Contact_Label = str_sub(feature, 2L, 2L)
  ) %>%
  spread(key = SITE, value = Site_Binary, fill = 0) %>%
  mutate(
    SM_Indicator = ifelse(STN_SM == 1, 1, 0)
  ) %>%
  select(-c(Other, STN_SM, feature)) %>%
  inner_join(groups) %>%
  spread(key = Contact_Label, value = SM_Indicator, fill = 0) %>%
  mutate(
    Left_SM = ifelse(is.na(L), 0, L),
    Right_SM = ifelse(is.na(R), 0, R),
    Either_Side_SM = ifelse(Left_SM == 1 | Right_SM == 1, 1, 0),
    Both_Sides_SM = ifelse(Left_SM == 1 & Right_SM == 1, 1, 0)
  ) %>%
  select(SUBJID, group, Left_SM, Right_SM, Either_Side_SM, Both_Sides_SM) %>%
  inner_join(asymmetry_data)

# Statistical tests: Asymmetry by SM location
sm_comparisons <- list(
  Both_Sites = compare_groups_wilcox(
    asym_by_location, "ON_ON", "group",
    GROUP_LABELS$Asym_to_Sym, GROUP_LABELS$Asym_to_Asym
  ),
  Either_Site = compare_groups_wilcox(
    asym_by_location, "ON_ON", "group",
    GROUP_LABELS$Asym_to_Sym, GROUP_LABELS$Asym_to_Asym
  )
)

print("Motor Asymmetry Comparison by SM Location:")
print(sm_comparisons)


# =====================================================================
# SECTION 4: AXIAL MOTOR SCORE ANALYSIS
# =====================================================================

# Calculate axial scores
axial_scores <- calculate_axial_scores()

# Prepare axial data with group classification
axial_by_group <- axial_scores %>%
  inner_join(asym_by_location) %>%
  rename(
    Pre_OP = AxialScore_PreOP,
    ON_DBS_OFF_Med = AxialScore_ONOFF,
    OFF_DBS_ON_Med = AxialScore_OFFON
  ) %>%
  mutate(Delta_Axial = ON_DBS_OFF_Med - Pre_OP)

# Statistical tests: Axial improvement by outcome group
axial_tests <- list(
  PreOP = compare_groups_wilcox(
    axial_by_group, "Pre_OP", "group",
    GROUP_LABELS$Asym_to_Sym, GROUP_LABELS$Asym_to_Asym
  ),
  PostOP = compare_groups_wilcox(
    axial_by_group, "ON_DBS_OFF_Med", "group",
    GROUP_LABELS$Asym_to_Sym, GROUP_LABELS$Asym_to_Asym
  ),
  Delta = compare_groups_wilcox(
    axial_by_group, "Delta_Axial", "group",
    GROUP_LABELS$Asym_to_Sym, GROUP_LABELS$Asym_to_Asym
  )
)

print("Axial Score Comparison by Outcome Group:")
print(axial_tests)


# =====================================================================
# SECTION 5: ELECTRODE COORDINATE ANALYSIS
# =====================================================================

# Extract X, Y, Z coordinates
coordinates <- vta_data %>%
  filter(grepl("X", feature) | grepl("Y", feature) | grepl("Z", feature)) %>%
  mutate(
    Contact_Side = str_sub(CONTACT, 2L, 2L),
    Axis = str_sub(feature, 4L, 4L),
    Coordinate = as.numeric(SITE)
  ) %>%
  select(SUBJID, Contact_Side, Axis, Coordinate) %>%
  inner_join(groups)

# Coordinate comparison by outcome group
coordinate_results <- data.frame(
  Axis = character(),
  Side = character(),
  W = numeric(),
  P_Value = numeric(),
  Group_A_Mean = numeric(),
  Group_A_SD = numeric(),
  Group_B_Mean = numeric(),
  Group_B_SD = numeric(),
  stringsAsFactors = FALSE
)

for (axis in COORDINATE_AXES) {
  for (side in c("L", "R")) {
    subset_data <- coordinates %>%
      filter(Axis == axis, Contact_Side == side)
    
    if (nrow(subset_data) > 0) {
      test_result <- compare_groups_wilcox(
        subset_data, "Coordinate", "group",
        GROUP_LABELS$Asym_to_Sym, GROUP_LABELS$Asym_to_Asym
      )
      
      coordinate_results <- rbind(coordinate_results, data.frame(
        Axis = axis,
        Side = side,
        W = test_result$W,
        P_Value = test_result$P_Value,
        Group_A_Mean = test_result$Group_A_Mean,
        Group_A_SD = test_result$Group_A_SD,
        Group_B_Mean = test_result$Group_B_Mean,
        Group_B_SD = test_result$Group_B_SD,
        stringsAsFactors = FALSE
      ))
    }
  }
}

print("Electrode Coordinate Analysis:")
print(coordinate_results)


# =====================================================================
# SECTION 6: VTA COVERAGE ANALYSIS - SUBSTANTIA NIGRA
# =====================================================================

# Load VTA coverage data
vta_coverage <- read_xlsx(VTA_COVERAGE_FILE, trim_ws = TRUE) %>%
  filter(Group == "B") %>%
  rename(SUBJID = ID) %>%
  select(
    SUBJID,
    contains("Zhang_Right_lateral_SNc"),
    contains("Zhang_Right_medial_SNc"),
    contains("Zhang_Right_ventral_SN"),
    contains("Zhang_Left_lateral_SNc"),
    contains("Zhang_Left_medial_SNc"),
    contains("Zhang_Left_ventral_SN")
  )

# VTA vs axial outcome
axial_delta <- axial_by_group %>%
  select(SUBJID, Delta_Axial) %>%
  inner_join(vta_coverage)

# Correlation analysis: VTA coverage and axial improvement
sn_predictors <- colnames(vta_coverage)[-1]
vta_snc_correlations <- correlation_batch(
  axial_delta,
  "Delta_Axial",
  sn_predictors,
  method = "spearman"
)

print("VTA Coverage (SNc) vs Axial Motor Improvement:")
print(vta_snc_correlations)


# =====================================================================
# SECTION 7: VTA COVERAGE ANALYSIS - SUBTHALAMIC NUCLEUS
# =====================================================================

# Load STN-specific VTA coverage
vta_stn <- read_xlsx(VTA_COVERAGE_FILE, trim_ws = TRUE) %>%
  select(
    ID,
    Group,
    contains("Ewert_Right_STN"),
    contains("Ewert_Left_STN")
  ) %>%
  rename(SUBJID = ID)

# Group comparison for STN coverage
stn_coverage_results <- data.frame(
  Region = character(),
  Side = character(),
  W = numeric(),
  P_Value = numeric(),
  Group_A_Mean = numeric(),
  Group_A_SD = numeric(),
  Group_B_Mean = numeric(),
  Group_B_SD = numeric(),
  stringsAsFactors = FALSE
)

stn_cols <- colnames(vta_stn)[-c(1, 2)]

for (col in stn_cols) {
  test_result <- compare_groups_wilcox(
    vta_stn, col, "Group",
    "A", "B"
  )
  
  stn_coverage_results <- rbind(stn_coverage_results, data.frame(
    Region = col,
    W = test_result$W,
    P_Value = test_result$P_Value,
    Group_A_Mean = test_result$Group_A_Mean,
    Group_A_SD = test_result$Group_A_SD,
    Group_B_Mean = test_result$Group_B_Mean,
    Group_B_SD = test_result$Group_B_SD,
    stringsAsFactors = FALSE
  ))
}

print("STN VTA Coverage Comparison:")
print(stn_coverage_results)


# =====================================================================
# SECTION 8: CLINICAL OUTCOMES SUMMARY
# =====================================================================

# Comprehensive summary by outcome group
clinical_summary <- analysis_cohort %>%
  left_join(site_summary %>% select(group, N_With_SM, N_Total)) %>%
  left_join(
    axial_by_group %>%
      group_by(group) %>%
      summarise(
        Mean_Delta_Axial = mean(Delta_Axial),
        SD_Delta_Axial = sd(Delta_Axial),
        N = n()
      )
  ) %>%
  distinct()

print("Clinical Outcome Summary:")
print(clinical_summary)

# =====================================================================
# END OF ANALYSIS
# =====================================================================
