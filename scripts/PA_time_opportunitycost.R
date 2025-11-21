# ============================================================
# Title: Spatial Covariate Extraction and Spearman Correlation Analysis
# Purpose:
#   Assess whether newer PAs are established in areas with higher pressures
# ============================================================

# -----------------------------
# 0. Setup
# -----------------------------
rm(list = ls())  # Clear workspace

# Set working directory

# Load required packages
library(dplyr)
library(ggplot2)

# -----------------------------
# 1. Load and Prepare Data
# -----------------------------
gridcell_sample <- read.csv("df_sampled_gridcells_with_covariates_Mar25.csv") %>%
  mutate(type = as.factor(type))  # 0=random, 1=PA, 2=buffer

# Covariates of interest
covariates <- c("Population", "Agri_03", "Access", "Elevation")

# Remove rows with NA in any covariate
gridcell_sample <- gridcell_sample %>%
  filter(if_all(all_of(covariates), ~ !is.na(.x)))

# -----------------------------
# 2. Analysis A: Within Protected Areas (Spearman correlations)
# -----------------------------
pa_df <- gridcell_sample %>%
  filter(type == "1", !is.na(Earliest_year), Earliest_year != 1600)

pa_df_thin <- pa_df %>%
  group_by(Name_PAs) %>%                     # preserve type balance
  sample_frac(0.02) %>%                        # keep 10% of each type
  ungroup()

table(pa_df$Name_PAs)
table(pa_df_thin$Name_PAs)

# Function to compute Spearman correlation per covariate
compute_spearman_info <- function(var, data) {
  x <- data$Earliest_year
  y <- data[[var]]
  test <- suppressWarnings(cor.test(x, y, method = "spearman"))
  
  # Significance label
  sig_label <- if (test$p.value < 0.001) "***" else if (test$p.value < 0.01) "**" else if (test$p.value < 0.05) "*" else "NS"
  
  # Direction
  direction <- ifelse(test$estimate > 0, "Increases with year", "Decreases with year")
  
  data.frame(
    Covariate = var,
    Spearman_rho = unname(test$estimate),
    P_value = test$p.value,
    Direction = direction,
    Signif = sig_label
  )
}

# Apply to all covariates
spearman_results <- lapply(covariates, compute_spearman_info, data = pa_df_thin) %>% bind_rows()

print(spearman_results)

