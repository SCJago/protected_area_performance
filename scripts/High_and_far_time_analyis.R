# ============================================================
# Title: Spatial Covariate Extraction and Spearman Correlation Analysis
# Purpose:
#   1. Compare environmental and socioeconomic covariates within Protected Areas (PAs)
#   2. Assess whether newer PAs are established in areas with higher pressures
#   3. Compare covariates between PAs and national random sites
# ============================================================

# -----------------------------
# 0. Setup
# -----------------------------
rm(list = ls())  # Clear workspace

# Set working directory
setwd("~/Papers/Ethiopia - PA history review/Write up - reviewing/Nat_Eco_Evo_Submission/Review/Ecoregion_additional_analysis")

# Load required packages
library(dplyr)
library(ggplot2)
library(tidyr)

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


# Custom labels for plots
covariate_labels <- c(
  Access = "Travel time to city (minutes)",
  Agri_03 = "Agricultural land (percentage)",
  Elevation = "Elevation (m)",
  Population = "Population"
)

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


#check autocorrelation
# Set spatial coordinates once
library(sp)
coordinates(pa_df_thin) <- ~ x_coord + y_coord


cutoff <- 30000     # keep only pairs within this distance
width  <- cutoff / 20          # ~12 lag bins within the cutoff

# If you prefer a hard distance (e.g., 50 km), just do:
# cutoff <- 50000; width <- cutoff / 12

# 3) Variograms using that cutoff (no per-PA splitting)
vgm_pop  <- variogram(Population ~ 1, data = pa_df_thin, cutoff = cutoff, width = width)
vgm_agri <- variogram(Agri_03    ~ 1, data = pa_df_thin, cutoff = cutoff, width = width)
vgm_acc  <- variogram(Access     ~ 1, data = pa_df_thin, cutoff = cutoff, width = width)
vgm_elev <- variogram(Elevation  ~ 1, data = pa_df_thin, cutoff = cutoff, width = width)

# 4) Quick base plots (same style as your residuals example)
par(mfrow = c(2, 2))
plot(vgm_pop,  main = "Semivariogram: Population (nearby pairs)")
plot(vgm_agri, main = "Semivariogram: Agri_03 (nearby pairs)")
plot(vgm_acc,  main = "Semivariogram: Access (nearby pairs)")
plot(vgm_elev, main = "Semivariogram: Elevation (nearby pairs)")
par(mfrow = c(1, 1))


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

# -----------------------------
# 2a. Convert to long format for plotting
# -----------------------------
pa_df_thin <- st_as_sf(pa_df_thin)    

pa_long <- pa_df_thin %>% 
  st_drop_geometry() %>%    
  select(Earliest_year, all_of(covariates)) %>%
  pivot_longer(
    cols = all_of(covariates),
    names_to = "Covariate",
    values_to = "Value"
  )



# -----------------------------
# 2b. Plot covariate trends over time (scatter with smoothed line)
# -----------------------------
ggplot(pa_long, aes(x = Earliest_year, y = Value)) +
  geom_smooth(method = "loess", color = "black", fill = "grey70") +
  facet_wrap(~Covariate, scales = "free_y", labeller = labeller(Covariate = covariate_labels)) +
  scale_x_continuous(
    breaks = seq(min(pa_df$Earliest_year), max(pa_df$Earliest_year), by = 20),
    minor_breaks = seq(min(pa_df$Earliest_year), max(pa_df$Earliest_year), by = 10)
  ) +
  labs(x = "Year of establishment", y = "") +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0)
  )

# -----------------------------
# 3. Analysis B: Compare PAs vs National Random Sites (unchanged)
# -----------------------------
set.seed(123)  # for reproducibility
gridcell_sample_thin <- gridcell_sample %>%
  group_by(Ecoregion) %>%                     # preserve type balance
  sample_frac(0.1) %>%                        # keep 10% of each type
  ungroup()

gridcell_sample_no_buffer <- gridcell_sample_thin %>%
  filter(type != "2") %>%
  mutate(type = droplevels(type))

response_matrix_bg <- as.matrix(gridcell_sample_no_buffer[, covariates])
manova_bg <- manova(response_matrix_bg ~ type, data = gridcell_sample_no_buffer)

cat("\nMANOVA: PAs vs National Random (type)\n")
print(summary(manova_bg, test = "Pillai"))
print(summary.aov(manova_bg))

# -----------------------------
# 4. Visualize Covariate Distributions by Site Type (unchanged)
# -----------------------------
for (v in covariates) {
  ggplot(gridcell_sample_no_buffer, aes(x = type, y = .data[[v]], fill = type)) +
    geom_boxplot(alpha = 0.6) +
    labs(title = paste("Distribution of", v, "by site type"),
         x = "Site type (0 = random, 1 = PA)", y = v) +
    theme_classic(base_size = 13) +
    theme(legend.position = "none") -> p
  print(p)
}

# -----------------------------
# 5. Summarize Group Means and SDs (unchanged)
# -----------------------------
group_means <- gridcell_sample_no_buffer %>%
  group_by(type) %>%
  summarise(across(all_of(covariates),
                   list(mean = \(x) mean(x, na.rm = TRUE),
                        sd   = \(x) sd(x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

print(group_means)

