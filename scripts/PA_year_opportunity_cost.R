# ======================================
# Full workflow: Opportunity Cost of PAs (pop + ag land + ag suitability)
# Decade aggregation + linear & GAM trends
# ======================================

# -----------------------------
# Load libraries
# -----------------------------
library(sf)
library(dplyr)
library(terra)
library(ggplot2)
library(viridis)
library(mgcv)

# -----------------------------
# Set working directory
# -----------------------------
setwd("~/Papers/Ethiopia - PA history review/Write up - reviewing/Nat_Eco_Evo_Submission/Review/Ecoregion_additional_analysis")

# -----------------------------
# Load PA shapefile and year data
# -----------------------------
PA_shp <- st_read("Eth_PAs_Sep24.shp")
PA_yr  <- read.csv("Earliest_PA_year.csv")
names(PA_yr)[names(PA_yr) == "Name_PAs"] <- "NAME_PAs"

# Fix invalid geometries
PA_shp <- st_make_valid(PA_shp)

# Join year info
PA_shp_joined <- PA_shp %>%
  left_join(PA_yr, by = "NAME_PAs")

# -----------------------------
# Load rasters
# -----------------------------
Ag_suit   <- rast("agri_suit_rain_clipped_1980-2009.tif")
Pop_2000  <- rast("eth_pd_2000_1km_UNadj.tif")
Ag_land   <- rast("covariate_baseline_agri_2003.tif")
Access   <- rast("updated_access_2000.tif")
Elev     <- rast("Elev.tif")

# Convert sf to SpatVector
PA_vect <- vect(PA_shp)

plot(PA_vect)

# -----------------------------
# Extract raster values per PA
# -----------------------------
vals_ag      <- terra::extract(Ag_suit, PA_vect)
vals_pop     <- terra::extract(Pop_2000, PA_vect)
vals_agland  <- terra::extract(Ag_land, PA_vect)
vals_elev    <- terra::extract(Elev, PA_vect)
vals_access    <- terra::extract(Access, PA_vect)
   

# -----------------------------
# Summarise raster values per PA (mean values)
# -----------------------------
vals_ag_summary <- vals_ag %>%
  mutate(NAME_PAs = PA_shp$NAME_PAs[ID]) %>%
  group_by(NAME_PAs) %>%
  summarise(mean_suit = mean(`agri_suit_rain_clipped_1980-2009`, na.rm = TRUE))

vals_pop_summary <- vals_pop %>%
  mutate(NAME_PAs = PA_shp$NAME_PAs[ID]) %>%
  group_by(NAME_PAs) %>%
  summarise(mean_pop = mean(`eth_pd_2000_1km_UNadj`, na.rm = TRUE))

vals_agland_summary <- vals_agland %>%
  mutate(NAME_PAs = PA_shp$NAME_PAs[ID]) %>%
  group_by(NAME_PAs) %>%
  summarise(mean_agland = mean(`covariate_baseline_agri_2003`, na.rm = TRUE))

vals_elev_summary <- vals_elev %>%
  mutate(NAME_PAs = PA_shp$NAME_PAs[ID]) %>%
  group_by(NAME_PAs) %>%
  summarise(mean_elev = mean(`Elev`, na.rm = TRUE))

vals_access_summary <- vals_access %>%
  mutate(NAME_PAs = PA_shp$NAME_PAs[ID]) %>%
  group_by(NAME_PAs) %>%
  summarise(mean_access = mean(`updated_access_2000`, na.rm = TRUE))

# -----------------------------
# Combine all data
# -----------------------------
PA_summary <- PA_yr %>%
  left_join(vals_ag_summary, by = "NAME_PAs") %>%
  left_join(vals_pop_summary, by = "NAME_PAs") %>%
  left_join(vals_agland_summary, by = "NAME_PAs") %>%
  left_join(vals_elev_summary, by = "NAME_PAs") %>%
  left_join(vals_access_summary, by = "NAME_PAs") %>%
  filter(Earliest_year != 1600)  # remove extreme outlier

# -----------------------------
# Compute opportunity cost (scaled)
# -----------------------------
PA_summary <- PA_summary %>%
  mutate(
    suit_scaled   = scale(mean_suit),
    pop_scaled    = scale(mean_pop),
    agland_scaled = scale(mean_agland),
    elev_scaled = scale(mean_elev),
    access_scaled = scale(mean_access),
    opportunity_cost = pop_scaled + agland_scaled
  )

colnames(PA_summary)

# List of response variables
response_vars <- c("mean_suit", "mean_pop", "mean_agland", "mean_elev", "mean_access")

# Create an empty list to store results
lm_results <- list()

# Loop through each response variable
for (var in response_vars) {
  # Create formula dynamically
  formula <- as.formula(paste(var, "~ Earliest_year"))
  
  # Fit linear model
  lm_model <- lm(formula, data = PA_summary)
  
  # Store summary in the list
  lm_results[[var]] <- summary(lm_model)
  
  # Print results
  cat("\n==============================\n")
  cat("Linear regression for:", var, "\n")
  print(lm_results[[var]])
}

# Fit linear model
lm_earliest <- lm(Earliest_year ~ suit_scaled + pop_scaled + agland_scaled + elev_scaled + access_scaled,
                  data = PA_summary)
options(na.action = "na.fail")

model_set <- dredge(lm_earliest)

# Extract the best model (lowest AICc)
best_model <- get.models(model_set, 1)[[1]]

# Summary of the best model
summary(best_model)
# -----------------------------
# Linear model & GAM
# -----------------------------
model_oc <- lm(opportunity_cost ~ Earliest_year, data = PA_summary)
summary(model_oc)

gam_oc <- gam(opportunity_cost ~ s(Earliest_year), data = PA_summary)
summary(gam_oc)

# -----------------------------
# Scatter plot with linear & GAM trends
# -----------------------------
ggplot(PA_summary, aes(x = Earliest_year, y = opportunity_cost)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", color = "darkgreen", se = TRUE) +
  #geom_line(aes(y = opportunity_cost_gam), color = "darkred", size = 1.2) +
  labs(
    x = "Year of PA establishment",
    y = "Opportunity cost (z-score)",
    title = "Opportunity costs of newly established PAs"
  ) +
  theme_minimal()

# -----------------------------
# Aggregate by decade
# -----------------------------
PA_summary_decade <- PA_summary %>%
  mutate(Decade = floor(Earliest_year / 10) * 10) %>%
  group_by(Decade) %>%
  summarise(
    mean_oc = mean(opportunity_cost, na.rm = TRUE),
    se_oc   = sd(opportunity_cost, na.rm = TRUE) / sqrt(n())
  )

# Linear model on decade-aggregated data
model_oc_decade <- lm(mean_oc ~ Decade, data = PA_summary_decade)
summary(model_oc_decade)

# -----------------------------
# Plot decade-aggregated opportunity cost
# -----------------------------
ggplot(PA_summary_decade, aes(x = Decade, y = mean_oc)) +
  geom_point(size = 3) +
  geom_line(color = "darkgreen", size = 1) +
  geom_errorbar(aes(ymin = mean_oc - se_oc, ymax = mean_oc + se_oc),
                width = 3, color = "darkgreen") +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  scale_x_continuous(breaks = seq(min(PA_summary_decade$Decade),
                                  max(PA_summary_decade$Decade), by = 10)) +
  
  labs(
    x = "Decade of PA establishment",
    y = "Mean opportunity cost (z-score)",
    title = "Decade-aggregated opportunity cost of PAs"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -----------------------------
# Map opportunity cost
# -----------------------------
PA_shp_plot <- PA_shp %>%
  left_join(PA_summary %>% select(NAME_PAs, opportunity_cost), by = "NAME_PAs") %>%
  mutate(geometry = st_make_valid(geometry)) %>%
  filter(st_is_valid(geometry))

ggplot(PA_shp_plot) +
  geom_sf(aes(fill = opportunity_cost), color = "black", size = 0.1) +
  scale_fill_viridis(option = "plasma", na.value = "grey80") +
  labs(title = "Opportunity cost of PAs (pop + ag land + ag suitability)") +
  theme_minimal()



# ======================================
# Test if PAs are located in more 'high and far' locations
# ======================================

# -----------------------------
# Load libraries
# -----------------------------
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)



# -----------------------------
# Get Ethiopia boundary locally from Natural Earth
# -----------------------------
ethiopia_boundary <- ne_countries(country = "Ethiopia", returnclass = "sf")
ethiopia_boundary <- st_make_valid(ethiopia_boundary)

# -----------------------------
# Load rasters
# -----------------------------
Ag_suit  <- rast("agri_suit_rain_clipped_1980-2009.tif")
Pop_2000 <- rast("eth_pd_2000_1km_UNadj.tif")
Ag_land  <- rast("covariate_baseline_agri_2003.tif")
Access   <- rast("updated_access_2000.tif")
Elev     <- rast("Elev.tif")

# -----------------------------
# Generate random points
# -----------------------------
set.seed(42)


# Rasterize polygons for sampling
r <- rast(PA_vect, res = 0.01)  # adjust resolution as needed
r[] <- NA
rasterized <- rasterize(PA_vect, r, field = "NAME_PAs")

# Extract coordinates inside polygons
coords <- xyFromCell(rasterized, which(!is.na(values(rasterized))))

# Sample 5000 points across all PAs
set.seed(42)
sample_idx <- sample(1:nrow(coords), size = 5000, replace = TRUE)
sample_points <- coords[sample_idx, ]

# Convert to sf points
rand_in_PA <- st_as_sf(vect(sample_points, type = "points"))
# Attach PA names
rand_in_PA$NAME_PAs <- extract(rasterized, vect(sample_points))[,2]
rand_in_PA <- st_as_sf(rand_in_PA)          # from your spatSample
st_crs(rand_in_PA) <- st_crs(PA_shp)            # assign same CRS as PA polygons
rand_in_PA$Source <- "Inside PA"
plot(rand_in_PA)

# -----------------------------
# 2. 5000 national random points
# -----------------------------
rand_nat <- st_sample(ethiopia_boundary, size = 5000, type = "random")
rand_nat <- st_sf(geometry = rand_nat)
rand_nat$Source <- "National Random"

# -----------------------------
# 3. Combine points
# -----------------------------
# Make sure both are in the same CRS
rand_in_PA <- st_transform(rand_in_PA, st_crs(rand_nat))
# Keep only geometry + Source
rand_in_PA <- rand_in_PA %>% select(Source, geometry)
rand_nat   <- rand_nat %>% select(Source, geometry)
all_points <- rbind(rand_in_PA, rand_nat)

# -----------------------------
# 4. Project to UTM for buffering
# -----------------------------
utm_crs <- 32637  # Ethiopia UTM zone
all_points_proj <- st_transform(all_points, crs = utm_crs)

# Buffer around points (1 km radius)
buffer_radius <- 1000
all_buffers <- st_buffer(all_points_proj, dist = buffer_radius)

# -----------------------------
# 5. Align rasters to reference raster
# -----------------------------
ref <- Access
Pop_2000 <- project(Pop_2000, ref)
Ag_land  <- project(Ag_land, ref)
Ag_suit  <- project(Ag_suit, ref)
Elev     <- project(Elev, ref)

cov_stack <- c(Ag_suit, Pop_2000, Ag_land, Access, Elev)
names(cov_stack) <- c("mean_suit", "mean_pop", "mean_agland", "mean_access", "mean_elev")

# -----------------------------
# 6. Extract mean raster values within buffers
# -----------------------------
all_buffers_vect <- vect(all_buffers)
vals_mean <- terra::extract(cov_stack, all_buffers_vect, fun = mean, na.rm = TRUE)

# Combine with buffer attributes
rand_df_buffer <- cbind(all_buffers, vals_mean[, -1])
rand_df_buffer_clean <- rand_df_buffer %>%
  filter(complete.cases(st_drop_geometry(.)))

# -----------------------------
# 7. T-tests for covariates
# -----------------------------
cat("\n==== Buffered point t-tests (Inside PA vs National Random) ====\n")
covariates <- c("mean_suit", "mean_pop", "mean_agland", "mean_access", "mean_elev")
for (v in covariates) {
  test <- t.test(rand_df_buffer_clean[[v]] ~ rand_df_buffer_clean$Source)
  cat("\n---", v, "---\n")
  print(test)
}
# -----------------------------
# Optional: logistic regression (multivariate)
# -----------------------------
 rand_df_buffer_clean$Source_bin <- ifelse(rand_df_buffer_clean$Source == "Inside PA", 1, 0)
 glm_model <- glm(Source_bin ~ mean_suit + mean_pop + mean_agland + mean_access + mean_elev,
                  data = rand_df_buffer_clean, family = binomial)
 summary(glm_model)

# -----------------------------
# Optional: map of points and buffers
# -----------------------------
 ggplot() +
   geom_sf(data = ethiopia_boundary, fill = "grey95", color = "black") +
      geom_sf(data = PA_shp, fill = NA, color = "darkgreen", size = 0.5) +
   geom_sf(data = all_buffers, aes(color = Source), alpha = 0.3, size = 0.5) +
   scale_color_viridis_d(option = "plasma") +
   labs(title = "Buffered random points inside PAs vs national") +
   theme_minimal()
 
 
 
 
 
 ##########################################################################################
 #do above analysis with YEAR
 ##########################################################################################
 # -----------------------------
 # 1. Load Ethiopia boundary
 # -----------------------------
 ethiopia_boundary <- ne_countries(country = "Ethiopia", returnclass = "sf")
 ethiopia_boundary <- st_make_valid(ethiopia_boundary)
 
 # -----------------------------
 # Load rasters
 # -----------------------------
 Ag_suit  <- rast("agri_suit_rain_clipped_1980-2009.tif")
 Pop_2000 <- rast("eth_pd_2000_1km_UNadj.tif")
 Ag_land  <- rast("covariate_baseline_agri_2003.tif")
 Access   <- rast("updated_access_2000.tif")
 Elev     <- rast("Elev.tif")
 
 # -----------------------------
 # Generate random points inside PAs
 # -----------------------------
 set.seed(42)
 
 # Rasterize polygons for sampling
 r <- rast(PA_vect, res = 0.01)  # adjust resolution
 r[] <- NA
 rasterized <- rasterize(PA_vect, r, field = "NAME_PAs")
 
 # Extract coordinates inside polygons
 coords <- xyFromCell(rasterized, which(!is.na(values(rasterized))))
 
 # Sample 5000 points
 sample_idx <- sample(1:nrow(coords), size = 5000, replace = TRUE)
 sample_points <- coords[sample_idx, ]
 
 # Convert to sf points
 rand_in_PA <- st_as_sf(vect(sample_points, type = "points"))
 
 # Attach PA names
 rand_in_PA$NAME_PAs <- extract(rasterized, vect(sample_points))[,2]
 
 # Attach Earliest_year from PA_yr
 rand_in_PA <- rand_in_PA %>%
   left_join(PA_yr %>% select(NAME_PAs, Earliest_year), by = "NAME_PAs")
 
 rand_in_PA <- st_as_sf(rand_in_PA)
 st_crs(rand_in_PA) <- st_crs(PA_shp)  # assign same CRS as PA polygons
 rand_in_PA$Source <- "Inside PA"
 plot(rand_in_PA)
 
 # -----------------------------
 # 2. 5000 national random points
 # -----------------------------
 rand_nat <- st_sample(ethiopia_boundary, size = 5000, type = "random")
 rand_nat <- st_sf(geometry = rand_nat)
 rand_nat$Source <- "National Random"
 
 # -----------------------------
 # 3. Combine points
 # -----------------------------
 rand_in_PA <- st_transform(rand_in_PA, st_crs(rand_nat))
 rand_in_PA <- rand_in_PA %>% select(Source, Earliest_year, geometry)
 rand_nat   <- rand_nat %>% select(Source, geometry)
 all_points <- rbind(rand_in_PA, rand_nat)
 
 # -----------------------------
 # 4. Project to UTM for buffering
 # -----------------------------
 utm_crs <- 32637  # Ethiopia UTM zone
 all_points_proj <- st_transform(all_points, crs = utm_crs)
 
 # Buffer around points (1 km radius)
 buffer_radius <- 1000
 all_buffers <- st_buffer(all_points_proj, dist = buffer_radius)
 
 # -----------------------------
 # 5. Align rasters to reference
 # -----------------------------
 ref <- Access
 Pop_2000 <- project(Pop_2000, ref)
 Ag_land  <- project(Ag_land, ref)
 Ag_suit  <- project(Ag_suit, ref)
 Elev     <- project(Elev, ref)
 
 cov_stack <- c(Ag_suit, Pop_2000, Ag_land, Access, Elev)
 names(cov_stack) <- c("mean_suit", "mean_pop", "mean_agland", "mean_access", "mean_elev")
 
 # -----------------------------
 # 6. Extract mean raster values within buffers
 # -----------------------------
 all_buffers_vect <- vect(all_buffers)
 vals_mean <- terra::extract(cov_stack, all_buffers_vect, fun = mean, na.rm = TRUE)
 
 # Combine with buffer attributes (including Earliest_year)
 rand_df_buffer <- cbind(all_buffers, vals_mean[, -1])
 rand_df_buffer_clean <- rand_df_buffer %>%
   filter(complete.cases(st_drop_geometry(.)))
 
 
 #############################################################################################
 
 
 
 # -----------------------------
 # Ecoregion analysis: remoteness vs PA coverage
 # -----------------------------

 
 library(dplyr)
 library(ggplot2)
 library(terra)
 library(viridis)
 library(tidyr)
 library(sf)
 
 # -----------------------------
 # Load ecoregions shapefile
 # -----------------------------
 ecoregions <- st_read("Eth_2017_ecoregions_141024.shp")
 ecoregions <- st_make_valid(ecoregions)
 eco_vect <- vect(ecoregions)  # convert to SpatVector
 
 # -----------------------------
 # Load ecoregion PA coverage with colour info
 # -----------------------------
 eco_coverage <- read.csv("Ecoregion_PA_coverage.csv")  # columns: ECO_NAME, proportion_covered, colour
 
 # -----------------------------
 # List of covariates (already in cov_stack)
 # -----------------------------
 covariates <- c("mean_access", "mean_agland", "mean_suit", "mean_elev", "mean_pop")
 
 # -----------------------------
 # Extract raster values per ecoregion
 # -----------------------------
 vals_eco <- terra::extract(cov_stack, eco_vect)
 vals_eco <- cbind(ECO_NAME = rep(ecoregions$ECO_NAME, times = table(vals_eco$ID)), vals_eco[, -1])
 
 # Convert to long format for ggplot
 vals_long <- vals_eco %>%
   pivot_longer(cols = all_of(covariates), names_to = "Covariate", values_to = "Value")
 
 # Compute mean ± SE per ecoregion and covariate
 eco_stats <- vals_long %>%
   group_by(ECO_NAME, Covariate) %>%
   summarise(
     mean_val = mean(Value, na.rm = TRUE),
     se_val   = sd(Value, na.rm = TRUE)/sqrt(n()),
     .groups = "drop"
   )
 
 # Add proportion_covered as another "covariate" with SE = 0
 eco_coverage_long <- eco_coverage %>%
   rename(mean_val = proportion_covered) %>%
   mutate(Covariate = "proportion_covered", se_val = 0)
 
 # Combine stats
 eco_stats_all <- bind_rows(eco_stats, eco_coverage_long)
 
 # Merge colour info
 eco_stats_all <- eco_stats_all %>%
   left_join(eco_coverage %>% select(ECO_NAME, colour), by = "ECO_NAME")
 
 # Order ECO_NAME by proportion_covered
 eco_order <- eco_coverage %>%
   arrange(proportion_covered) %>%
   pull(ECO_NAME)
 eco_stats_all$ECO_NAME <- factor(eco_stats_all$ECO_NAME, levels = eco_order)
 
 # -----------------------------
 # Plot multi-panel figure
 # -----------------------------
 ggplot(eco_stats_all, aes(x = ECO_NAME, y = mean_val, fill = colour.y)) +
   geom_col() +
   geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), width = 0.2) +
   coord_flip() +
   facet_wrap(~Covariate, scales = "free_x") +
   scale_fill_identity() +  # use colours from CSV
   labs(
     x = "Ecoregion",
     y = "Mean ± SE",
     title = "Ecoregion covariates and PA coverage across Ethiopia"
   ) +
   theme_minimal() +
   theme(legend.position = "none")
 
 library(sf)
 library(terra)
 library(dplyr)
 library(ggplot2)
 
 # -----------------------------
 # Threshold for underrepresented
 # -----------------------------
 threshold <- 9.4  # e.g., 9.4% coverage
 
 # Merge coverage info with ecoregions
 ecoregions <- st_read("Eth_2017_ecoregions_141024.shp") %>% st_make_valid()
 eco_coverage <- read.csv("Ecoregion_PA_coverage.csv")
 ecoregions <- ecoregions %>%
   left_join(eco_coverage, by = "ECO_NAME") %>%
   mutate(coverage_group = ifelse(proportion_covered < threshold, 
                                  "Underrepresented", "Well-represented"))
 
 # -----------------------------
 # Sample random points
 # -----------------------------
 set.seed(42)
 n_points_per_group <- 5000
 
 points_list <- list()
 for (grp in c("Underrepresented", "Well-represented")) {
   grp_polys <- ecoregions %>% filter(coverage_group == grp)
   points <- st_sample(grp_polys, size = n_points_per_group, type = "random")
   points_sf <- st_sf(geometry = points)
   points_sf$coverage_group <- grp
   points_list[[grp]] <- points_sf
 }
 
 all_points <- do.call(rbind, points_list)
 
 ggplot() +
   geom_sf(data = ecoregions, fill = "grey95", color = "grey70", size = 0.2) +
   geom_sf(data = all_points, aes(color = coverage_group), size = 0.1, alpha = 0.6) +
   scale_color_manual(values = c(
     "Underrepresented" = "#E63946",  # red
     "Well-represented" = "#1D3557"   # blue
   )) +
   theme_minimal() +
   labs(color = "Coverage group",
        title = "Random sample of underrepresented vs well-represented ecoregions")
 
 # -----------------------------
 # Project to raster CRS
 # -----------------------------
 all_points_vect <- vect(st_transform(all_points, crs(cov_stack)))
 plot(all_points_vect)
 
 # -----------------------------
 # Extract covariate values
 # -----------------------------
 vals_points <- terra::extract(cov_stack, all_points_vect, fun = mean, na.rm = TRUE)
 rand_df <- cbind(all_points, vals_points[, -1])
 rand_df <- rand_df[complete.cases(rand_df), ]
 
 # -----------------------------
 # Compare underrepresented vs well-represented
 # -----------------------------
 covariates <- c("mean_suit", "mean_pop", "mean_agland", "mean_access", "mean_elev")
 for (v in covariates) {
   cat("\n---", v, "---\n")
   t_res <- t.test(rand_df[[v]] ~ rand_df$coverage_group)
   print(t_res)
   
   # Optional: boxplot
   p <- ggplot(rand_df, aes_string(x = "coverage_group", y = v, fill = "coverage_group")) +
     geom_boxplot(alpha = 0.7) +
     labs(x = "Ecoregion coverage group", y = v, title = paste(v, "by coverage group")) +
     scale_fill_viridis_d(option = "plasma") +
     theme_minimal()
   print(p)
 }

 
 library(ggplot2)
 library(dplyr)
 library(tidyr)
 library(gridExtra)
 
 # -----------------------------
 # Prepare data for plotting
 # -----------------------------
 covariates <- c("mean_suit", "mean_pop", "mean_agland", "mean_access", "mean_elev")
 
 # Convert to long format for multi-panel plotting
 plot_df <- rand_df %>%
   select(all_of(covariates), coverage_group) %>%
   pivot_longer(cols = all_of(covariates),
                names_to = "Covariate",
                values_to = "Value")
 
 # Compute mean and SE for each covariate and coverage group
 summary_df <- plot_df %>%
   group_by(Covariate, coverage_group) %>%
   summarise(
     mean_val = mean(Value, na.rm = TRUE),
     se_val   = sd(Value, na.rm = TRUE)/sqrt(n()),
     .groups = "drop"
   )
 
 # -----------------------------
 # Plot: multi-panel bar plots with SE
 # -----------------------------
 p <- ggplot(summary_df, aes(x = coverage_group, y = mean_val, fill = coverage_group)) +
   geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
   geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                 width = 0.2, position = position_dodge(0.9)) +
   facet_wrap(~Covariate, scales = "free_y") +
   scale_fill_viridis_d(option = "plasma", end = 0.8) +
   labs(x = "Coverage group", y = "Mean value", title = "Covariates by PA coverage group") +
   theme_minimal() +
   theme(legend.position = "none", strip.text = element_text(face = "bold"))
 
 print(p)
 
 # -----------------------------
 # GLM: Predict underrepresentation
 # -----------------------------
 # Binary outcome: Underrepresented = 1, Well-represented = 0
 rand_df$coverage_bin <- ifelse(rand_df$coverage_group == "Underrepresented", 1, 0)
 
 glm_model <- glm(coverage_bin ~ mean_suit + mean_pop + mean_agland + mean_access + mean_elev,
                  data = rand_df,
                  family = binomial)
 
 summary(glm_model) 
 

 # --- Threshold ---
 threshold <- 9.4  # %
 
 # --- Read and prepare ecoregions ---
 ecoregions <- st_read("Eth_2017_ecoregions_141024.shp") %>%
   st_make_valid()
 
 eco_coverage <- read.csv("Ecoregion_PA_coverage.csv")
 
 ecoregions <- ecoregions %>%
   left_join(eco_coverage, by = "ECO_NAME") %>%
   mutate(
     coverage_group = ifelse(proportion_covered < threshold, 
                             "Underrepresented", "Well-represented"),
     # Distance from threshold (positive value)
     distance_from_threshold = abs(proportion_covered - threshold),
     # Separate weights for under- and over-represented
     underrep_weight = ifelse(proportion_covered < threshold, threshold - proportion_covered, 0),
     overrep_weight  = ifelse(proportion_covered > threshold, proportion_covered - threshold, 0)
   )
 
 # --- Sampling setup ---
 set.seed(42)
 n_points_total <- 500  # total points (split by group)
 ecoregions <- ecoregions %>%
   mutate(area_km2 = as.numeric(st_area(.)) / 1e6)
 
 # helper to allocate counts that sum exactly
 allocate_counts <- function(w, total) {
   w <- ifelse(is.na(w) | w < 0, 0, w)
   if (sum(w) == 0) rep(0, length(w)) else as.vector(rmultinom(1, size = total, prob = w / sum(w)))
 }
 
 # --- UNDERREPRESENTED ---
 eco_under <- ecoregions %>%
   filter(coverage_group == "Underrepresented") %>%
   mutate(weight_area = underrep_weight * area_km2)
 
 eco_under$n_points <- allocate_counts(eco_under$weight_area, n_points_total / 2)
 
 under_points_list <- vector("list", nrow(eco_under))
 for (i in seq_len(nrow(eco_under))) {
   n <- eco_under$n_points[i]
   if (n > 0) {
     pts <- st_sample(eco_under[i, ], size = n, type = "random")
     if (length(pts) > 0) {
       under_points_list[[i]] <- st_sf(
         geometry = pts,
         coverage_group = "Underrepresented",
         weight_value = eco_under$underrep_weight[i],
         area_km2 = eco_under$area_km2[i]
       )
     }
   }
 }
 under_points <- do.call(rbind, under_points_list)
 
 # --- WELL-REPRESENTED ---
 # --- WELL-REPRESENTED ---
 eco_well <- ecoregions %>%
   filter(coverage_group == "Well-represented") %>%
   mutate(weight_area = overrep_weight * area_km2)
 
 eco_well$n_points <- allocate_counts(eco_well$weight_area, n_points_total / 2)
 
 well_points_list <- vector("list", nrow(eco_well))
 for (i in seq_len(nrow(eco_well))) {
   n <- eco_well$n_points[i]
   if (n > 0) {
     pts <- st_sample(eco_well[i, ], size = n, type = "random")
     if (length(pts) > 0) {
       well_points_list[[i]] <- st_sf(
         geometry = pts,
         coverage_group = "Well-represented",
         weight_value = eco_well$overrep_weight[i],
         area_km2 = eco_well$area_km2[i]
       )
     }
   }
 }
 well_points <- do.call(rbind, well_points_list)
 
 # --- Combine all points ---
 all_points <- rbind(under_points, well_points)
 
 # --- Plot ---
 ggplot() +
   geom_sf(data = ecoregions, fill = "grey95", color = "grey70", size = 0.2) +
   geom_sf(data = all_points, aes(color = coverage_group, alpha = weight_value), size = 0.6) +
   scale_color_manual(values = c(
     "Underrepresented" = "#E63946",  # red
     "Well-represented" = "#1D3557"   # blue
   )) +
   scale_alpha(range = c(0.3, 1)) +
   theme_minimal() +
   labs(
     color = "Coverage group",
     alpha = "Degree from threshold",
     title = "Weighted sampling by distance from 9.4% coverage threshold"
   )
 
 # ======================================================
 # 🧩 CONTINUE FROM YOUR LAST LINE
 # ======================================================
 
 # --- Project to raster CRS ---
 all_points_vect <- terra::vect(st_transform(all_points, crs(cov_stack)))
 
 # --- Extract covariate values ---
 vals_points <- terra::extract(cov_stack, all_points_vect, fun = mean, na.rm = TRUE)
 
 # Combine covariates with sampled points
 rand_df <- cbind(all_points, vals_points[, -1])
 rand_df <- rand_df[complete.cases(rand_df), ]
 
 # ======================================================
 # 📊 Compare underrepresented vs well-represented
 # ======================================================
 
 # Define the covariates you extracted
 covariates <- names(cov_stack)  # or manually: c("mean_suit", "mean_pop", ...)
 
 for (v in covariates) {
   cat("\n---", v, "---\n")
   t_res <- t.test(rand_df[[v]] ~ rand_df$coverage_group)
   print(t_res)
   
   # Boxplot
   p <- ggplot(rand_df, aes_string(x = "coverage_group", y = v, fill = "coverage_group")) +
     geom_boxplot(alpha = 0.7, outlier.shape = NA) +
     labs(x = "Coverage group", y = v, title = paste(v, "by coverage group")) +
     scale_fill_manual(values = c("Underrepresented" = "#E63946", "Well-represented" = "#1D3557")) +
     theme_minimal()
   print(p)
 }
 
 # ======================================================
 # 📈 Multi-panel bar plot of mean ± SE
 # ======================================================
 plot_df <- rand_df %>%
   select(all_of(covariates), coverage_group) %>%
   pivot_longer(cols = all_of(covariates), names_to = "Covariate", values_to = "Value")
 
 summary_df <- plot_df %>%
   group_by(Covariate, coverage_group) %>%
   summarise(
     mean_val = mean(Value, na.rm = TRUE),
     se_val   = sd(Value, na.rm = TRUE) / sqrt(n()),
     .groups = "drop"
   )
 
 p <- ggplot(summary_df, aes(x = coverage_group, y = mean_val, fill = coverage_group)) +
   geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
   geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                 width = 0.2, position = position_dodge(0.9)) +
   facet_wrap(~Covariate, scales = "free_y") +
   scale_fill_manual(values = c("Underrepresented" = "#E63946", "Well-represented" = "#1D3557")) +
   labs(x = "Coverage group", y = "Mean value", title = "Covariates by coverage group (weighted samples)") +
   theme_minimal() +
   theme(legend.position = "none", strip.text = element_text(face = "bold"))
 
 print(p)
 
 # ======================================================
 # 🔍 Logistic regression: What predicts underrepresentation?
 # ======================================================
 rand_df$coverage_bin <- ifelse(rand_df$coverage_group == "Underrepresented", 1, 0)
 
 # Optional: standardize covariates for comparability
 num_covs <- covariates[sapply(rand_df[covariates], is.numeric)]
 rand_df[num_covs] <- scale(rand_df[num_covs])
 
 glm_formula <- as.formula(
   paste("coverage_bin ~", paste(covariates, collapse = " + "))
 )
 
 glm_model <- glm(glm_formula, data = rand_df, family = binomial)
 summary(glm_model)
 
 # ======================================================
 # 🗺️ Optional: Map covariate patterns (if needed)
 # ======================================================
 # Example — visualize one covariate’s mean per ecoregion:
 eco_cov_summary <- rand_df %>%
   st_drop_geometry() %>%
   group_by(coverage_group) %>%
   summarise(across(all_of(covariates), mean, na.rm = TRUE))
 
 eco_cov_summary
 
 
 # BRING IN PERFORMANCE BUDGET DAT
 
 env_dat <- read.csv("env_performance_budget_data.csv")
 env_dat <- env_dat %>% rename(NAME_PAs = Name_PAs)
 pov_dat <- read.csv("pov_performance_budget_data.csv")
 pov_dat <- pov_dat %>% rename(NAME_PAs = Name_PAs)
 
 # =====================================================
 # 🧩 Combine PA covariates, env & pov data, and main ecoregion
 # =====================================================
 
 library(sf)
 library(dplyr)
 
 crs(ecoregions)
 ext(ecoregions)
 ecoregions_vect <- vect(ecoregions)
 
 # Reproject PAs to match ecoregions
 PA_vect_utm <- project(PA_vect, crs(ecoregions))
 
 crs(ecoregions_vect)
 crs(PA_vect_utm)
 ext(ecoregions_vect)
 ext(PA_vect_utm)
 
 
 invalid_idx <- !is.valid(PA_vect_utm)
 PA_vect_utm[invalid_idx, ]  # inspect these
 # 2. Attempt to fix invalid geometries
PA_vect_fixed <- makeValid(PA_vect_utm)
 # 3. Check if all are now valid
 all(is.valid(PA_vect_fixed))  # should return TRUE
 
 # Then intersect
 PA_main_eco_vect <- terra::intersect(PA_vect_fixed, ecoregions_vect)

 PA_main_eco_vect$intersect_area <- expanse(PA_main_eco_vect)  # area in map units
 
 # Get majority ecoregion per PA
 PA_main_eco_df <- as.data.frame(PA_main_eco_vect) %>%
   group_by(NAME_PAs) %>%
   slice_max(order_by = intersect_area, n = 1, with_ties = FALSE) %>%
   ungroup() %>%
   select(NAME_PAs, Main_ecoregion = ECO_NAME)
 
 
 # --- Optional: check results ---
 head(PA_main_eco_df)
 
 # -----------------------------
 # Load rasters
 # -----------------------------
 Ag_suit   <- rast("agri_suit_rain_clipped_1980-2009.tif")
 Pop_2000  <- rast("eth_pd_2000_1km_UNadj.tif")
 Ag_land   <- rast("covariate_baseline_agri_2003.tif")
 Access   <- rast("updated_access_2000.tif")
 Elev     <- rast("Elev.tif")
 Temp <- rast("Temp.tif")
 Ppt <- rast("ppt.tif")


 # Set reference raster
 ref <- Access
 
 # Align all rasters to the reference
 Ag_suit   <- project(Ag_suit, ref)
 Pop_2000  <- project(Pop_2000, ref)
 Ag_land   <- project(Ag_land, ref)
 Elev      <- project(Elev, ref)
 Temp      <- project(Temp, ref)
 Ppt       <- project(Ppt, ref)
 cov_stack2 <- c(Ag_suit, Pop_2000, Ag_land, Access, Elev, Temp, Ppt)
 names(cov_stack2) <- c("mean_suit", "mean_pop", "mean_agland",
                       "mean_access", "mean_elev", "mean_temp", "mean_ppt")
 
 # -----------------------------
 # 2. Extract raster values per PA
 # -----------------------------
 vals_all <- terra::extract(cov_stack2, PA_vect, fun = mean, na.rm = TRUE)
 
 # Add PA names
 vals_all_df <- vals_all %>%
   mutate(NAME_PAs = PA_shp$NAME_PAs[ID]) %>%
   group_by(NAME_PAs) %>%
   summarise(across(mean_suit:mean_ppt, mean, na.rm = TRUE)) %>%
   ungroup()
 
 # -----------------------------
 # 3. Combine with PA-year data
 # -----------------------------
 PA_summary2 <- PA_yr %>%
   left_join(vals_all_df, by = "NAME_PAs")
 
 
 
 # --- Combine with environmental performance data ---
 env_full <- env_dat %>%
   left_join(PA_main_eco_df, by = "NAME_PAs") %>%
   left_join(PA_summary2,  by = "NAME_PAs")
 
 # --- Combine with poverty / budget data ---
 pov_full <- pov_dat %>%
   left_join(PA_main_eco_df, by = "NAME_PAs") %>%
   left_join(PA_summary2,  by = "NAME_PAs")
 
 # --- (optional) Inspect results ---
 glimpse(env_full)
 glimpse(pov_full)

 library(vegan)
 library(dplyr)
 
 # -----------------------------
 # 1. Numeric PA covariates
 # -----------------------------
 PA_covs <- c("mean_suit", "mean_pop", "mean_agland",
              "mean_access", "mean_elev", "mean_temp", "mean_ppt", "AREA_km2") 

 
 # Ensure ecoregion is a factor
 env_full$Main_ecoregion <- as.factor(env_full$Main_ecoregion)
 pov_full$Main_ecoregion <- as.factor(pov_full$Main_ecoregion)
 
 
 # Select only numeric predictors
 numeric_vars <- env_full[, c("budget_resid", "mean_suit", "mean_pop", 
                              "mean_agland", "mean_access", "mean_elev", 
                              "mean_temp", "mean_ppt", "AREA_km2")]
 
 # Compute correlation matrix
 cor_matrix <- cor(numeric_vars, use = "complete.obs")
 round(cor_matrix, 2)
 
 # Optional: visualize correlations
 library(corrplot)
 corrplot(cor_matrix, method = "circle", type = "upper")
 
 plot(env_full$budget_resid, env_full$biodiversity)
 
 lm_model <- lm(biodiversity ~ budget_resid + mean_temp + mean_elev + mean_ppt, 
                data = env_full)
 
 
 plot(env_full$budget_resid, env_full$biodiversity,
      xlab = "Budget Residual", ylab = "Biodiversity",
      main = "Scatter plot for linearity check")
 abline(lm(biodiversity ~ budget_resid, data = env_full), col = "red")
 
 shapiro.test(env_full$biodiversity)
 shapiro.test(env_full$budget_resid)
 
 cor.test(env_full$biodiversity, env_full$budget_resid, method = "pearson")
 cor.test(env_full$biodiversity, env_full$budget_resid, method = "spearman")
 library(ggplot2)
 
 ggplot(env_full, aes(x = budget_resid, y = biodiversity)) +
   geom_point() +
   geom_smooth(method = "loess", se = TRUE, color = "blue") +
   labs(x = "Budget Residual", y = "Biodiversity",
        title = "Scatterplot with LOESS smoothing")

 # Basic boxplot of biodiversity by ecoregion
 ggplot(pov_full, aes(x = Main_ecoregion, y = poverty, fill = Main_ecoregion)) +
   geom_boxplot(alpha = 0.7) +
   geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +  # Points
   geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +  # Reference line at 0
   theme_classic() +
   labs(
     x = "Ecoregion",
     y = "Wellbeing performance"
   ) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1),
         legend.position = "none")
 
 ggplot(env_full, aes(x = Main_ecoregion, y = biodiversity, fill = Main_ecoregion)) +
   geom_boxplot(alpha = 0.7) +
   geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.6) +  # Points
   geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1) +  # Reference line at 0
   theme_classic() +
   labs(
     x = "Ecoregion",
     y = "Biodiversity performance"
   ) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1),
         legend.position = "none")
 

 # Use dplyr explicitly to avoid conflicts
 env_full_plot <- dplyr::select(env_full, Main_ecoregion, biodiversity) %>%
   dplyr::rename(Performance = biodiversity) %>%
   dplyr::mutate(Metric = "Biodiversity")
 
 pov_full_plot <- dplyr::select(pov_full, Main_ecoregion, poverty) %>%
   dplyr::rename(Performance = poverty) %>%
   dplyr::mutate(Metric = "Wellbeing")
 
 combined_plot <- bind_rows(env_full_plot, pov_full_plot)
 
 # Plot
 ggplot(combined_plot, aes(x = Main_ecoregion, y = Performance, fill = Metric)) +
   geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
   geom_jitter(aes(color = Metric),
               size = 2, alpha = 0.6,
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
   theme_classic() +
   labs(x = "Ecoregion", y = "Performance", title = "Biodiversity vs Wellbeing Across Ecoregions") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 
 summary(lm_model)
 anova(lm_model)
 library(car)
 vif(lm_model)
 
 # Example for elevation
 boxplot(biodiversity ~ Main_ecoregion, data = env_full)
 
 # Example for temperature
 boxplot(mean_temp ~ Main_ecoregion, data = env_full)
 library(lsr)
 etaSquared(aov(mean_elev ~ Main_ecoregion, data = env_full))
 etaSquared(aov(mean_temp ~ Main_ecoregion, data = env_full))
 
 library(MuMIn)
 
 # Make sure options are set for LM
 options(na.action = "na.fail")  
 
 # Dredge the model
 dredge_results <- dredge(lm_model)
 
 # View top models
 head(dredge_results)
 
 # Optionally, get model averaging
 avg_model <- model.avg(dredge_results, subset = delta < 2)
 summary(avg_model) 
 
 best_model <- lm(biodiversity ~ budget_resid,
                  data = env_full)
 
 summary(best_model)
 anova(best_model)
 
 
 top_models <- dredge_results[1:3, ]  # top 3 models
 model_avg <- model.avg(get.models(dredge_results, subset = delta < 2), 
                        revised.var = TRUE, fit = TRUE)
 summary(model_avg)
 par(mfrow = c(1,1))
 plot(best_model)
 car::vif(best_model)
 
 
 
 # -----------------------------
 # 2. Environmental outcomes RDA
 # -----------------------------
 Y_env <- env_full %>% select(Forest_ATE, Grass_ATE, Agri_ATE)
 
 X_env <- env_full %>%
   select(budget_resid, all_of(PA_covs), Main_ecoregion)
 
 lm_model <- lm(biodiversity ~ budget_resid +  mean_ppt +  Main_ecoregion, 
                data = env_full)
 summary(lm_model)
 anova(lm_model)
 
 # Full RDA: budget_resid + PA covariates + ecoregion
 rda.env_full <- rda(env_full$biodiversity ~ ., data = X_env, scale = FALSE)
 set.seed(200)
 anova.cca(rda.env_full, step = 1000, by = "term")
 
 # -----------------------------
 # 3. Poverty outcomes RDA
 # -----------------------------
 Y_pov <- pov_full %>% select(MAHFP_ATE, HDDS_ATE, asset_ATE)
 
 X_pov <- pov_full %>%
   select(budget_resid, all_of(PA_covs), Main_ecoregion)
 
 # Full RDA: budget_resid + PA covariates + ecoregion
 rda.pov_full <- rda(Y_pov ~ ., data = X_pov, scale = FALSE)
 set.seed(200)
 anova.cca(rda.pov_full, step = 1000, by = "term")
 
 # -----------------------------
 # 4. Adjusted R² for each model
 # -----------------------------
 RsquareAdj(rda.env_full)
 RsquareAdj(rda.pov_full)
 
 # Environmental
 rda.env_null <- rda(Y_env ~ 1, data = X_env)
 rda.env_best <- ordiR2step(rda.env_null, scope = formula(rda.env_full),
                            direction = "both", R2scope = TRUE, perm.max = 999)#
 # Check permutation significance for each term
 anova.cca(rda.env_best, by = "term", permutations = 999)

 
 # Poverty
 rda.pov_null <- rda(Y_pov ~ 1, data = X_pov)
 rda.pov_best <- ordiR2step(rda.pov_null, scope = formula(rda.pov_full),
                            direction = "both", R2scope = TRUE, perm.max = 999)
 # Check permutation significance for each term
 anova.cca(rda.env_best, by = "term", permutations = 999)
 
 rda.env_sig <- rda(Y_env ~ mean_access + Main_ecoregion + mean_ppt, data = X_env, scale = FALSE)
 anova.cca(rda.env_sig, permutations = 999, by = "term")
 
 ##### budget residuals explained by area and accessibility
 
 # Environmental
 lm_budget_env <- lm(Avr_budg_scale_2014_infla_USD ~ AREA_km2 + mean_access + Main_ecoregion, data = env_full)
 env_full$budget_resid2 <- residuals(lm_budget_env)
 
 # Poverty
 lm_budget_pov <- lm(Avr_budg_scale_2014_infla_USD ~ AREA_km2 + mean_access + Main_ecoregion, data = pov_full)
 pov_full$budget_resid2 <- residuals(lm_budget_pov)

 rda.env <- rda(Y_env~ budget_resid2,
                data = env_full,
                 scale= F) 
 set.seed(40)
 anova.cca(rda.env, step = 1000, by = "term")
 RsquareAdj(rda.env)
 
 rda.env_best <- rda(Y_env ~ budget_resid + mean_access + mean_ppt + Main_ecoregion,
                     data = X_env, scale = FALSE)
 anova.cca(rda.env_best, permutations = 999, by = "term")
 RsquareAdj(rda.env_best)
 
 library(car)
 vif(lm(budget_resid ~ mean_access + mean_ppt + Main_ecoregion, data = X_env))
 
 rda.env1 <- rda(Y_env ~ budget_resid, data = X_env, scale = F)
 anova.cca(rda.env1, step = 1000, by = "term")
 
 rda.env2 <- rda(Y_env ~ budget_resid + Main_ecoregion, data = X_env)
 anova.cca(rda.env2, by = "term", permutations = 999)
 
 rda.env3 <- rda(Y_env ~ budget_resid + mean_access + Main_ecoregion, data = X_env)
 anova.cca(rda.env3, by = "term", permutations = 999)
 
 
 library(ggplot2)
 library(dplyr)
 library(tidyr)
 library(broom)
 
 library(ggplot2)
 library(dplyr)
 
 # Median values for numeric predictors not being plotted
 median_budget <- median(env_full$budget_resid)
 median_ppt    <- median(env_full$mean_ppt)
 median_area   <- median(env_full$AREA_km2)
 
 # -------------------------------
 # A) Ecoregion effect
 # -------------------------------
 eco_df <- data.frame(
   Main_ecoregion = unique(env_full$Main_ecoregion),
   budget_resid = median_budget,
   mean_ppt = median_ppt,
   AREA_km2 = median_area
 )
 eco_pred <- cbind(eco_df, predict(best_model, eco_df, interval = "confidence"))
 
 ggplot(eco_pred, aes(x = Main_ecoregion, y = fit)) +
   geom_point(size = 3, color = "steelblue") +
   geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2, color = "steelblue") +
   theme_minimal() +
   labs(x = "Ecoregion", y = "Predicted biodiversity",
        title = "Effect of Main Ecoregion on Biodiversity") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 # -------------------------------
 # B) Precipitation effect
 # -------------------------------
 ppt_df <- data.frame(
   mean_ppt = seq(min(env_full$mean_ppt), max(env_full$mean_ppt), length.out = 50),
   budget_resid = median_budget,
   AREA_km2 = median_area,
   Main_ecoregion = factor("Ethiopian montane grasslands and woodlands", 
                           levels = levels(env_full$Main_ecoregion))
 )
 ppt_pred <- cbind(ppt_df, predict(best_model, ppt_df, interval = "confidence"))
 
 ggplot(ppt_pred, aes(x = mean_ppt, y = fit)) +
   geom_line(color = "darkgreen", size = 1.2) +
   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "darkgreen") +
   theme_minimal() +
   labs(x = "Mean precipitation", y = "Predicted biodiversity",
        title = "Effect of Precipitation on Biodiversity")
 
 # -------------------------------
 # C) Budget effect
 # -------------------------------
 budget_df <- data.frame(
   budget_resid = seq(min(env_full$budget_resid), max(env_full$budget_resid), length.out = 50),
   mean_ppt = median_ppt,
   AREA_km2 = median_area,
   Main_ecoregion = factor("Ethiopian montane grasslands and woodlands", 
                           levels = levels(env_full$Main_ecoregion))
 )
 budget_pred <- cbind(budget_df, predict(best_model, budget_df, interval = "confidence"))
 
 ggplot(budget_pred, aes(x = budget_resid, y = fit)) +
   geom_line(color = "purple", size = 1.2) +
   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "purple") +
   theme_minimal() +
   labs(x = "Residual budget", y = "Predicted biodiversity",
        title = "Effect of Budget on Biodiversity")
 
 # -------------------------------
 # D) Area effect
 # -------------------------------
 area_df <- data.frame(
   AREA_km2 = seq(min(env_full$AREA_km2), max(env_full$AREA_km2), length.out = 50),
   budget_resid = median_budget,
   mean_ppt = median_ppt,
   Main_ecoregion = factor("Ethiopian montane grasslands and woodlands", 
                           levels = levels(env_full$Main_ecoregion))
 )
 area_pred <- cbind(area_df, predict(best_model, area_df, interval = "confidence"))
 
 ggplot(area_pred, aes(x = AREA_km2, y = fit)) +
   geom_line(color = "orange", size = 1.2) +
   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, fill = "orange") +
   theme_minimal() +
   labs(x = "Area (km²)", y = "Predicted biodiversity",
        title = "Effect of Area on Biodiversity")
 
 library(nlme)
 
 # Global model with ecoregion as a random effect
 global_model_re <- lme(
   biodiversity ~ 
     budget_resid + 
   +mean_suit + mean_pop + mean_agland + 
     mean_access + mean_elev +  mean_ppt + AREA_km2,
     random = ~ 1 | Main_ecoregion, 
   method = "ML",
   data = env_full
 )
 
 plot(global_model_re)
 plot(global_model_re)             # residuals vs fitted
 qqnorm(resid(global_model_re)); qqline(resid(global_model_re))  # normality
 
 # Inspect summary
 summary(global_model_re)
 VarCorr(global_model_re)  
 
 library(MuMIn)
 options(na.action = "na.fail")  # required
 
 dredge_results_nlme <- dredge(global_model_re)
 avg_model <- model.avg(dredge_results_nlme, subset = delta < 2)
 summary(avg_model)

 # Get the best model (lowest AICc)
 best_model_ml <- get.models(dredge_results_nlme, subset = 1)[[1]]
 
 # Refit best model with REML for final estimates
 best_model_reml <- update(best_model_ml, method = "REML")
 summary(best_model_reml)
 
 plot(best_model_reml)             # residuals vs fitted
 qqnorm(resid(best_model_reml)); qqline(resid(best_model_reml))  # normality
 
 
 library(nlme)
 library(ggplot2)
 library(dplyr)
 
 # Extract random effects (intercepts) for each ecoregion
 re_effects <- ranef(global_model_re)
 re_df <- data.frame(
   Main_ecoregion = rownames(re_effects),
   RandomIntercept = re_effects$`(Intercept)`
 )
 
 # Merge with overall fixed-effect intercept to get adjusted means
 fixed_intercept <- fixef(global_model_re)["(Intercept)"]
 re_df <- re_df %>%
   mutate(AdjustedMean = fixed_intercept + RandomIntercept)
 
 # Order by adjusted mean for plotting
 re_df <- re_df %>% arrange(AdjustedMean)
 re_df$Main_ecoregion <- factor(re_df$Main_ecoregion, levels = re_df$Main_ecoregion)
 
 # Plot
 ggplot(re_df, aes(x = Main_ecoregion, y = AdjustedMean)) +
   geom_col(fill = "steelblue") +
   geom_text(aes(label = round(AdjustedMean,2)), vjust = -0.5) +
   theme_minimal() +
   labs(
     x = "Ecoregion",
     y = "Predicted biodiversity (adjusted for fixed effects)",
     title = "Random intercepts: biodiversity variation across ecoregions"
   ) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
 
 # One-way ANOVA: biodiversity ~ ecoregion
 anova_ecoregion <- aov(biodiversity ~ Main_ecoregion, data = env_full)
 
 # Summary of ANOVA
 summary(anova_ecoregion)
 # Tukey HSD test
 tukey_results <- TukeyHSD(anova_ecoregion)
 
 # View results
 tukey_results
 
 table(env_full$Main_ecoregion)
 