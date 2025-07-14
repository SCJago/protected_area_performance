########################################
# Ethiopia PA Environment Representativeness Analysis
# Authors: Sophie Jago, Bezawit Genanaw, James Borrell
# Description: Assess representativeness of Ethiopia's protected area network across environmental conditions
########################################

# ============ 1. Load Libraries ============
required_packages <- c(
  "raster", "sp", "rgdal", "rgeos", "scales", "ggplot2", "dplyr",
  "factoextra", "vegan", "png", "alphahull", "virtualspecies",
  "ENMTools", "ENMeval", "usdm", "geosphere"
)

lapply(required_packages, require, character.only = TRUE)

# ============ 2. Load and Prepare Data ============

# Load Ethiopia boundary
ETH <- getData("GADM", country = "ETH", level = 0)

# Load climate rasters and crop/mask to Ethiopia
all_raster <- stack("data/chelsa_bioclim_stack.grd")
names(all_raster) <- paste0("bio", 1:19)

all_raster_ETH <- crop(all_raster, ETH)
all_raster_ETH <- mask(all_raster_ETH, ETH)

# Load and reproject protected areas shapefile
pa_shape <- readOGR("data/ETH-WPAs-2023-ed.shp")
pa_shape <- spTransform(pa_shape, CRS(proj4string(all_raster)))

# Crop/mask climate data to PA extent
pa_data <- crop(all_raster, pa_shape)
pa_data <- mask(pa_data, pa_shape)

# ============ 3. Sample Data and Run PCA ============

# Sample background environmental space and PA space
data_for_pca2 <- sampleRandom(all_raster_ETH, 100000, na.rm = TRUE, xy = TRUE)
pa_data_extracted <- sampleRandom(pa_data, 9400, na.rm = TRUE, xy = TRUE)

# PCA on background environmental variables
pca2 <- prcomp(data_for_pca2[, 3:21], center = TRUE, scale. = TRUE)

# Project PAs into PCA space
PAs <- predict(pca2, newdata = pa_data_extracted)

# ============ 4. Plot PCA ============

png("output/Allraster_PAs_PCA_300425.png", units = "cm", width = 15, height = 15, res = 300)
plot(pca2$x[,1], pca2$x[,2], pch = 16, cex = 1, col = alpha("grey", 0.2),
     xlab = "PC1 (40.1%)", ylab = "PC2 (22.5%)")
points(PAs, pch = 16, cex = 1, col = alpha("darkgreen", 0.2))
dev.off()

# ============ 5. Variable Contributions ============

png("output/Contribution_PC1_300425.png", units = "cm", width = 15, height = 15, res = 300)
fviz_contrib(pca2, choice = "var", axes = 1)
dev.off()

png("output/Contribution_PC2_300425.png", units = "cm", width = 15, height = 15, res = 300)
fviz_contrib(pca2, choice = "var", axes = 2)
dev.off()

# ============ 6. Alpha Hull and Coverage ============

unique_PAs <- unique(PAs[, 1:2])
unique_pca2 <- unique(pca2$x[, 1:2])

hull <- ahull(x = unique_PAs[,1], y = unique_PAs[,2], alpha = 0.5)
full_hull <- ahull(unique_pca2, alpha = 0.5)

png("output/alphahulls.png", units = "cm", width = 15, height = 15, res = 300)
plot(pca2$x[,1], pca2$x[,2], pch = 16, cex = 1, col = alpha("grey", 0.2),
     xlab = "PC1 (40.1%)", ylab = "PC2 (22.5%)")
points(PAs, pch = 16, cex = 1, col = alpha("darkgreen", 0.2))
plot(full_hull, add = TRUE, col = "black", wpoints = FALSE)
plot(hull, add = TRUE, col = "darkgreen", wpoints = FALSE)
dev.off()

# ============ 7. Representativeness Area Calculations ============

area_PA_space <- areaahull(hull, timeout = 5)
area_env_space <- areaahull(full_hull, timeout = 5)
percent_overlap <- (area_PA_space / area_env_space) * 100
cat("Percent environmental space covered by PAs:", round(percent_overlap, 2), "%\n")

# ============ 8. Map Missing Environmental Regions ============

# Points outside alpha hull
vector <- inahull(hull, pca2$x[,1:2])
missing_regions <- cbind.data.frame(data_for_pca2[,1:2], vector)

# Plot map of Ethiopia with representativeness
png("output/Representativeness_map.png", units = "cm", width = 15, height = 15, res = 300)
plot(ETH)
missing_regions_ordered <- missing_regions %>% arrange(vector)

colours <- c("grey", "darkgreen")[factor(missing_regions_ordered$vector, levels = c(FALSE, TRUE))]
sizes <- c(2, 0.5)[factor(missing_regions_ordered$vector, levels = c(FALSE, TRUE))]

points(missing_regions_ordered[, 1:2], col = colours, pch = 20, cex = sizes)
dev.off()
