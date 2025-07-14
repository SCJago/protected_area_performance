# =============================================================================
# Protected Area Network Expansion Timeline in Ethiopia
# =============================================================================

# --- Load Required Libraries ---
library(ggplot2)
library(scales)
library(grid)
library(ggprism)
library(sf)
library(dplyr)
library(raster)

# =============================================================================
# 1. Import Data
# =============================================================================

# Load PA shapefile
ETH_PAs <- st_read("Eth_PAs_sep24.shp") %>%
  st_zm(drop = TRUE, what = "ZM") %>%
  st_transform(crs = 32637) %>%  # UTM Zone 37N
  st_make_valid()

# Calculate area in km²
ETH_PAs$area_km2 <- as.numeric(st_area(ETH_PAs) / 1e6)

# Join PA year info
PA_year <- read.csv("Earliest_PA_year.csv")
ETH_PAs <- ETH_PAs %>%
  left_join(PA_year, by = c("NAME_PAs" = "Name_PAs"))

# Replace missing IUCN categories with VI
ETH_PAs <- ETH_PAs %>%
  mutate(IUCN_MG_CA = if_else(is.na(IUCN_MG_CA), "VI", IUCN_MG_CA))

# Create supplement table
supplement_table <- ETH_PAs %>%
  select(NAME_PAs, DESIGNATIO, IUCN_MG_CA, area_km2, Earliest_year) %>%
  arrange(NAME_PAs) %>%
  st_drop_geometry()

write.csv(supplement_table, "Eth_PAs_supplement_table.csv", row.names = FALSE)

# =============================================================================
# 2. Calculate Cumulative Area and Number of PAs Over Time
# =============================================================================

cum_area_df <- data.frame(
  Year = sort(unique(ETH_PAs$Earliest_year)),
  Cum_Area = 0,
  Cum_Num = 0
)

cumulative_polygons <- st_sfc(crs = st_crs(ETH_PAs))
cum_num <- 0

for (i in seq_along(cum_area_df$Year)) {
  current_year_polygons <- ETH_PAs %>% filter(Earliest_year == cum_area_df$Year[i])
  cum_num <- cum_num + nrow(current_year_polygons)
  cum_area_df$Cum_Num[i] <- cum_num
  
  unioned <- st_union(current_year_polygons)
  new_area <- st_difference(unioned, cumulative_polygons)
  area_km2 <- as.numeric(st_area(new_area)) / 1e6
  
  cumulative_polygons <- st_union(cumulative_polygons, unioned)
  cum_area_df$Cum_Area[i] <- ifelse(i == 1, area_km2, cum_area_df$Cum_Area[i - 1] + area_km2)
}

# Add % national coverage
Eth_area <- 1127095  # Total Ethiopia area in km²
cum_area_df$Percent_area <- (cum_area_df$Cum_Area / Eth_area) * 100

# Add row for 2024 (if unchanged since 2018)
cum_area_df <- rbind(cum_area_df, data.frame(
  Year = 2024,
  Cum_Area = cum_area_df$Cum_Area[cum_area_df$Year == 2018],
  Cum_Num = cum_area_df$Cum_Num[cum_area_df$Year == 2018],
  Percent_area = cum_area_df$Percent_area[cum_area_df$Year == 2018]
))

# =============================================================================
# 3. Plot: PA Network Expansion Over Time
# =============================================================================

numberColor <- "black"
areaColor <- "#999999"

PA_expansion_plot <- ggplot(cum_area_df, aes(x = Year)) +
  # Period shading
  geom_rect(aes(xmin = 1958, xmax = 1980, ymin = -10, ymax = -1), fill = "#CC79A7") +
  geom_rect(aes(xmin = 1980, xmax = 2000, ymin = -10, ymax = -1), fill = "#E69F00") +
  geom_rect(aes(xmin = 2000, xmax = 2024, ymin = -10, ymax = -1), fill = "#009E73") +

  # Data: Area & PA count
  geom_area(aes(y = Percent_area * 10), fill = areaColor) +  # rescaled for dual axis
  geom_line(aes(y = Cum_Num), size = 0.8, color = numberColor) +
  geom_point(aes(y = Cum_Num), size = 1, color = numberColor) +

  # Policy annotations
  geom_segment(aes(x = 1962, y = -22, xend = 1962, yend = -29), arrow = arrow(length = unit(0.2, 'cm'))) +
  geom_segment(aes(x = 1994, y = -13, xend = 1994, yend = -29), arrow = arrow(length = unit(0.2, 'cm'))) +
  geom_segment(aes(x = 2000, y = -22, xend = 2000, yend = -29), arrow = arrow(length = unit(0.2, 'cm'))) +
  geom_segment(aes(x = 2010, y = -13, xend = 2010, yend = -29), arrow = arrow(length = unit(0.2, 'cm'))) +
  geom_segment(aes(x = 2015, y = -22, xend = 2015, yend = -29), arrow = arrow(length = unit(0.2, 'cm'))) +
  geom_segment(aes(x = 2022, y = -13, xend = 2022, yend = -29), arrow = arrow(length = unit(0.2, 'cm'))) +

  geom_label(label = "Partnered with UNESCO", x = 1966, y = -22, size = 7, fill = "lightgrey") +
  geom_label(label = "Ratified CBD", x = 1994, y = -15, size = 7, fill = "lightgrey") +
  geom_label(label = "Aligned with MDGs", x = 2000, y = -22, size = 7, fill = "lightgrey") +
  geom_label(label = "Adopted Aichi targets", x = 2005, y = -15, size = 7, fill = "lightgrey") +
  geom_label(label = "Aligned with SDGs", x = 2016, y = -22, size = 7, fill = "lightgrey") +
  geom_label(label = "Adopted GBF targets", x = 2017.5, y = -15, size = 7, fill = "lightgrey") +

  # Axes
  scale_y_continuous(
    limits = c(-29, 105),
    breaks = seq(0, 100, 10),
    name = "               Cumulative number of PAs",
    sec.axis = sec_axis(~ . / 10, name = "Cumulative PA coverage (%)               ", breaks = seq(0, 10, 1))
  ) +
  scale_x_continuous(
    limits = c(1958, 2024),
    breaks = seq(1945, 2024, 5),
    guide = "prism_minor",
    minor_breaks = seq(1945, 2024, 1)
  ) +
  theme_classic() +
  theme(
    axis.title.y = element_text(color = numberColor, size = 22),
    axis.title.y.right = element_text(color = areaColor, size = 22),
    axis.text.y.left = element_text(color = numberColor, size = 21),
    axis.text.y.right = element_text(color = areaColor, size = 21),
    axis.text.x = element_text(size = 22),
    axis.title.x = element_text(size = 21),
    axis.line = element_line(linewidth = 1),
    axis.ticks.length = unit(0.2, "cm")
  )

ggsave("PA_expansion_percent_timeline_fig.png", PA_expansion_plot, dpi = 300, h = 27, w = 44, units = "cm")

# =============================================================================
# 4. Strict Protection (IUCN Category II) Cumulative Area
# =============================================================================

catII <- ETH_PAs %>%
  filter(IUCN_MG_CA == "II")

catII_cum_df <- data.frame(Year = sort(unique(catII$Year)), Cum_Area = 0, Cum_Num = 0)
catII_cumulative_polygons <- st_sfc(crs = st_crs(catII))
cum_num <- 0

for (i in seq_along(catII_cum_df$Year)) {
  this_year <- catII %>% filter(Year == catII_cum_df$Year[i])
  cum_num <- cum_num + nrow(this_year)
  catII_cum_df$Cum_Num[i] <- cum_num

  unioned <- st_union(this_year)
  new_area <- st_difference(unioned, catII_cumulative_polygons)
  area_km2 <- as.numeric(st_area(new_area)) / 1e6

  catII_cumulative_polygons <- st_union(catII_cumulative_polygons, unioned)
  catII_cum_df$Cum_Area[i] <- ifelse(i == 1, area_km2, catII_cum_df$Cum_Area[i-1] + area_km2)
}

catII_cum_df$Percent_area <- (catII_cum_df$Cum_Area / Eth_area) * 100
print(catII_cum_df)
