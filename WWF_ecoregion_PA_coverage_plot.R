library(ggplot2) #load ggplot
library(dplyr)
library(sf)
library(raster)
library(geodata)

setwd("~/Papers/Ethiopia - PA history review/Representativeness/WWF_ecoregions") #set WD

# Get the boundaries for Ethiopia
ETH <- gadm(country="ETH", level=0, path=tempdir())
ETH <- st_as_sf(ETH)
plot(ETH)

#get 2017 ecoregions
ecoregion_2017 <- st_read("Ecoregions2017.shp")
head(ecoregion_2017)

# Ensure both datasets are in the same CRS
# Get the CRS of each dataset
eth_crs <- st_crs(ETH)
ecoregion_crs <- st_crs(ecoregion_2017)

# If the CRS are different, transform ecoregion data to match ETH CRS
if (!identical(eth_crs, ecoregion_crs)) {
  ecoregion_2017 <- st_transform(ecoregion_2017, eth_crs)
}

# Check for invalid geometries
invalid_geoms <- st_is_valid(ecoregion_2017)

# If there are invalid geometries, try to make them valid
if (any(!invalid_geoms)) {
  ecoregion_2017 <- st_make_valid(ecoregion_2017)
}

# Clip the ecoregion data to the boundaries of Ethiopia
clipped_ecoregions <- st_intersection(ecoregion_2017, ETH)

# Plot the clipped ecoregions to verify the results
plot(clipped_ecoregions)
head(clipped_ecoregions)

selected_ecoregions <- clipped_ecoregions %>%
  dplyr::select(ECO_NAME, COLOR)

ggplot(data = selected_ecoregions) +
  geom_sf(aes(fill = ECO_NAME), color = NA) +  # No border color
  scale_fill_manual(values = selected_ecoregions$COLOR) +  # Set fill colors from the COLOR column
  labs(title = "Ecoregions of Ethiopia",
       fill = "ECO_NAME") +  # Title and legend label
  theme_minimal() +  # Use a minimal theme
  theme(legend.position = "right")  # Position the legend to the right

############################################################################################

# bring in PA shapefile
ETH_PAs <- st_read("Eth_PAs_sep24.shp")
ETH_PAs <- st_zm(ETH_PAs, drop = TRUE, what = "ZM")

ETH_PAs <- st_transform(ETH_PAs, crs = 32637) # UTM Zone 37N for Ethiopia

# Check for invalid geometries
invalid_geom <- ETH_PAs %>% 
  filter(!st_is_valid(ETH_PAs))

# If there are invalid geometries, make them valid
if(nrow(invalid_geom) > 0) {
  ETH_PAs <- st_make_valid(ETH_PAs)
}


######################################################################################
#get percentage overlap of each ecoreiong in a PA

selected_ecoregions <- st_transform(selected_ecoregions, st_crs(ETH_PAs))
selected_ecoregions <- st_make_valid(selected_ecoregions)

ETH_PAs_union <- st_union(ETH_PAs)

selected_ecoregions <- selected_ecoregions %>%
  mutate(total_area = as.numeric(st_area(geometry)))

# Intersect the ecoregions with the protected areas
intersected <- st_intersection(selected_ecoregions, ETH_PAs_union)

# Calculate the area of the intersected (overlapping) regions (convert to numeric)
intersected <- intersected %>%
  mutate(intersect_area = as.numeric(st_area(geometry)))

# Summarize the total area and the intersected area for each ecoregion (ECO_NAME)
proportion_covered <- intersected %>%
  group_by(ECO_NAME) %>%
  summarise(
    total_intersect_area = sum(intersect_area),
    total_area = unique(total_area)  # Each ecoregion has one total area value
  ) %>%
  mutate(proportion_covered = 100*(total_intersect_area / total_area))

# View the results
print(proportion_covered)

# Calculate the total area for "Horn of Africa xeric bushlands"
horn_of_africa_area <- selected_ecoregions %>%
  filter(ECO_NAME == "Horn of Africa xeric bushlands") %>%
  summarise(total_area = as.numeric(st_area(geometry)))

# Create a new row with the calculated area and 0 for intersected area
new_row <- tibble(
  ECO_NAME = "Horn of Africa xeric bushlands",
  total_area = horn_of_africa_area$total_area,
  total_intersect_area = 0,  # since there is no intersection
  proportion_covered = 0       # proportion is 0 since it doesn't intersect
)

# Combine this new row with the existing proportion_covered data frame
proportion_covered <- bind_rows(proportion_covered, new_row)

Ecoregion_represent <- as.data.frame(proportion_covered)
Ecoregion_represent$ECO_NAME_clean <- c("Djibouti xeric shrublands", "East Sudanian savanna (ESS)",
                                        "Ethiopian montane forests", "Ethiopian montane grasslands and woodlands (EMG)",
                                        "Ethiopian montane moorlands (EMM)", "Masai xeric grasslands and shrublands",
                                        "Northern Acacia-Commiphora bushlands and thickets (NAC)", "Sahelian Acacia savanna",
                                        "Somali Acacia-Commiphora bushlands and thickets (SAC)","Saharan flooded grasslands (SFG)",
                                        "Horn of Africa xeric bushlands")

Ecoregion_represent$colour <- c("#D7D79E", "#FFA77F", "#267300", "#70A800", 
                                "#BED2FF", "#FFFFBE", "#98E600", "#FCC369",
                                "#FCD135", "#BEFFE8", "#FFAA00")

mean(Ecoregion_represent$proportion_covered, na.rm = TRUE)

###########################################################################################################

map_cols = Ecoregion_represent$colour

ECO_NAME_linebreak <- c("Saharan flooded grasslands",
                        "Ethiopian montane moorlands",
                        "Ethiopian montane forests\n(EMF)",
                        "Djibouti xeric shrublands",
                        "Masai xeric grasslands and\nshrublands",
                        "East Sudanian savanna",
                        "Sahelian Acacia savanna",
                        "Somali Acacia-Commiphora\nbushlands and thickets (SAC)",
                        "Ethiopian montane grasslands\nand woodlands (EMG)",
                        "Northern Acacia-Commiphora\nbushlands and thickets",
                        "Horn of Africa xeric bushlands\n(HAB)")

summary(Ecoregion_represent$proportion_covered)

plot <- Ecoregion_represent %>%
  ggplot() +
  labs(title = "",
       x = "",
       y = "Proportion of ecoregion area protected (%)") + 
  geom_bar(aes(x = reorder(as.factor(ECO_NAME_clean), -proportion_covered), y = proportion_covered), stat = "identity", show.legend = T, fill = map_cols, col = c("black", "black", "red", "red", "black","black","black","black","red", "black", "black")) +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(colour=c("black", "black", "red", "black", "black","black","black","red","red", "black", "red"),
                                   hjust = 1), axis.title.x = element_text(size = 12))+
  scale_y_continuous(limits = c(0,105), expand = c(0,0), breaks = seq(0, 105, by = 10))+
  geom_hline(yintercept = 9.4, linetype = "dashed")+
  annotate("text", x = 11.3, y = 12.6, label = "9.4% \ncurrent") +
  geom_hline(yintercept = 17, linetype = "dashed") +
  annotate("text", x=11.3, y=20, label="17% \nAichi") +
  geom_hline(yintercept = 30, linetype = "dashed")+
  annotate("text", x=11.3, y=33, label="30% \nGBF")  + 
  scale_x_discrete(labels = c (ECO_NAME_linebreak))
plot

plot <- Ecoregion_represent %>%
  ggplot() +
  labs(title = "",
       x = "",
       y = "Proportion of ecoregion area protected (%)") + 
  geom_bar(aes(x = reorder(as.factor(ECO_NAME_clean), -proportion_covered), y = proportion_covered), stat = "identity", show.legend = T, fill = map_cols, col = c("black", "black", "red", "red", "black","black","black","black","red", "black", "black")) +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(axis.text.y = element_text(colour=c("black", "black", "red", "black", "black","black","black","red","red", "black", "red"),
                                   hjust = 1), axis.title.x = element_text(size = 12))+
  scale_y_continuous(limits = c(0,105), expand = c(0,0), breaks = seq(0, 105, by = 10))+
  geom_hline(yintercept = 9.4, linetype = "dashed")+
  geom_hline(yintercept = 17, linetype = "dashed") +
  geom_hline(yintercept = 30, linetype = "dashed")+
  scale_x_discrete(labels = c (ECO_NAME_linebreak))
plot


ggsave("WWF_ecoregion_representativeness_9125.png", plot = plot, dpi=300, h=6, w = 9)


##############################################################################################
clipped_ecoregions$area_within_eth <- st_area(clipped_ecoregions)

# Step 1: Identify unique ecoregions in clipped_ecoregions
unique_ecoregions_in_eth <- unique(clipped_ecoregions$ECO_NAME)

# Step 2: Filter the global ecoregions dataset to include only those present in Ethiopia
filtered_ecoregions <- ecoregion_2017 %>%
  filter(ECO_NAME %in% unique_ecoregions_in_eth)

# Step 3: Calculate the global area for only those filtered ecoregions
filtered_ecoregions$global_area <- st_area(filtered_ecoregions)

# Step 4: Join the global area information back to the clipped_ecoregions dataset
clipped_ecoregions_df <- as.data.frame(clipped_ecoregions)
filtered_ecoregions_df <- as.data.frame(filtered_ecoregions)

clipped_ecoregions_df <- clipped_ecoregions_df %>%
  left_join(filtered_ecoregions_df %>% dplyr::select(ECO_NAME, global_area), by = "ECO_NAME")

# Step 5: Calculate the proportion of the global area within Ethiopia for each ecoregion
clipped_ecoregions_df$proportion_in_eth <- 100*(clipped_ecoregions_df$area_within_eth / clipped_ecoregions_df$global_area)

# View the results
clipped_ecoregions_df <- subset(clipped_ecoregions_df, select= c(ECO_NAME, proportion_in_eth, NNH_NAME))

nature_imperiled <-  clipped_ecoregions_df %>%
  filter(NNH_NAME == "Nature Imperiled")

nature_imperiled <- nature_imperiled %>%
  mutate(ECO_NAME = case_when(
    ECO_NAME == "Ethiopian montane forests" ~ "EMF",
    ECO_NAME == "Ethiopian montane grasslands and woodlands" ~ "EMG",
    ECO_NAME == "Horn of Africa xeric bushlands" ~ "HAB",
    ECO_NAME == "Somali Acacia-Commiphora bushlands and thickets" ~ "SAC",
    TRUE ~ ECO_NAME  # keep the original value if no match
  ))
levels = c("EMF", "SAC", "EMG", "HAB")
nature_imperiled_cols <- c("#267300","#70A800","#FFAA00","#FCD135" )

nature_imperiled$proportion_in_eth <- as.numeric(nature_imperiled$proportion_in_eth)


plot2 <- nature_imperiled %>%
  ggplot() +
  labs(title = "",
       x = "",
       y = "Global extent in Ethiopia (%)") + 
  geom_bar(aes(x = ECO_NAME, y = proportion_in_eth), stat = "identity", show.legend = F, fill = nature_imperiled_cols, col= "red") +
  coord_flip() +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  scale_x_discrete(limits = levels) +
  theme_classic(base_size = 10)+
  theme(axis.text.y = element_text(colour = c("black"), size = 8), axis.text.x = element_text(size = 8, colour = "black"),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(color = "red", fill = "transparent", linewidth = .6), panel.background = element_rect(fill = "transparent"),
        plot.margin = margin(-0.2, 0.2, 0.2, -0.2, unit = "cm"))

plot2
ggsave("WWF_ecoregion_global_extent_eth_141024.png", plot = plot2, dpi=300)


library(patchwork)
comb_plot <- plot + 
  inset_element(p = plot2,
                left = 0.61,
                right = 0.979,
                top = 0.45,
                bottom = 0.02)
#inset_element(legend, left = 0.1, bottom = 0.6, right = 0.4, top = 1, align_to = "full")
comb_plot
ggsave("PA_ecoregion_2017_representativeness_9125.png", plot = comb_plot, width = 9, height = 6, dpi = 300)

####################################################################################
#map plot
selected_ecoregions$colour <- c("#D7D79E", "#FFA77F", "#267300", "#70A800", 
                                "#BED2FF", "#FCBE51","#FFFFBE", "#98E600", "#FFAA00",
                                "#FCD135", "#BEFFE8")

st_write(selected_ecoregions, "Eth_2017_ecoregions_141024.shp")

cols <- selected_ecoregions$colour 

library(ggspatial)
# Plot selected ecoregions with colors and protected areas outlines in black
map_plot <- ggplot() +
  # Plot the ecoregions with custom fill colors (use your "map_cols" color vector)
  geom_sf(data = selected_ecoregions, aes(fill = ECO_NAME), color = NA) +
  # Overlay the protected areas with black outlines
  geom_sf(data = ETH, fill = NA, color = "black", size = 0.7) +
  # Custom color scale for ecoregions
  scale_fill_manual(values = cols) +
  # Adjust the theme
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "transparent"),
        panel.background = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank()) 


# Show the plot
map_plot


agri_potential <- raster("Agri_potential_aggregated_10m_average_10.tif")
agri_potential_crop <- crop(agri_potential, ETH)

agri_diversity <- raster("act_D.tif")
agri_diversity_crop <- crop(agri_diversity, ETH)

# Mask the cropped raster using the shapefile boundary
Eth_agri_potential <- mask(agri_potential_crop, ETH)
Eth_agri_diversity <- mask(agri_diversity, ETH)

plot(Eth_agri_potential)
plot(Eth_agri_diversity)
