########################################
# Ethiopia PA Ecoregion Representativeness Analysis
# Authors: Sophie Jago
# Description: Assess representativeness of Ethiopia's protected area network across ecoregions 
########################################
# =============================================================================
# Ecoregion Representation in Ethiopia's Protected Areas
# =============================================================================

# --- Load Libraries ---
library(ggplot2)
library(dplyr)
library(sf)
library(raster)
library(geodata)
library(patchwork)

# =============================================================================
# 1. Load and Prepare Spatial Data
# =============================================================================

## --- Ethiopia Boundary ---
ETH <- gadm(country = "ETH", level = 0, path = tempdir()) %>%
  st_as_sf()

## --- Ecoregions 2017 ---
ecoregion_2017 <- st_read("Ecoregions2017.shp")

## --- Ensure Matching CRS and Valid Geometries ---
if (!identical(st_crs(ETH), st_crs(ecoregion_2017))) {
  ecoregion_2017 <- st_transform(ecoregion_2017, st_crs(ETH))
}
ecoregion_2017 <- st_make_valid(ecoregion_2017)

## --- Clip Ecoregions to Ethiopia ---
clipped_ecoregions <- st_intersection(ecoregion_2017, ETH)

## --- Select Relevant Fields ---
selected_ecoregions <- clipped_ecoregions %>%
  select(ECO_NAME, COLOR) %>%
  st_make_valid()

# =============================================================================
# 2. Load and Prepare Protected Areas Data
# =============================================================================

ETH_PAs <- st_read("Eth_PAs_sep24.shp") %>%
  st_zm(drop = TRUE) %>%
  st_transform(crs = 32637) %>%
  st_make_valid()

# =============================================================================
# 3. Calculate Ecoregion Coverage in PAs
# =============================================================================

selected_ecoregions <- st_transform(selected_ecoregions, st_crs(ETH_PAs)) %>%
  mutate(total_area = as.numeric(st_area(geometry)))

ETH_PAs_union <- st_union(ETH_PAs)

## --- Intersection and Area Calculation ---
intersected <- st_intersection(selected_ecoregions, ETH_PAs_union) %>%
  mutate(intersect_area = as.numeric(st_area(geometry)))

## --- Summarise by Ecoregion ---
proportion_covered <- intersected %>%
  group_by(ECO_NAME) %>%
  summarise(
    total_intersect_area = sum(intersect_area),
    total_area = unique(total_area)
  ) %>%
  mutate(proportion_covered = 100 * total_intersect_area / total_area)

## --- Add Non Intersecting Ecoregion ---
horn_area <- selected_ecoregions %>%
  filter(ECO_NAME == "Horn of Africa xeric bushlands") %>%
  summarise(total_area = as.numeric(st_area(geometry)))

missing_row <- tibble(
  ECO_NAME = "Horn of Africa xeric bushlands",
  total_area = horn_area$total_area,
  total_intersect_area = 0,
  proportion_covered = 0
)

proportion_covered <- bind_rows(proportion_covered, missing_row)

# =============================================================================
# 4. Prepare Data for Plotting
# =============================================================================

Ecoregion_represent <- proportion_covered %>%
  as.data.frame() %>%
  mutate(
    ECO_NAME_clean = c("Djibouti xeric shrublands", "East Sudanian savanna (ESS)",
                       "Ethiopian montane forests", "Ethiopian montane grasslands and woodlands (EMG)",
                       "Ethiopian montane moorlands (EMM)", "Masai xeric grasslands and shrublands",
                       "Northern Acacia-Commiphora bushlands and thickets (NAC)", "Sahelian Acacia savanna",
                       "Somali Acacia-Commiphora bushlands and thickets (SAC)", "Saharan flooded grasslands (SFG)",
                       "Horn of Africa xeric bushlands"),
    colour = c("#D7D79E", "#FFA77F", "#267300", "#70A800", "#BED2FF", "#FFFFBE",
               "#98E600", "#FCC369", "#FCD135", "#BEFFE8", "#FFAA00")
  )

mean(Ecoregion_represent$proportion_covered, na.rm = TRUE)

# =============================================================================
# 5. Bar Plot: % Ecoregion Covered by PAs
# =============================================================================

ECO_NAME_linebreak <- c(
  "Saharan flooded grasslands", "Ethiopian montane moorlands",
  "Ethiopian montane forests\n(EMF)", "Djibouti xeric shrublands",
  "Masai xeric grasslands and\nshrublands", "East Sudanian savanna",
  "Sahelian Acacia savanna", "Somali Acacia-Commiphora\nbushlands and thickets (SAC)",
  "Ethiopian montane grasslands\nand woodlands (EMG)",
  "Northern Acacia-Commiphora\nbushlands and thickets",
  "Horn of Africa xeric bushlands\n(HAB)"
)

plot1 <- Ecoregion_represent %>%
  ggplot(aes(x = reorder(ECO_NAME_clean, -proportion_covered), y = proportion_covered)) +
  geom_bar(stat = "identity", fill = Ecoregion_represent$colour,
           col = c("black", "black", "red", "red", "black", "black", "black", "black", "red", "black", "black")) +
  coord_flip() +
  labs(y = "Proportion of ecoregion area protected (%)", x = "") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.y = element_text(color = c("black", "black", "red", "black", "black", "black", "black", "red", "red", "black", "red")),
    axis.title.x = element_text(size = 12)
  ) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 105, 10), expand = c(0, 0)) +
  scale_x_discrete(labels = ECO_NAME_linebreak) +
  geom_hline(yintercept = c(9.4, 17, 30), linetype = "dashed")

ggsave("WWF_ecoregion_representativeness.png", plot1, dpi = 300, height = 6, width = 9)

# =============================================================================
# 6. Plot: % Global Ecoregion Area in Ethiopia
# =============================================================================

clipped_ecoregions$area_within_eth <- st_area(clipped_ecoregions)

ecoregion_subset <- ecoregion_2017 %>%
  filter(ECO_NAME %in% unique(clipped_ecoregions$ECO_NAME)) %>%
  mutate(global_area = st_area(.))

clipped_ecoregions_df <- as.data.frame(clipped_ecoregions) %>%
  left_join(st_drop_geometry(ecoregion_subset) %>% select(ECO_NAME, global_area), by = "ECO_NAME") %>%
  mutate(proportion_in_eth = 100 * area_within_eth / global_area)

nature_imperiled <- clipped_ecoregions_df %>%
  filter(NNH_NAME == "Nature Imperiled") %>%
  mutate(ECO_NAME = recode(ECO_NAME,
                           "Ethiopian montane forests" = "EMF",
                           "Ethiopian montane grasslands and woodlands" = "EMG",
                           "Horn of Africa xeric bushlands" = "HAB",
                           "Somali Acacia-Commiphora bushlands and thickets" = "SAC"))

nature_cols <- c("#267300", "#FCD135", "#70A800", "#FFAA00")
levels <- c("EMF", "SAC", "EMG", "HAB")

plot2 <- ggplot(nature_imperiled, aes(x = ECO_NAME, y = as.numeric(proportion_in_eth))) +
  geom_bar(stat = "identity", fill = nature_cols, col = "red") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  scale_x_discrete(limits = levels) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    plot.background = element_rect(color = "red", fill = "transparent", linewidth = .6),
    panel.background = element_rect(fill = "transparent"),
    plot.margin = margin(-0.2, 0.2, 0.2, -0.2, "cm")
  ) +
  labs(y = "Global extent in Ethiopia (%)", x = "")

ggsave("WWF_ecoregion_global_extent_eth_141024.png", plot2, dpi = 300)

# =============================================================================
# 7. Combined Plot with Inset
# =============================================================================

comb_plot <- plot1 +
  inset_element(plot2, left = 0.61, right = 0.979, top = 0.45, bottom = 0.02)

ggsave("PA_ecoregion_2017_representativeness_9125.png", comb_plot, width = 9, height = 6, dpi = 300)


