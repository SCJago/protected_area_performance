# Load libraries
library(sf)
library(tidyverse)
library(giscoR)
library(rCAT)
library(ggpp)
library(DescTools)
library(ggplot2)
library(ggpp)
sf_use_s2(F)

# Set working directory
setwd("~/Papers/Ethiopia - PA history review/Representativeness/Species_representativeness/New_species_represenativeness_method")

########################### ---- IUCN RED LIST RANGE DATA ---- ###########################
# Read in the species data and filter PRESENCE and ORIGIN based on p. 204 of the KBA 
# guidelines
species_dat <- rbind(st_read("New_Data/IUCN_Red_List/data_0.shp"),
                     st_read("New_Data/IUCN_Red_List/data_1.shp"))

species_dat_prepared <- species_dat %>% 
  filter(PRESENCE %in% c(1, 2), ORIGIN %in% c(1, 2, 6)) %>% 
  dplyr::select(ID_NO, SCI_NAME, SUBSPECIES) %>% # retain relevant columns
  mutate(
    FULL_NAME = ifelse(
      is.na(SUBSPECIES) | str_detect(SCI_NAME, "ssp\\.|var\\."),
      SCI_NAME,  # If SUBSPECIES is NA or SCI_NAME contains "ssp." or "var.", keep SCI_NAME as is
      paste(SCI_NAME, SUBSPECIES)  # Otherwise, combine SCI_NAME with SUBSPECIES
    )
  ) %>% 
  mutate(
    # Remove exact matches of "ssp.", "subsp.", or "var." and any extra spaces around them
    FULL_NAME = str_replace_all(FULL_NAME, "\\s+(ssp\\.|subsp\\.|var\\.)\\s+", " "), # removes if between words
    FULL_NAME = str_replace_all(FULL_NAME, "\\s+(ssp\\.|subsp\\.|var\\.)$", ""),     # removes if at the end
    
    # Remove extra spaces (including leading/trailing spaces and multiple spaces between words)
    SCI_NAME = str_squish(SCI_NAME)
  ) %>% 
  st_transform(., crs = "EPSG:32637") # reproject to UTM Zone 37N 

# Identify and fix any invalid geometries
invalid_geom <- !st_is_valid(species_dat_prepared) # identify invalid geometries
table(invalid_geom, useNA = "ifany") # there is one NA value, which complicates things slightly as we have to account for that in the code below
na_index <- which(is.na(invalid_geom)) # find the index of the row(s) where invalid_geom is NA
species_dat_prepared <- species_dat_prepared[-na_index, ] # remove the row with the NA value from species_dat_prepared
invalid_geom <- invalid_geom[-na_index] # remove the NA element from invalid_geom
species_dat_prepared$geometry[invalid_geom] <- st_make_valid(species_dat_prepared$geometry[invalid_geom]) # fix only the invalid geometries

# Group by FULL_NAME and merge geometries while retaining the values of non-geometry columns
species_dat_prepared_union <- species_dat_prepared %>%
  group_by(FULL_NAME) %>%
  summarise(
    across(-geometry, first),        # apply `first()` to all non-geometry columns - in reality it wouldn't matter which you selected as the taxon and status columns should be the same for each record of the same species
    geometry = st_union(geometry),   # union the geometries for each species
  ) %>% 
  relocate(ID_NO) # put ID column back in the first position

# Read in and prepare the summary data
summary_dat <- read.csv("New_Data/IUCN_Red_List/simple_summary.csv") %>% 
  rename(SCI_NAME = scientificName, ID_NO = internalTaxonId) %>% # rename columns to align with polygon column names
mutate(
  FULL_NAME = ifelse(
    infraName == "" | str_detect(SCI_NAME, "ssp\\.|var\\.|subsp\\."),
    SCI_NAME,  # If infraName (subspecies) is NA or SCI_NAME contains "ssp.", "var." or "subsp.", keep SCI_NAME as is
    paste(SCI_NAME, infraName)  # Otherwise combine SCI_NAME with infraName
  )
) %>% 
  mutate(
    # Remove exact matches of "ssp.", "subsp.", or "var." and any extra spaces around them
    FULL_NAME = str_replace_all(FULL_NAME, "\\s+(ssp\\.|subsp\\.|var\\.)\\s+", " "), # removes if between words
    FULL_NAME = str_replace_all(FULL_NAME, "\\s+(ssp\\.|subsp\\.|var\\.)$", ""), # removes if at the end
    
    # Remove extra spaces (including leading or trailing spaces and multiple spaces between words)
    SCI_NAME = str_squish(SCI_NAME)
  ) %>% 
  dplyr::select(ID_NO, FULL_NAME, kingdomName, className, redlistCategory)

# Combine the species data with the summary data by the mutual "FULL_NAME" column
name_matched_dat <- left_join(species_dat_prepared_union, summary_dat, by = "FULL_NAME")

# If any rows have not matched successfully, match them by "ID_NO" instead
unmatched_dat <- name_matched_dat %>% 
  filter(is.na(redlistCategory))
unmatched_dat <- species_dat_prepared_union[species_dat_prepared_union$FULL_NAME %in% unmatched_dat$FULL_NAME,]

id_matched_dat <- left_join(unmatched_dat, summary_dat, by = "ID_NO")

matched_dat <- rbind(subset(name_matched_dat, !is.na(redlistCategory)) %>% 
                        dplyr::select(ID_NO.x, FULL_NAME, kingdomName, className, redlistCategory) %>% 
                        rename(ID_NO = "ID_NO.x"),
                      id_matched_dat %>% 
                        dplyr::select(ID_NO, FULL_NAME.x, kingdomName, className, redlistCategory) %>% 
                        rename(FULL_NAME = "FULL_NAME.x"))

# Save the prepared range data
write_sf(matched_dat, "New_Data/Output_Data/IUCN_Ranges.gpkg")

###################### ---- POINT DATA (PLANTS WITHOUT RANGES) ---- ######################
# Read the prepared IUCN Red List range data back in
matched_dat <- st_read("New_Data/Output_Data/IUCN_Ranges.gpkg")

# Read in and prepare the IUCN Red List points data (filtering for plants only)
species_RLpoint_dat <- read.csv("New_Data/IUCN_Red_List/points_data.csv") %>% 
  filter(presence %in% c(1, 2), origin %in% c(1, 2, 6)) %>% 
  mutate(
    FULL_NAME = ifelse(
      is.na(subspecies) | str_detect(sci_name, "ssp\\.|var\\.|subsp\\."),
      sci_name,  # If subspecies is NA or SCI_NAME contains "ssp." or "var.", keep sci_name as is
      paste(sci_name, subspecies)  # Otherwise, combine sci_name with subspecies
    )
  ) %>% mutate(
    # Remove exact matches of "ssp.", "subsp.", or "var." and any extra spaces around them
    FULL_NAME = str_replace_all(FULL_NAME, "\\s+(ssp\\.|subsp\\.|var\\.)\\s+", " "), # removes if between words
    FULL_NAME = str_replace_all(FULL_NAME, "\\s+(ssp\\.|subsp\\.|var\\.)$", ""),     # removes if at the end
    
    # Remove extra spaces (including leading/trailing spaces and multiple spaces between words)
    FULL_NAME = str_squish(FULL_NAME)
  ) %>% 
  left_join(., summary_dat, by = "FULL_NAME") %>% # add the summary data in
  filter(kingdomName == "PLANTAE") %>% # filter to only retain plants
  dplyr::select(ID_NO, FULL_NAME, kingdomName, className, redlistCategory, longitude, latitude) # retain relevant columns

# Read in and prepare the BRAHMS point data
species_BRAHMSpoint_dat <- read.csv("New_Data/BRAHMS/BRAHMS_Occurrence_Data.csv") %>% 
  left_join(., summary_dat, by = c("SPECIES" = "FULL_NAME")) %>% # add the summary data in
  drop_na(redlistCategory) %>% # remove species that were not matched and are thus assumed to be not evaluated on the IUCN Red List
  dplyr::select(ID_NO, SPECIES, kingdomName, className, redlistCategory, LONGDEC, LATDEC) %>% # retain relevant columns
  rename(longitude = "LONGDEC", latitude = "LATDEC", FULL_NAME = "SPECIES") # rename columns to match species_RLpoint_dat

# Combine the two point dataframes into one
point_dat <- rbind(species_RLpoint_dat, species_BRAHMSpoint_dat) %>% 
  distinct() # remove duplicate rows, i.e. those with identical values in all columns

# Remove occurrence records of species that already have an IUCN Red List range
index <- point_dat$FULL_NAME[point_dat$FULL_NAME %in% matched_dat$FULL_NAME]
point_dat <- point_dat %>%
  filter(!FULL_NAME %in% index)

# Generate Rapoport ranges for species with 3 or more occurrences and add buffers to species
# with 1/2 occurences
species_with_point_data <- unique(point_dat$FULL_NAME) # vector of plant species names that DO NOT have ranges but DO have occurrence data

rap_list <- list()
for (i in seq_along(species_with_point_data)) {
  # Select rows in turn for each species 
  species <- point_dat %>% 
    filter(FULL_NAME == species_with_point_data[i])
  
  # Create a dataframe with the latitude and longitude of each distinct occurrence
  ll <- data.frame(lat = species$latitude,
                   long = species$longitude)
  
  # If the species has 3 or more distinct occurrences, generate Rapoport ranges
  if (nrow(ll) >= 3) {
    thepoints <- simProjWiz(ll, returnV = "S")
    # sfs <- subLocRapoport(thepoints, bufferDis = 5000, returnV = "SF") # uncomment this to manually assign a 5km buffer distance
    sfs <- subLocRapoport(thepoints, returnV = "SF") # uncomment this to assign the default buffer distance
    sfs <- sfs$buffers %>% 
      st_transform(., crs = "EPSG:32637") %>% # reproject to UTM Zone 37N 
      st_as_sf() %>% 
      mutate(ID_NO = species$ID_NO[1],
             FULL_NAME = species$FULL_NAME[1],
             kingdomName = species$kingdomName[1],
             className = species$className[1],
             redlistCategory = species$redlistCategory[1]) %>% 
      rename(geometry = x)
  }
  
  # If the species has only 1/2 distinct occurrences, add a buffer around each point
  else {
    sfs <- species %>% 
      st_as_sf(., coords = c("longitude", "latitude"), crs = "EPSG:4326") %>% # add spatial features to the dataframe 
      st_transform(., crs = "EPSG:32637") # reproject to UTM Zone 37N 
    sfs <- st_buffer(sfs, dist = 5000) # add 5 km buffer
    
    # If there are 2 buffered points for the same species, dissolve their geometries into one
    if (nrow(sfs) > 1) {
      sfs <- sfs %>% 
        group_by(FULL_NAME) %>%
        summarise(
          across(-geometry, first),        # apply `first()` to all non-geometry columns - in reality it wouldn't matter which you selected as the columns should be the same for each record of the same species
          geometry = st_union(geometry),   # dissolve the geometries for the species
        ) %>% 
        relocate(ID_NO) %>%  # put ID column back in the first position
        st_as_sf()
    }
  }
  rap_list[[i]] <- sfs
}

# Combine each species' range into a dataframe and crop by Ethiopia's boundary
dat <- do.call(rbind, rap_list)

# Save the Rapoport and buffered range data
write_sf(dat, "New_Data/Output_Data/Rapoport_Ranges.gpkg") # to overwrite the 5km subLocRapoport method data
write_sf(dat, "New_Data/Output_Data/Rapoport_Ranges_auto.gpkg") # to overwrite the automatic/default subLocRapoport method data



############################### ---- OVERLAP ANALYSIS ---- ###############################
# Get Ethiopia's national boundary
eth <- gisco_get_countries(country = "Ethiopia", year = "2020", resolution = "01") %>% # country boundary at high resolution 
  dplyr::select(geometry) %>% 
  st_transform(., crs = "EPSG:32637") # reproject to UTM Zone 37N 

# Read in Ethiopia's protected areas
protected_areas <- st_read("New_Data/Protected_Areas/Eth_PAs_Updated/final _ETH_PAS.shp") %>% 
  st_zm(drop = T) %>% 
  st_transform(., crs = "EPSG:32637") %>% # reproject to UTM Zone 37N 
  dplyr::select(geometry)

# Find and fix invalid geometries in PA data
invalid_geom <- !st_is_valid(protected_areas) # identify invalid geometries
table(invalid_geom) 
protected_areas$geometry[invalid_geom] <- st_make_valid(protected_areas$geometry[invalid_geom]) # fix only the invalid geometries

# Read in and combine the IUCN range data and the Rapoport ranges we generated for plants
# without an IUCN Red List range
species_dat <- rbind(st_read("New_Data/Output_Data/IUCN_Ranges.gpkg"),
                     st_read("New_Data/Output_Data/Rapoport_Ranges_auto.gpkg"))

# Find and fix invalid geometries in species_dat
invalid_geom <- !st_is_valid(species_dat) # identify invalid geometries
table(invalid_geom, useNA = "ifany") # there are NA values
na_index <- which(is.na(invalid_geom)) # find the index of the rows where invalid_geom is NA
species_dat <- species_dat[-na_index, ] # remove rows with NA values from species_dat
invalid_geom <- invalid_geom[-na_index] # remove NA elements from invalid_geom
species_dat$geom[invalid_geom] <- st_make_valid(species_dat$geom[invalid_geom]) # fix only the invalid geometries

# Plot the data to make sure everything looks right (spoiler... it does!)
terra::plot(eth %>% as_Spatial())
terra::plot(species_dat[7,] %>% as_Spatial(), add = T, col = "green")
terra::plot(protected_areas %>% as_Spatial(), add = T, col = "red")

# Remove any species that have none of their range inside Ethiopia
species_dat <- species_dat %>% 
  filter(st_intersects(., eth, sparse = FALSE))

# Add two new columns, one for taxon and one for binary threat status
species_dat <- species_dat %>% 
  mutate(
    threat_status = case_when(
      redlistCategory %in% c("Critically Endangered", "Endangered", "Vulnerable") ~ "Threatened",
      redlistCategory %in% c("Least Concern", "Near Threatened", "Lower Risk/near threatened") ~ "Not Threatened"
    ),
    taxon = case_when(
      className %in% c("REPTILIA", "AMPHIBIA") ~ "Herptile",
      className == "MAMMALIA" ~ "Mammal",
      className == "AVES" ~ "Bird",
      kingdomName == "PLANTAE" ~ "Plant"
    )
  ) 

# View the number of species by threat status and taxon
species_dat %>% 
  st_drop_geometry() %>% 
  group_by(threat_status, taxon) %>% 
  summarise(n = n())

# Calculate the area of each species' global range
species_dat <- species_dat %>% 
  mutate(global_range_area_km2 = st_area(.) %>% 
           units::set_units(km2) %>% 
           units::drop_units())

# Calculate the area of each species' Ethiopian range
species_dat <- species_dat %>% 
  st_intersection(., eth) %>% # find intersection with Ethiopia
  mutate(national_range_area_km2 = st_area(.) %>% 
           units::set_units(km2) %>% 
           units::drop_units())

# Find and calculate the area of the part of each species' range that is inside PAs
unionised_protected_areas <- st_union(protected_areas)
protected_species_ranges <- st_intersection(species_dat, unionised_protected_areas) %>% 
  mutate(protected_range_area_km2 = st_area(.) %>% 
           units::set_units(km2) %>% 
           units::drop_units()) %>% 
  st_drop_geometry() %>% 
  dplyr::select(FULL_NAME, protected_range_area_km2)

# Join protected_species_ranges to species_dat. If the species has a missing protected_range_area_km2
# value, i.e. the species does not fall in any PA, change the value to 0
species_dat <- left_join(species_dat, protected_species_ranges, by = "FULL_NAME") %>%
  st_drop_geometry() %>% 
  mutate(protected_range_area_km2 = ifelse(is.na(protected_range_area_km2), 0, protected_range_area_km2))


# Find the proportion of each species' total national range area that is protected, and
# the proportion of each species' global range that is inside Ethiopia
species_dat <- species_dat %>% 
  mutate(national_global_pct = national_range_area_km2/global_range_area_km2*100,
         range_protection_pct = protected_range_area_km2/national_range_area_km2*100)

# Save the results data
write.csv(species_dat, "New_Data/Output_Data/species_dat.csv", row.names = FALSE) # to overwrite the 5km subLocRapoport method data
write.csv(species_dat, "New_Data/Output_Data/species_dat_auto_rapoport.csv", row.names = FALSE) # to overwrite the automatic/default subLocRapoport method data

############################## ---- STATISTICAL TESTS ---- ###############################
# Read in the prepared output data
species_dat <- read.csv("New_Data/Output_Data/species_dat_auto_rapoport.csv")

# 1) Kruskal-Wallis test: differences in coverage of all species by PAs, in regards to taxon
out_all <- kruskal.test(range_protection_pct ~ taxon, 
                        data = species_dat)
print(out_all) # H = 37.119, df = 3, p < 0.001

DunnTest(range_protection_pct ~ taxon, 
         data = species_dat,
         method = "bonferroni")

# Dunn's test of multiple comparisons using rank sums : bonferroni  
# 
#                 mean.rank.diff    pval    
# Herptile-Bird       -122.14539 0.03344 *  
# Mammal-Bird          163.93174 0.00057 ***
# Plant-Bird           -58.79120 0.31387    
# Mammal-Herptile      286.07713 3.4e-07 ***
# Plant-Herptile        63.35419 0.89634    
# Plant-Mammal        -222.72294 6.3e-07 ***

# 2) Kruskal-Wallis test: differences in coverage of threatened species by PAs, in regards to taxon
out_threatened <- kruskal.test(range_protection_pct ~ taxon, 
                               data = subset(species_dat, threat_status == "Threatened"))
print(out_threatened) # H = 51.69, df = 3, p < 0.001

DunnTest(range_protection_pct ~ taxon, 
         data = subset(species_dat, threat_status == "Threatened"),
         method = "bonferroni")

# Dunn's test of multiple comparisons using rank sums : bonferroni  
# 
#                 mean.rank.diff    pval    
# Herptile-Bird         6.708333 1.00000    
# Mammal-Bird          15.708333 1.00000    
# Plant-Bird          -65.206174 0.00012 ***
# Mammal-Herptile       9.000000 1.00000    
# Plant-Herptile      -71.914508 0.00886 ** 
# Plant-Mammal        -80.914508 9.1e-09 ***

################################### ---- FIGURES ---- ####################################
# Plot the main figure
ggplot(species_dat) +
  geom_boxplot(aes(x = range_protection_pct, 
                   y = factor(taxon, level = c('Plant', 'Bird', 'Mammal', 'Herptile')),
                   fill = threat_status),
               outlier.shape = NA) +
  geom_point(data = subset(species_dat, redlistCategory != "Critically Endangered"),
             aes(x = range_protection_pct,
                 y = factor(taxon, level = c('Plant', 'Bird', 'Mammal', 'Herptile')),
                 shape = threat_status),
             position = position_jitterdodge(jitter.width = 0.3, seed = 002),
             alpha = 0.6, size = 1)+
  geom_point(data = subset(species_dat, redlistCategory == "Critically Endangered"),
             aes(x = range_protection_pct,
                 y = factor(taxon, level = c('Plant', 'Bird', 'Mammal', 'Herptile')),
                 color = "red"),
             position = position_jitternudge(height = 0.16, y = 0.19, nudge.from = "jittered.y",
                                             seed = 010),
             size = 1.25)+
  scale_shape_manual(values = c(17, 16))+
  scale_fill_manual(values = c("grey90", "grey50")) +
  scale_color_manual(values = c("red")) +
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0.2),breaks = seq(0, 105, by = 10)) +
  theme_classic(base_size = 9) +
  labs(x = "Percentage of range protected",
       y = NULL) +
  theme(axis.text.y = element_text(colour = c("black")),
        # axis.line.y = element_blank(),
        legend.spacing = unit(0.105, "cm"),
        legend.margin = margin(0,0.15,0.17,0.15, "cm"),
        legend.position = c(.6089, .378),
        legend.background = element_rect(fill = NULL, color = "black", linewidth = 0.55))+
  guides(shape = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         fill = guide_legend(title = NULL)) 

# Create a new dataframe for just CR species - this will be used for the inset figure
critically_endangered_species <- subset(species_dat, 
                                        redlistCategory == "Critically Endangered")

# View number of CR species in each taxonomic group
critically_endangered_species %>% 
  group_by(taxon) %>% 
  summarise(n = n())

# Count rows where kingdomName is "PLANTAE" and range_protection_pct is 0
num_rows_CE <- nrow(subset(critically_endangered_species, kingdomName == "PLANTAE" & range_protection_pct == 0))
print(num_rows_CE)

# Plot the inset figure
critically_endangered_species %>% 
  ggplot(aes(x = national_global_pct, 
             y = factor(taxon, level = c('Plant', 'Bird', 'Mammal', 'Herptile')),
             fill = taxon))+
  labs(x = "Mean global extent in Ethiopia (%)",
       y = "Critically Endangered species") +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, show.legend = F)+
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) + 
  scale_y_discrete(labels = c('Plants\n(n = 31)', 'Birds\n(n = 7)', 'Mammals\n(n = 5)',
                              'Herptiles\n(n = 2)'))+
  theme_classic(base_size = 9)+
  theme(axis.text.y = element_text(colour = c("black")))+
  theme(plot.background = element_rect(color = "black", fill = NA, linewidth = 1))

############################# SOPHIE FIGURES ##############################################

p1 <- ggplot() +
  geom_boxplot(data = species_dat,
               aes(x = range_protection_pct, 
                   y = factor(taxon, level = c('Plant', 'Bird', 'Mammal', 'Herptile')),
                   fill = threat_status),
               outlier.shape = NA) +
  geom_point(data = subset(species_dat, redlistCategory != "Critically Endangered"),
             aes(x = range_protection_pct,
                 y = factor(taxon, level = c('Plant', 'Bird', 'Mammal', 'Herptile')),
                 shape = threat_status),
             position = position_jitterdodge(jitter.width = 0.2, seed = 002),
             alpha = 0.3, size = 1.25) +
  geom_point(data = subset(species_dat, redlistCategory == "Critically Endangered"),
             aes(x = range_protection_pct,
                 y = factor(taxon, level = c('Plant', 'Bird', 'Mammal', 'Herptile')),
                 color = "Critically Endangered"),
             position = position_jitternudge(height = 0.1, y = 0.19, nudge.from = "jittered.y",
                                             seed = 010),
             size = 1.5)+
  scale_shape_manual(values = c(16, 17), breaks = c("Threatened", "Not Threatened"))+
  scale_fill_manual(values = c("grey50", "grey90"), breaks = c("Threatened", "Not Threatened")) +
  scale_color_manual(values = c("Critically Endangered" = "red")) +  # Updated for correct legend label
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0.2),breaks = seq(0, 105, by = 10)) +
  scale_y_discrete(labels = c( "       ", "       ", "       ", "       "))+
  theme_classic(base_size = 14) +
  labs(x = "Percentage of range protected",
       y = NULL) +
  theme(axis.text.y = element_text(colour = c("black")), axis.title.x  = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.key.width = unit(1.5, "cm"),
        legend.spacing = unit(0, "cm"),
        legend.margin = margin(0,0,0,0, "cm"),
        legend.position = c(.8, .1),
        legend.background = element_rect(fill = "transparent"))+
  guides(shape = guide_legend(title = NULL, override.aes = list(size = 3)),
         color = guide_legend(title = NULL, override.aes = list(size = 3)),
         fill = guide_legend(title = NULL)) +
  geom_rect(aes(xmin = 67, xmax = 103, ymin = 0.5, ymax = 1.15),
            color = "black", alpha = 0, fill = "grey", size = 0.3) 

p1

#CREATE AN INSET FIGURE SHOWING THE MEAN GLOBAL EXTENT OF CRITICALLY ENDANGERED
#SPECIES IN ETHIOPIA
p2 <- critically_endangered_species %>% 
  ggplot(aes(x = national_global_pct, y = factor(taxon, level = c('Plant', 'Bird', 'Mammal', 'Herptile')), fill = taxon,
  ))+
  labs(x = "Mean proportion of range\n found within Ethiopia (%)",
       y = "Critically endangered\nspp. subset") +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, show.legend = F, fill = c("#009E73", "#0072B2", "#E69F00", "#CC79A7"))+
  scale_x_continuous(limits = c(0, 105), expand = c(0, 0)) + 
  scale_y_discrete(labels = c('        (n=31)', '        (n=7)', '        (n=5)',
                              '        (n=2)'))+
  theme_classic(base_size = 11)+
  theme(axis.text.y = element_text(colour = c("black"), size = 8), axis.text.x = element_text(size = 11, colour = "black"),
        axis.ticks.y = element_blank(),
        plot.background = element_rect(color = "red", fill = "transparent", linewidth = .6), panel.background = element_rect(fill = "transparent"))

p2

#ADD THE TWO FIGURES TOGETHER
library(patchwork)
p3 <- p1 +
  inset_element(p = p2,
                left = 0.636,
                right = 0.979,
                top = 0.63,
                bottom = 0.21)

p3

ggsave("ETH_species_representativeness9125.png", plot = p3, width = 9, height = 6, dpi = 300)



############################## ---- SPECIES IN EACH PA ---- ###############################
# Get Ethiopia's national boundary
eth <- gisco_get_countries(country = "Ethiopia", year = "2020", resolution = "01") %>% # country boundary at high resolution 
  dplyr::select(geometry) %>% 
  st_transform(., crs = "EPSG:32637") # reproject to UTM Zone 37N 

# Read in Ethiopia's protected areas
protected_areas <- st_read("New_Data/Protected_Areas/Eth_PAs_Updated/final _ETH_PAS.shp") %>% 
  st_zm(drop = T) %>% 
  st_transform(., crs = "EPSG:32637") %>% # reproject to UTM Zone 37N 
  group_by(ORIG_NAME) %>%
  summarize(geometry = st_union(geometry)) # merge the geometries of separate rows with same PA name
  
# Find and fix invalid geometries in PA data
invalid_geom <- !st_is_valid(protected_areas) # identify invalid geometries
table(invalid_geom) 
protected_areas$geometry[invalid_geom] <- st_make_valid(protected_areas$geometry[invalid_geom]) # fix only the invalid geometries

# Read in and combine the IUCN range data and the Rapoport ranges we generated for plants
# without an IUCN Red List range
species_dat <- rbind(st_read("New_Data/Output_Data/IUCN_Ranges.gpkg"),
                     st_read("New_Data/Output_Data/Rapoport_Ranges_auto.gpkg"))

# Find and fix invalid geometries in species_dat
invalid_geom <- !st_is_valid(species_dat) # identify invalid geometries
table(invalid_geom, useNA = "ifany") # there are NA values
na_index <- which(is.na(invalid_geom)) # find the index of the rows where invalid_geom is NA
species_dat <- species_dat[-na_index, ] # remove rows with NA values from species_dat
invalid_geom <- invalid_geom[-na_index] # remove NA elements from invalid_geom
species_dat$geom[invalid_geom] <- st_make_valid(species_dat$geom[invalid_geom]) # fix only the invalid geometries

# Remove any species that have none of their range inside Ethiopia
species_dat <- species_dat %>% 
  filter(st_intersects(., eth, sparse = FALSE))

# Add two new columns, one for taxon and one for binary threat status
species_dat <- species_dat %>% 
  mutate(
    threat_status = case_when(
      redlistCategory %in% c("Critically Endangered", "Endangered", "Vulnerable") ~ "Threatened",
      redlistCategory %in% c("Least Concern", "Near Threatened", "Lower Risk/near threatened") ~ "Not Threatened"
    ),
    taxon = case_when(
      className %in% c("REPTILIA", "AMPHIBIA") ~ "Herptile",
      className == "MAMMALIA" ~ "Mammal",
      className == "AVES" ~ "Bird",
      kingdomName == "PLANTAE" ~ "Plant"
    )
  ) 

# Create a spatial intersection matrix with TRUE/FALSE values depending on whether a given
# species is inside a given PA. Each row will correspond to a PA, and each column to a species
intersects_matrix <- st_intersects(protected_areas, species_dat, sparse = FALSE)

# Apply a function to each protected area, iterating over the range from 1 to the number of 
# rows in protected_areas. `pa_index` represents the current protected_area row being processed.
pa_species <- lapply(1:nrow(protected_areas), function(pa_index) {
  
  # For the current PA, find the species ranges that intersect it
  intersected_species <- which(intersects_matrix[pa_index, ])
  
  if (length(intersected_species) > 0) {
    
    # Create a data frame for the current PA
    data.frame(
      PA_ID = protected_areas$ORIG_NAME[pa_index], # PA's original name (ID)
      species_id = species_dat$FULL_NAME[intersected_species], # name of intersecting species
      taxon = species_dat$taxon[intersected_species], # taxon of intersecting species
      threat_status = species_dat$threat_status[intersected_species] # threat status of intersecting species
    )
    
  } else {
    NULL
  }
})

# Combine into a single dataframe
pa_species_df <- do.call(rbind, pa_species)

# Count the number of species in each PA by taxon and threat status
pa_species_df <- pa_species_df %>%
  group_by(PA_ID, taxon, threat_status) %>%
  summarise(species_count = n(), .groups = 'drop')

# Pivot wider so we have one row per PA and one column per taxon/threat_status
pa_species_df <- pa_species_df %>% 
  pivot_wider(names_from = c("taxon", "threat_status"), values_from = "species_count") %>% 
  mutate(across(2:9, ~ if_else(is.na(.), 0, .))) # assign 0 to NAs

# Save this count dataframe
write.csv(pa_species_df, "New_Data/Output_Data/n_species_per_pa.csv", row.names = FALSE)
