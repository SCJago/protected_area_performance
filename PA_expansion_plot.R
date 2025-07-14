library(ggplot2)
library(scales)
library(vistime) #to ggplot timeline data
library(egg) #for ggarrange
library(grid) #for plotting only legend
library(gridExtra) #for plotting only legend
library(ggprism) #for minor ticks
library(patchwork) #for graph inset
library(png) # bring in legend png file
library(sf)
library(dplyr)
library(raster)

setwd("~/Papers/Ethiopia - PA history review/EWCA_updated_PAs_May23/Expansion_PA_network")


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

# Calculate area in square meters and convert to square kilometers
ETH_PAs$area_km2 <- as.numeric(st_area(ETH_PAs) / 1e6)

# Check the data
head(ETH_PAs)

#bring in PA year
PA_year <- read.csv("Earliest_PA_year.csv")

# Join the shapefile data with the CSV data by matching the names
ETH_PAs <- ETH_PAs %>%
  left_join(PA_year, by = c("NAME_PAs" = "Name_PAs"))

print(ETH_PAs$NAME_PAs)

nrow(ETH_PAs)
num_IUCN_CAT_II <- ETH_PAs %>%
  filter(IUCN_MG_CA == "II") %>%
  nrow()

# Print the result
print(num_IUCN_CAT_II)

ETH_PAs <- ETH_PAs %>%
  mutate(IUCN_MG_CA = if_else(is.na(IUCN_MG_CA), "VI", IUCN_MG_CA))

num_IUCN_CAT_LS <- ETH_PAs %>%
  filter(IUCN_MG_CA != "II") %>%
  nrow()

# Print the result
print(num_IUCN_CAT_LS)

#get and save csv with PA names and info
Eth_PAs_supplement_table <- subset(ETH_PAs, select = c("NAME_PAs", "DESIGNATIO", "IUCN_MG_CA", "area_km2", "Earliest_year"))
Eth_PAs_supplement_table <- Eth_PAs_supplement_table[order(Eth_PAs_supplement_table$NAME_PAs), ]
Eth_PAs_supplement_table <- st_drop_geometry(Eth_PAs_supplement_table)
str(Eth_PAs_supplement_table)
write.csv(Eth_PAs_supplement_table, "Eth_PAs_supplement_table.csv")

# Initialize a dataframe to store the cumulative area and number of PAs for each year
cum_area_df <- data.frame(Year = sort(unique(ETH_PAs$Earliest_year)), Cum_Area = 0, Cum_Num = 0)

# Create an empty object to store the cumulative protected areas
cumulative_polygons <- st_sfc(crs = st_crs(ETH_PAs))  # Start with an empty spatial feature

# Initialize a counter for the cumulative number of PAs
cum_num <- 0

# Loop through each year, union the polygons, subtract overlaps, calculate area, and update PA count
for (i in seq_along(cum_area_df$Year)) {
  # Subset polygons established in the current year
  current_year_polygons <- ETH_PAs %>% 
    filter(Earliest_year == cum_area_df$Year[i])
  
  print(current_year_polygons$NAME_PAs)
  
  # Update the cumulative number of protected areas
  cum_num <- cum_num + nrow(current_year_polygons)  # Count the number of PAs in the current year
  print(cum_num)
  cum_area_df$Cum_Num[i] <- cum_num  # Store the cumulative number of PAs
  
  # Union the polygons for the current year (this merges any overlapping PAs within the same year)
  unioned_current_year <- st_union(current_year_polygons)
  
  # Subtract any overlapping area with previously protected areas (cumulative polygons)
  new_non_overlapping_area <- st_difference(unioned_current_year, cumulative_polygons)
  
  # Calculate the area of the new non-overlapping area (in km^2)
  new_area <- as.numeric(st_area(new_non_overlapping_area)) / 1e6
  
  # Add the new area to the cumulative polygons
  cumulative_polygons <- st_union(cumulative_polygons, unioned_current_year)
  
  # Update the cumulative area for the current year
  if (i == 1) {
    cum_area_df$Cum_Area[i] <- new_area  # First year, just the new area
  } else {
    cum_area_df$Cum_Area[i] <- cum_area_df$Cum_Area[i-1] + new_area  # Add new area to previous total
  }
  

}

# View the result
print(cum_area_df)

###########################################################################################

numberColor <- "black" #set colours for number line
areaColor <- "#999999" #set colours for area block
Eth_area <- 1127095 #constant area of Eth to get percentage PA coverage
cum_area_df$Percent_area <- (cum_area_df$Cum_Area/Eth_area)*100 #calculate percent coverage as new column

#make it go to 2024
new_row <- data.frame(
  Year = 2024,
  Cum_Area = cum_area_df$Cum_Area[cum_area_df$Year == 2018],
  Cum_Num = cum_area_df$Cum_Num[cum_area_df$Year == 2018],
  Percent_area = cum_area_df$Percent_area[cum_area_df$Year == 2018]
)

# Combine the new row with the existing dataframe
cum_area_df <- rbind(cum_area_df, new_row)

############################################################################################################################################################################
#create plot for PA expansion over time using area of PAs and number of PAs 
############################################################################################################################################################################
PA_expansion_percent_timeline <- ggplot(cum_area_df, aes(x=Year)) +
  
  geom_rect(aes(xmin = 1958, xmax = 1965, ymin = -10, ymax = -1),
            fill = "#999999") +
  geom_rect(aes(xmin = 1965, xmax = 1980, ymin = -10, ymax = -1),
            fill = "#E69F00")+
  geom_rect(aes(xmin = 1980, xmax = 1993, ymin = -10, ymax = -1),
            fill = "#56B4E9") +
  geom_rect(aes(xmin = 1993, xmax = 1995, ymin = -10, ymax = -1),
            fill = "#009E73") +
  geom_rect(aes(xmin = 1995, xmax = 1998, ymin = -10, ymax = -1),
            fill = "#F0E442") +
  geom_rect(aes(xmin = 1998, xmax = 2003, ymin = -10, ymax = -1),
            fill = "#0072B2")+
  geom_rect(aes(xmin = 2003, xmax = 2008, ymin = -10, ymax = -1),
            fill = "#D55E00") +
  geom_rect(aes(xmin = 2008, xmax = 2024, ymin = -10, ymax = -1),
            fill = "#CC79A7") +
  geom_area( aes(y=Percent_area*10), size=2, fill = areaColor) +   #the *10 is to put on same scale as other axis so they overlap properly
  geom_line( aes(y=Cum_Num), stat="identity", size=.8, color=numberColor, alpha=1) +
  geom_point( aes(y=Cum_Num), stat="identity", size=1, fill=numberColor, color = numberColor, alpha=1) +
  
  
  geom_segment(aes(x=1962, y=-22, xend=1962, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=1994, y=-13, xend=1994, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=2000, y=-22, xend=2000, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=2010, y=-13, xend=2010, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=2015, y=-22, xend=2015, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=2022, y=-13, xend=2022, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  
  geom_label(label="Partnered with UNESCO", x=1962, y=-22,label.padding = unit(0.3, "lines"), label.size = 0.35,color = "black",fill="lightgrey") +
  geom_label(label="Ratified CBD", x=1994, y=-13,label.padding = unit(0.3, "lines"), label.size = 0.35,color = "black",fill="lightgrey") +
  geom_label(label="Aligned with MDGs", x=2000, y=-22,label.padding = unit(0.3, "lines"), label.size = 0.35,color = "black",fill="lightgrey") +
  geom_label(label="Adopted Aichi targets", x=2008, y=-13,label.padding = unit(0.3, "lines"), label.size = 0.35,color = "black",fill="lightgrey") +
  geom_label(label="Aligned with SDGs", x=2015, y=-22,label.padding = unit(0.3, "lines"), label.size = 0.35,color = "black",fill="lightgrey") +
  geom_label(label="Adopted GBF targets", x=2019, y=-13,label.padding = unit(0.3, "lines"), label.size = 0.35,color = "black",fill="lightgrey") +
  
  scale_y_continuous(expand= c(0,0), limits = c(-29,105),
                     name = "               Cumulative number of PAs",
                     sec.axis = sec_axis(~./10, name="Cumulative PA coverage (%)               ", #here have to re-divide it by 10 to make axis reflect actual values
                                         breaks = seq(0,10,1)),
                     breaks = seq(0,100,10)) +
  
  theme_classic() +
  
  theme(
    axis.title.y = element_text(color = numberColor, size=16),
    axis.title.y.right = element_text(color = areaColor, size=16), 
    axis.text.y.left = element_text(color = numberColor, size = 15),
    axis.text.y.right = element_text(color = areaColor, size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    axis.title.x = element_text(color = "black", size = 16),
    axis.line.y = element_line(linewidth = 1),
    axis.line.x = element_line(linewidth = 1), 
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"))+
  
  scale_x_continuous(limits = c(1958, 2024), breaks = seq(from = 1945, to = 2024, by = 5),guide = "prism_minor", minor_breaks = seq(1945, 2024, 1), expand = c(0,0)) 

PA_expansion_percent_timeline

ggsave("PA_expansion_percent_timeline_fig_280425.png", plot = PA_expansion_percent_timeline, dpi = 300, h = 27, w = 44, unit = "cm")



############################################################################################################################################################################
#create plot for PA expansion over time using area of PAs and number of PAs 
############################################################################################################################################################################
PA_expansion_percent_timeline <- ggplot(cum_area_df, aes(x=Year)) +
  
  geom_rect(aes(xmin = 1958, xmax = 1980, ymin = -10, ymax = -1),
            fill = "#CC79A7") +
  geom_rect(aes(xmin = 1980, xmax = 2000, ymin = -10, ymax = -1),
            fill = "#E69F00")+
  geom_rect(aes(xmin = 2000, xmax = 2024, ymin = -10, ymax = -1),
            fill = "#009E73") +
  geom_area( aes(y=Percent_area*10), size=2, fill = areaColor) +   #the *10 is to put on same scale as other axis so they overlap properly
  geom_line( aes(y=Cum_Num), stat="identity", size=.8, color=numberColor, alpha=1) +
  geom_point( aes(y=Cum_Num), stat="identity", size=1, fill=numberColor, color = numberColor, alpha=1) +
  
  
  geom_segment(aes(x=1962, y=-22, xend=1962, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=1994, y=-13, xend=1994, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=2000, y=-22, xend=2000, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=2010, y=-13, xend=2010, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=2015, y=-22, xend=2015, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  geom_segment(aes(x=2022, y=-13, xend=2022, yend=-29), arrow = arrow(length=unit(0.2, 'cm'))) +
  
  geom_label(label="Partnered with UNESCO", x=1966, y=-22,label.padding = unit(0.3, "lines"), label.size = 0.45,color = "black",fill="lightgrey", size = 7) +
  geom_label(label="Ratified CBD", x=1994, y=-15,label.padding = unit(0.3, "lines"), label.size = 0.45,color = "black",fill="lightgrey", size = 7) +
  geom_label(label="Aligned with MDGs", x=2000, y=-22,label.padding = unit(0.3, "lines"), label.size = 0.45,color = "black",fill="lightgrey", size = 7) +
  geom_label(label="Adopted Aichi targets", x=2005, y=-15,label.padding = unit(0.3, "lines"), label.size = 0.45,color = "black",fill="lightgrey",  size = 7) +
  geom_label(label="Aligned with SDGs", x=2016, y=-22,label.padding = unit(0.3, "lines"), label.size = 0.45,color = "black",fill="lightgrey", size = 7) +
  geom_label(label="Adopted GBF targets", x=2017.5, y=-15,label.padding = unit(0.3, "lines"), label.size = 0.45,color = "black",fill="lightgrey", size = 7) +
  
  scale_y_continuous(expand= c(0,0), limits = c(-29,105),
                     name = "               Cumulative number of PAs",
                     sec.axis = sec_axis(~./10, name="Cumulative PA coverage (%)               ", #here have to re-divide it by 10 to make axis reflect actual values
                                         breaks = seq(0,10,1)),
                     breaks = seq(0,100,10)) +
  
  theme_classic() +
  
  theme(
    axis.title.y = element_text(color = numberColor, size=22),
    axis.title.y.right = element_text(color = areaColor, size=22), 
    axis.text.y.left = element_text(color = numberColor, size = 21),
    axis.text.y.right = element_text(color = areaColor, size = 21),
    axis.text.x = element_text(color = "black", size = 22),
    axis.title.x = element_text(color = "black", size = 21),
    axis.line.y = element_line(linewidth = 1),
    axis.line.x = element_line(linewidth = 1), 
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"))+
  
  scale_x_continuous(limits = c(1958, 2024), breaks = seq(from = 1945, to = 2024, by = 5),guide = "prism_minor", minor_breaks = seq(1945, 2024, 1), expand = c(0,0)) 

PA_expansion_percent_timeline

ggsave("PA_expansion_percent_timeline_fig_280425.png", plot = PA_expansion_percent_timeline, dpi = 300, h = 27, w = 44, unit = "cm")




############################################################################################################################################################################
#add comparison of rate of change from 2010
############################################################################################################################################################################

global_2010 <- 14.1 #from maxwell et al 2020
global_2019 <- 15.3 #from maxwell et al 2020
global_yearly_change <- (global_2019-global_2010)/9 #0.13

Eth_2010 <- 6.7647863
Eth_2019 <- 9.4371098
Eth_yearly_change <- (Eth_2019 - Eth_2010)/10 #0.2672324

diff_change <- Eth_yearly_change/global_yearly_change #2

compare_glob_eth <- data.frame( Year = c(2010,2019),
                                change = c(0, 1.2, 0, 2.67),
                                group = c("Global", "Global", "Ethiopia", "Ethiopia"))

compare_glob_eth$group <- as.factor(compare_glob_eth$group)

library(tidyverse)
library(ggrepel)

rates_compare <- ggplot(compare_glob_eth, aes(x = Year, y = change, group = group, colour = group)) + 
  geom_line(size = 1, show.legend = F) + 
  geom_point(size = 3, show.legend = F)+
  scale_color_manual(values = c("#38A800", "darkgrey")) + 
  theme_classic() +
  scale_x_continuous(limits = c(2009.8,2019.5), expand = c(0,0), breaks = c(2010, 2019)) +
  scale_y_continuous(limits = c(-0.05, 2.8), expand=c(0,0), breaks = seq(0, 2.7, .5)) +
  ylab("Increase in % PA coverage \nsince 2010") + 
  theme(axis.line = element_line(size = 1), axis.title = element_text(size = 13), 
        axis.text = element_text(size = 12), plot.background = element_rect(color = "black", 
                                                                         fill = NA, 
                                                                         size = 0.5)) +
  geom_label(label="Ethiopia", x=2013.75, y=1.7,label.padding = unit(0.3, "lines"), label.size = 0.35, color = "#38A800", fill = "white") +
  geom_label(label="Global", x=2015, y=0.96,label.padding = unit(0.3, "lines"), label.size = 0.35, color = "darkgrey", fill = "white") 
  
 
library(patchwork)
comb_plot <- PA_expansion_percent_timeline + inset_element(rates_compare, left = 0.46, bottom = 0.62, right = 0.71, top =1) 
  #inset_element(legend, left = 0.1, bottom = 0.6, right = 0.4, top = 1, align_to = "full")
comb_plot
ggsave("PA_expansion_percent_timeline_fig_inset2_131124.png", plot = comb_plot, dpi = 300, h = 9, w = 14.7)

##################################################################################################
#simple plot for presentation
################################################################################################

zero_row <- data.frame(
    Year = 1958,
    Cum_Area = 0,
    Cum_Num = 0,
    Percent_area = 0
  )

# Append the row to the data frame
cum_area_df_zero <- rbind(cum_area_df, zero_row)


PA_expansion_percent_simple <- ggplot(cum_area_df, aes(x=Year)) +

  geom_area( aes(y=Percent_area*10), size=2, fill = areaColor) +   #the *10 is to put on same scale as other axis so they overlap properly
  geom_line( aes(y=Cum_Num), stat="identity", size=.8, color=numberColor, alpha=1) +
  geom_point( aes(y=Cum_Num), stat="identity", size=1, fill=numberColor, color = numberColor, alpha=1) +
  
   scale_y_continuous(expand= c(0,0), limits = c(0,105),
                     name = "               Cumulative number of PAs",
                     sec.axis = sec_axis(~./10, name="Cumulative PA coverage (%)               ", #here have to re-divide it by 10 to make axis reflect actual values
                                         breaks = seq(0,10,1)),
                     breaks = seq(0,100,10)) +
  
  theme_classic() +
  
  theme(
    axis.title.y = element_text(color = numberColor, size=16),
    axis.title.y.right = element_text(color = areaColor, size=16), 
    axis.text.y.left = element_text(color = numberColor, size = 15),
    axis.text.y.right = element_text(color = areaColor, size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    axis.title.x = element_text(color = "black", size = 16),
    axis.line.y = element_line(linewidth = 1),
    axis.line.x = element_line(linewidth = 1), 
    axis.ticks = element_line(linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"))+
  
  scale_x_continuous(limits = c(1958, 2024), breaks = seq(from = 1958, to = 2024, by = 10),guide = "prism_minor", minor_breaks = seq(1958, 2024, 2), expand = c(0,0)) 

PA_expansion_percent_simple

ggsave("PA_expansion_percent_simple_fig_291124.png", plot = PA_expansion_percent_simple, dpi = 300, h = 13.5, w = 22, unit = "cm")





############################################################################################################################################################################
#make dummy legend in ggplot for gg_vistime fig
######################################################################################################################

leg_dat <- data.frame(authority = c("Ministry of Agriculture", "Ethiopian Wildlife Conservation Organisation", "Forest and Wildlife Authority", "Ministry of Natural Resources Development and Protection", "Federal Bureau of Wildlife Conservation Authority", "Biodiversity Conservation Research Institute", "Ministry of Agriculture and Rural Development", "Ethiopian Wildlife Conservation Authority"),x=c(1,2,3,4,5,6,7,8),y=c(1,1,1,1,1,1,1,1))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
leg_dat$authority <- factor(leg_dat$authority, levels = c("Ministry of Agriculture", "Ethiopian Wildlife Conservation Organisation", "Forest and Wildlife Authority", "Ministry of Natural Resources Development and Protection", "Federal Bureau of Wildlife Conservation Authority", "Biodiversity Conservation Research Institute", "Ministry of Agriculture and Rural Development", "Ethiopian Wildlife Conservation Authority"))
dummy <- ggplot(data=leg_dat, aes(authority ,y, fill = factor(authority))) + geom_bar(stat="identity") +scale_fill_manual(values=cbPalette) + theme(legend.title = element_blank(), plot.background = element_blank(), p)
dummy
library(cowplot)
leg <- get_legend(dummy)


grid.newpage()   # Create new plot window                           

grid.draw(leg)  # Draw Only legend 

png("legend.png", width     = 10,
    height    = 6,
    units     = "cm",
    res       = 1200)

## Create a graphical object g here
grid.draw(leg) # print it

## Stop writing to the PDF file
dev.off()

legend <- readPNG("legend.png", native = T, info = T)
plot(legend)
############################################################################################################################################################################

#get strict area coverage
IUCN_CAT_II <- ETH_PAs %>%
  filter(IUCN_MG_CA == "II") 

# Print the result
print(IUCN_CAT_II)

# Initialize a dataframe to store the cumulative area and number of PAs for each year
catII_cum_area_df <- data.frame(Year = sort(unique(IUCN_CAT_II$Year)), Cum_Area = 0, Cum_Num = 0)

# Create an empty object to store the cumulative protected areas
catII_cumulative_polygons <- st_sfc(crs = st_crs(IUCN_CAT_II))  # Start with an empty spatial feature

# Initialize a counter for the cumulative number of PAs
cum_num <- 0

# Loop through each year, union the polygons, subtract overlaps, calculate area, and update PA count
for (i in seq_along(catII_cum_area_df$Year)) {
  # Subset polygons established in the current year
  catII_current_year_polygons <- IUCN_CAT_II %>% 
    filter(Year == catII_cum_area_df$Year[i])
  
  print(catII_current_year_polygons$NAME_PAs)
  
  # Update the cumulative number of protected areas
  cum_num <- cum_num + nrow(catII_current_year_polygons)  # Count the number of PAs in the current year
  print(cum_num)
  catII_cum_area_df$Cum_Num[i] <- cum_num  # Store the cumulative number of PAs
  
  # Union the polygons for the current year (this merges any overlapping PAs within the same year)
  catII_unioned_current_year <- st_union(catII_current_year_polygons)
  
  # Subtract any overlapping area with previously protected areas (cumulative polygons)
  catII_new_non_overlapping_area <- st_difference(catII_unioned_current_year, catII_cumulative_polygons)
  
  # Calculate the area of the new non-overlapping area (in km^2)
  catII_new_area <- as.numeric(st_area(catII_new_non_overlapping_area)) / 1e6
  
  # Add the new area to the cumulative polygons
  catII_cumulative_polygons <- st_union(catII_cumulative_polygons, catII_unioned_current_year)
  
  # Update the cumulative area for the current year
  if (i == 1) {
    catII_cum_area_df$Cum_Area[i] <- catII_new_area  # First year, just the new area
  } else {
    catII_cum_area_df$Cum_Area[i] <- catII_cum_area_df$Cum_Area[i-1] + catII_new_area  # Add new area to previous total
  }
  
  
}

# View the result
print(catII_cum_area_df)


catII_cum_area_df$Percent_area <- (catII_cum_area_df$Cum_Area/Eth_area)*100 #calculate percent coverage as new column

