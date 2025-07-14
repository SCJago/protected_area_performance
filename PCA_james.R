library("virtualspecies")
library("raster")
library("ENMTools")
library("ENMeval")
library("rgeos")
library("rgdal")
library("usdm")
library("geosphere")
library("scales")
library("factoextra")
library("vegan")
library("rgeos")
library("sp")
library("png")
library("ggplot2")
library("alphahull")
library("rgeos")
library("scales")
library('rgdal')
library('factoextra')
library('alphahull')
library('dplyr')

#Setting working directory 
setwd("~/Papers/Ethiopia - PA history review/Representativeness")

#load ethiopia boundary
ETH<-getData("GADM", country= "ETH",level=0)


#extract environmental rasters for ethiopia
all_raster<-stack("chelsa_bioclim_stack.grd") 
names(all_raster) <- c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")
all_raster_ETH <- crop(all_raster, ETH)
all_raster_ETH <-  mask(all_raster, ETH)

#Then take  your projected areas and import as a shape file 
pa_shape = readOGR("~/Papers/Ethiopia - PA history review/EWCA_updated_PAs_May23/ETH-WPAs2023 (2)/ETH-WPAs2023/ETH-WPAs-2023-ed.shp")
pa_shape<- spTransform(pa_shape, CRS(proj4string(all_raster)))

#Crop and mask the Ethiopian climate data to just this PA 
pa_data <- crop(all_raster,pa_shape)
pa_data <- mask(pa_data,pa_shape) 

#Randomly sample some data from for background env and PA space 
data_for_pca2 <- sampleRandom(all_raster_ETH, 100000, na.rm=TRUE, xy=T)
pa_data_extracted <- sampleRandom(pa_data, 9400, na.rm=TRUE, xy=T) 

#run pca
pca2 <- prcomp(data_for_pca2[,3:21], center = T, scale = T) 

#put PAs in pca space
PAs <- predict(pca2, newdata = pa_data_extracted) 

#plot pca with pa points
png("~/Papers/Ethiopia - PA history review/Representativeness/Allraster_PAs_PCA_300425.png", units = "cm", 
    width = 15, height = 15, res = 300)
plot(pca2$x[,1],pca2$x[,2], pch=16, cex=1, col = alpha("grey", 0.2), 
     xlab="PC1 (40.1%)", ylab="PC2 (22.5%)") 
points (PAs, pch=16, cex=1, col = alpha("darkgreen", 0.2) )
dev.off()


#Contribution of variables
fviz_pca_ind(pca2, geom.ind="point",
             addEllipses = T, ellipse.level = 0.95)

#plot contribution of variables to pc1
png("~/Papers/Ethiopia - PA history review/Representativeness/Contribution_PC1_300425png", units = "cm", 
    width = 15, height = 15, res = 300)
fviz_contrib(pca2, choice = "var", axes = 1)
dev.off()


#plot contribution of variables to pc1
png("~/Papers/Ethiopia - PA history review/Representativeness/Contribution_PC2_300425.png", units = "cm", 
    width = 15, height = 15, res = 300)
fviz_contrib(pca2, choice = "var", axes = 2)
dev.off()

###############################################################################

#get alpha hull
unique_PAs <- unique(PAs[, 1:2])
unique_pca2 <- unique(pca2$x[,1:2])
hull <- ahull(x = unique_PAs[,1], y = unique_PAs[,2], alpha = 0.5)
full_hull <- ahull(unique_pca2, alpha = 0.5)

#plot with alpha hull
png("~/Papers/Ethiopia - PA history review/Representativeness/alphahulls.png", units = "cm", 
    width = 15, height = 15, res = 300)
plot(pca2$x[,1],pca2$x[,2], pch=16, cex=1, col = alpha("grey", 0.2), 
     xlab="PC1 (40.1%)", ylab="PC2 (22.5%)") 
points (PAs, pch=16, cex=1, col = alpha("darkgreen", 0.2)) 
plot(full_hull, add=T, col = "black", wpoints=F )
plot(hull, add = T, col = "darkgreen", wpoints = F)
dev.off()

area_PA_space <- areaahull(hull, timeout = 5)
area_env_space <- areaahull(full_hull, timeout = 5)

percent_overlap <- (area_PA_space/area_env_space)*100 #54% with 0.5 alpha
percent_overlap

# now lets find the points not inside the hull
vector <- inahull(hull, pca2$x[,1:2])

# and then bind that to the original dataframe to get the coordinates
missing_regions <- cbind.data.frame(data_for_pca2[,1:2],vector)


# then we take a map of Ethiopia and we can colour in places that were covered and places that were not

png("~/Papers/Ethiopia - PA history review/Representativeness/Representativeness_map.png", units = "cm", 
    width = 15, height = 15, res = 300)
plot(ETH)
missing_regions_ordered <- missing_regions %>%
  arrange(vector)
colours <- c('grey', 'darkgreen')[factor(missing_regions_ordered$vector, levels = c('FALSE', 'TRUE'))]
sizes <- c(2,.5)[factor(missing_regions_ordered$vector, levels = c("FALSE", "TRUE"))]
points(missing_regions_ordered[,1:2], col=colours, pch=20, cex=sizes)
dev.off()



ggplot(missing_regions_ordered, aes(x=x, y=y, col = vector, size = vector))+geom_point()+
  scale_colour_manual(values=c("grey", "darkgreen"))+
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 2))+
  theme_classic()+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())



