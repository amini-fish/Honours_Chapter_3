install.packages(c("cowplot", 
                   "googleway", 
                   "ggplot2", 
                   "ggrepel", 
                   "ggspatial", 
                   "libwgeom", 
                   "sf", 
                   "rnaturalearth", 
                   "rnaturalearthdata", 
                   "ozmaps"))


install.packages("ggplot2", version='3.4.9')

install.packages("leaflet")
install.packages("rgdal")
install.packages("raster")

library(ozmaps)   
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(leaflet) #for making the interactive map
library(rgdal) #for importing vector data
library(raster)
library(RColorBrewer)

## alrighty now were goot to go

theme_set(theme_bw()) 

world <- ne_countries(scale = "medium", returnclass = "sf")

class(world)

## Load in our spatial data ##

setwd("C:/Users/samue/Desktop/Honours/Maps")

sf_oz <- ozmap("states")

NT <- st_read("D:/Research/Neotrygon/Spatial/NT_STATE_POLYGON_shp.shp")

Daly <- st_read("C:/Users/samue/Desktop/Honours/Maps/WaterResources_daly_500/Daly_Catchment_500/data/daly_500_bnd.shp")

crs(Daly)

points <- read.csv("C:/Users/samue/Desktop/Honours/analysis/Daly_meta.csv")

range <-st_read("C:/Users/samue/Desktop/Honours/Maps/pristis_range/data_0.shp")

## Dist map


PPrange <- ggplot(data = ozmap_country) +
  geom_sf(fill = NA)+
  geom_sf(data = range, fill = "#ADEBB3") +
  coord_sf(xlim = c(110, 158), ylim = c(-5, -45), expand = FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

print(PPrange)

## 
## Now we can do a rough plot

Daly_map <- ggplot(data = NT) +
  geom_sf(fill = NA)+
  geom_sf(data = Daly, fill = "lightgrey") +
  geom_point(data = points, 
             mapping = aes(x = Long, y = Lat), col = "Green") +
  coord_sf(xlim = c(129, 138), ylim = c(-10.5, -18), expand = FALSE) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

  print(Daly_map)

## Try zooming in on the study area 
Daly_map <- ggplot(data = NT) +
    geom_sf(fill = NA)+
    geom_sf(data = Daly, fill = "lightgrey") +
    geom_point(data = points, 
               mapping = aes(x = Long, y = Lat), col = "Darkgreen", size = 4) +
    coord_sf(xlim = c(130, 132), ylim = c(-13, -15), expand = FALSE) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  print(Daly_map)
  