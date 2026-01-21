
##### Script 04: Map biodiversity patterns #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Map current species richness across sub-basins
# Map mean current net diversification rates from BAMM across sub-basins

###

### Inputs

# Sub-basins sf object
# BAMM output for diversification dynamics
# Morpho data summarized per sub-basins

###

### Sources 

## BAMM analyses and sub-basins data

# Cerezer, F.O., Dambros, C.S., Coelho, M.T.P. et al. Accelerated body size evolution in upland environments is correlated with recent speciation in South American freshwater fishes.
# Nature Communications. 14, 6070 (2023). https://doi.org/10.1038/s41467-023-41812-7

### Outputs

# Map current species richness across sub-basins
# Map mean current net diversification rates from BAMM across sub-basins
# Map mean MBL from BAMM across sub-basins

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(sf)

### 1.2/ Load data per sub-basins ####

# Load 
SubBasins_estimates <- read_csv("./input_data/Cerezer_2023_Neotropical_freshwater_fishes/All_data_combined/SubBasins_estimates.csv")
SubBasins_estimates$HYBAS_ID <- as.character(SubBasins_estimates$HYBAS_ID)

# Check correlation between speciaiton rates (as used in Cerezer et al., 2023) and net diversification rates (as used in this study)
plot(SubBasins_estimates$BAMM_speciation, SubBasins_estimates$BAMM_NetDiv)
cor.test(x = SubBasins_estimates$BAMM_speciation, y = SubBasins_estimates$BAMM_NetDiv, method = "spearman")
# Spearman's rho = 0.945

### 1.3/ Load shapefiles for sub-basins ####

Basin_shp <- read_sf(dsn = "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Shapefiles/Hydroatlas/", layer = "BasinATLAS_v10_lev05")
Basin_shp <- st_make_valid(Basin_shp)
st_is_valid(Basin_shp)

plot(Basin_shp[, "MAIN_BAS"])
plot(Basin_shp[, "HYBAS_ID"])

### 1.4/ Load map stuff ####

# Load Morrone subregions
Morrone_Subregions_shp <- readRDS(file = "./input_data/Morrone_2022/Morrone_Subregions_shp.rds")

# Load color palette
pal_bl_red_Mannion <- readRDS(file = "./input_data/Maps/pal_bl_red_Mannion.rds")


##### 2/ Retrieve data per sub-basins ####

### 2.1/ Retrieve BAMM rates ####

Basin_data_sf <- Basin_shp %>% 
  left_join(SubBasins_estimates) %>%
  select("HYBAS_ID", "MAIN_BAS", "BAMM_NetDiv", "MBL_evol", "diversity")

### 2.2/ Retrieve mean MBL data (not rates!) ####

# Load trait data
morpho_tip_data_df <- readRDS(file = "./outputs/Morpho/morpho_tip_data_df.rds")

# Load binary table of sub-basin occurrence ####
Basin_occurrence_binary_table_df <- read.csv("./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Occurrences/PresAbs_SubBasin_subset_2638spp.csv")

## Filter data to intersect species lists

row.names(Basin_occurrence_binary_table_df) <- Basin_occurrence_binary_table_df$HYBAS_ID

Basin_occurrence_binary_table_df <- Basin_occurrence_binary_table_df %>% 
  select(-HYBAS_ID)

# Reorder binary table as in morpho_tip_data
Basin_occurrence_binary_table_df <- Basin_occurrence_binary_table_df[, morpho_tip_data_df$Taxa]

table(morpho_tip_data_df$Taxa == names(Basin_occurrence_binary_table_df))

## Compute mean MBL per sub-basins

dim(as.matrix(morpho_tip_data_df$MBL))
dim(Basin_occurrence_binary_table_df)

# Compute sum of MBL
MBL_per_subbasins <- t(as.matrix(morpho_tip_data_df$MBL)) %*% t(as.matrix(Basin_occurrence_binary_table_df))

# Divide by species richness to get mean
mean_MBL_per_subbasins <- as.vector(MBL_per_subbasins / rowSums(Basin_occurrence_binary_table_df))
names(mean_MBL_per_subbasins) <- row.names(Basin_occurrence_binary_table_df)

# Replace NaN due to division by zero with NA
mean_MBL_per_subbasins[is.nan(mean_MBL_per_subbasins)] <- NA

hist(mean_MBL_per_subbasins)

## Store in Sub-basin data object
mean_MBL_per_subbasins_df <- data.frame(mean_MBL = mean_MBL_per_subbasins, HYBAS_ID = names(mean_MBL_per_subbasins))

Basin_data_sf <- Basin_data_sf %>% 
  left_join(mean_MBL_per_subbasins_df)

## Save subbasin sf with associated data
saveRDS(object = Basin_data_sf, file = "./outputs/Morpho/Basin_data_sf.rds")


##### 3/ Plot maps with ggplot #####

### 3.1/ Build contour maps ####

## Load subbasin sf with associated data
Basin_data_sf <- readRDS(file = "./outputs/Morpho/Basin_data_sf.rds")

## Build contour map of SA
South_America_contour_sf <- st_union(Basin_data_sf)

plot(South_America_contour_sf)

# Save contour map of SA
saveRDS(object = South_America_contour_sf, file = "./outputs/South_America_contour_sf.rds")

## Clean and simplify Morrone subregions

Morrone_Subregions_clean_sf <- sf::st_crop(x = Morrone_Subregions_shp, y = South_America_contour_sf) %>% 
  st_simplify(dTolerance = 10000) %>% # Simply shape with tolerance in map units (e.g., meters)
  st_intersection(y = South_America_contour_sf) # Keep within overall contours

plot(Morrone_Subregions_clean_sf)

# Save Morrone subregions
saveRDS(object = Morrone_Subregions_clean_sf, file = "./outputs/Morrone_Subregions_clean_sf.rds")


### 3.2/ Plot Morrone's subregions ####

colors_per_ranges <- c("limegreen", "gold", "dodgerblue3")
names(colors_per_ranges) <- c("Brazilian", "Chacoan", "South American Transition Zone")

## Plot Morrone's subregions in WGS84
Morrone_subregions_ggplot_WGS84 <- ggplot(data = South_America_contour_sf) +
  
  # Plot contours
  geom_sf(data = South_America_contour_sf, 
          colour = "black",
          fill = "#EDEDED",
          alpha = 1.0) +
  
  # Plot Morrone subregions
  geom_sf(data = Morrone_Subregions_clean_sf,
          aes(fill = Subregion),
          colour = "black",
          linewidth = 1.5,
          alpha = 1.0) +
  
  # Plot sub-basins
  geom_sf(data = Basin_data_sf, 
          colour = "black",
          fill = NA,
          alpha = 1.0) +
  
  # Plot contours
  geom_sf(data = South_America_contour_sf, 
          colour = "black",
          fill = NA,
          alpha = 1.0) +
  
  # Adjust color scheme and legend
  scale_fill_manual("Subregions", breaks = names(colors_per_ranges),
                    labels = c("Brazilian", "Chacoan", "SATZ"), values = colors_per_ranges) +
  
  # # Adjust CRS
  # # coord_sf(default_crs = sf::st_crs(4326)) +
  # coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
  #          clip = "off", # To allow plotting arrow outside of World map
  #          expand = FALSE) +
  
  # Add title
  ggtitle(label =  paste0("Morrone's Subregions")) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 10)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 12, face = "bold"),
        legend.box.margin = margin(l = 5),
        legend.position = c(0.9, 0.20),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())

pdf(file = paste0("./outputs/Biogeo/Morrone_Subregions_ggplot_WGS84.pdf"),
    width = 5, height = 8)

print(Morrone_subregions_ggplot_WGS84)

dev.off()



### 3.2/ Plot mean BAMM net diversification rates ####

# Make BAMM color gradient
BAMM_color_palette <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))

# Inspect CRS
sf::st_crs(Basin_data_sf) # WGS84

hist(Basin_data_sf$BAMM_NetDiv)
hist(log(Basin_data_sf$BAMM_NetDiv)) # Use log-transformation

## Plot BAMM net diversification rates in WGS84
BAMM_NetDiv_ggplot_WGS84 <- ggplot(data = South_America_contour_sf) +
  
  # Plot contours
  geom_sf(data = South_America_contour_sf, 
          colour = "black",
          fill = "#EDEDED",
          alpha = 1.0) +
  
  # Plot bioregion sf maps
  geom_sf(data = Basin_data_sf, 
          mapping = aes(fill = BAMM_NetDiv),
          # colour = "black",
          colour = NA,
          alpha = 1.0) +
  
  # Plot contours
  geom_sf(data = South_America_contour_sf, 
          colour = "black",
          fill = NA,
          alpha = 1.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Net div. rates", colors = BAMM_color_palette, na.value = "#EDEDED",
                       transform = "log", labels = scales::label_number(accuracy = 0.01)) +
  
  # # Adjust CRS
  # # coord_sf(default_crs = sf::st_crs(4326)) +
  # coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
  #          clip = "off", # To allow plotting arrow outside of World map
  #          expand = FALSE) +
  
  # Add title
  ggtitle(label =  paste0("BAMM net diversification rates")) +
  
  # Adjust legend aesthetics
  guides(fill = guide_colorbar(barwidth = 2, barheight = 10)) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 10)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 12, face = "bold"),
        legend.box.margin = margin(l = 5),
        legend.position = c(0.9, 0.25),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())

pdf(file = paste0("./outputs/Biogeo/BAMM_NetDiv_ggplot_WGS84.pdf"),
    width = 5, height = 8)

print(BAMM_NetDiv_ggplot_WGS84)

dev.off()


### 3.3/ Plot Species Richness ####

# Inspect CRS
sf::st_crs(Basin_data_sf) # WGS84

hist(Basin_data_sf$diversity)

## Plot Species richness in WGS84
SR_ggplot_WGS84 <- ggplot(data = South_America_contour_sf) +
  
  # Plot contours
  geom_sf(data = South_America_contour_sf, 
          colour = "black",
          fill = "#EDEDED",
          alpha = 1.0) +
  
  # Plot bioregion sf maps
  geom_sf(data = Basin_data_sf, 
          mapping = aes(fill = diversity),
          # colour = "black",
          colour = NA,
          alpha = 1.0) +
  
  # # Plot Morrone subregions
  # geom_sf(data = Morrone_Subregions_clean_sf, 
  #         colour = "white",
  #         fill = NA,
  #         linewidth = 1.5,
  #         alpha = 1.0) +
  
  # Plot contours
  geom_sf(data = South_America_contour_sf, 
          colour = "black",
          fill = NA,
          alpha = 1.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("Species\nrichness", colors = pal_bl_red_Mannion, na.value = "#EDEDED") +
  
  # # Adjust CRS
  # # coord_sf(default_crs = sf::st_crs(4326)) +
  # coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
  #          clip = "off", # To allow plotting arrow outside of World map
  #          expand = FALSE) +
  
  # Add title
  ggtitle(label =  paste0("Species Richness")) +
  
  # Adjust legend aesthetics
  guides(fill = guide_colorbar(barwidth = 2, barheight = 10)) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 10)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 12, face = "bold"),
        legend.box.margin = margin(l = 5),
        legend.position = c(0.9, 0.25),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())

pdf(file = paste0("./outputs/Biogeo/SR_ggplot_WGS84.pdf"),
    width = 5, height = 8)

print(SR_ggplot_WGS84)

dev.off()


### 3.4/ Plot mean MBL ####

# Make Spectral color gradient for MBL
Spectral_color_palette <- rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral"))

# Inspect CRS
sf::st_crs(Basin_data_sf) # WGS84

hist(Basin_data_sf$mean_MBL)

## Plot mean MBL in WGS84
mean_MBL_ggplot_WGS84 <- ggplot(data = South_America_contour_sf) +
  
  # Plot contours
  geom_sf(data = South_America_contour_sf, 
          colour = "black",
          fill = "#EDEDED",
          alpha = 1.0) +
  
  # Plot bioregion sf maps
  geom_sf(data = Basin_data_sf, 
          mapping = aes(fill = mean_MBL),
          # colour = "black",
          colour = NA,
          alpha = 1.0) +
  
  # Plot contours
  geom_sf(data = South_America_contour_sf, 
          colour = "black",
          fill = NA,
          alpha = 1.0) +
  
  # Adjust color scheme and legend
  scale_fill_gradientn("MBL", colors = Spectral_color_palette, na.value = "#EDEDED") +
  
  # # Adjust CRS
  # # coord_sf(default_crs = sf::st_crs(4326)) +
  # coord_sf(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
  #          clip = "off", # To allow plotting arrow outside of World map
  #          expand = FALSE) +
  
  # Add title
  ggtitle(label =  paste0("Mean Maximum Body Length")) +
  
  # Adjust legend aesthetics
  guides(fill = guide_colorbar(barwidth = 2, barheight = 10)) +
  
  # Adjust aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill = NA, colour = NA),
        # panel.grid.major = element_line(colour = "grey70", linetype = "dashed", linewidth = 0.5), # Plot graticules
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(b = 10)),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), # trbl
        legend.title = element_text(size = 14, face = "bold", margin = margin(b = 15)),
        legend.text = element_text(size = 12, face = "bold"),
        legend.box.margin = margin(l = 5),
        legend.position = c(0.9, 0.25),
        # axis.ticks = element_line(linewidth = 1.0),
        # axis.ticks.length = unit(10, "pt"),
        # axis.text = element_text(size = 21, color = "black", face = "bold"),
        # axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        # axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank())

pdf(file = paste0("./outputs/Morpho/mean_MBL_ggplot_WGS84.pdf"),
    width = 5, height = 8)

print(mean_MBL_ggplot_WGS84)

dev.off()
