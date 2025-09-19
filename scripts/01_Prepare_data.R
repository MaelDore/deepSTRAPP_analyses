
##### Script 01: Prepare data for analyses #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Curate taxonomic information from species list to prepare for analyses
# Extract species-level inforamtion

###

### Inputs

#

###

### Sources 

# Cerezer, F.O., Dambros, C.S., Coelho, M.T.P. et al. Accelerated body size evolution in upland environments is correlated with recent speciation in South American freshwater fishes.
# Nature Communications. 14, 6070 (2023). https://doi.org/10.1038/s41467-023-41812-7

# Zenodo records: https://zenodo.org/records/8301082

###

### Outputs

# 

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(tidyverse)
library(readxl)
library(xlsx)      # Need the Java Development Kit (JDK) installed
library(openxlsx)  # Use Rccp. No need of Java
library(ape)
library(deepSTRAPP)

### 1.2/ Load Morphological trait data ####

# BEL = Body Elongation
# MBL = Maximum Body Length => The one that is significant
# OGP = Oral Gape Position
# RES = Relative Eye Size
# RML = Relative Maxillary Length

morpho_folder <- "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Biotic_factors/Morphological_evolution/"

BEL_df <- read.table(file = paste0(morpho_folder, "BEL/BAMM_outputs/BEL.txt"))
names(BEL_df) <- c("Taxa", "BEL")
MBL_df <- read.table(file = paste0(morpho_folder, "MBL/BAMM_outputs/MBL.txt"))
names(MBL_df) <- c("Taxa", "MBL")
OGP_df <- read.table(file = paste0(morpho_folder, "OGP/BAMM_outputs/OGP.txt"))
names(OGP_df) <- c("Taxa", "OGP")
RES_df <- read.table(file = paste0(morpho_folder, "RES/BAMM_outputs/RES.txt"))
names(RES_df) <- c("Taxa", "RES")
RML_df <- read.table(file = paste0(morpho_folder, "RML/BAMM_outputs/RMI.txt"))
names(RML_df) <- c("Taxa", "RML")

morpho_tip_data_df <- BEL_df %>% 
  left_join(y = MBL_df) %>%
  left_join(y = OGP_df) %>%
  left_join(y = RES_df) %>%
  left_join(y = RML_df)

### 1.2/ Load the time-calibrated phylogeny

Fish_tree <- read.tree(file = "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Phylogeny/NFF_subset_2638spp.tre")

plot(Fish_tree)

table(Fish_tree$tip.label %in% morpho_tip_data_df$Taxa)

### 1.3/ Load BAMM results for diversification

# Generate an bammdata object from BAMM output
BAMMdata_object_diversification <- BAMMtools::getEventData(
    phy = Fish_tree, 
    eventdata = "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Speciation_rates/BAMM/event_data.txt",
    burnin = 0.2) 

BAMMdata_object_diversification <- BAMMdata_object
head(BAMMdata_object_diversification$eventData)

# Save BAMMdata_object_diversification
saveRDS(object = BAMMdata_object_diversification, file = "./outputs/BAMM/BAMMdata_object_diversification.rds")

### 1.4/ Load BAMM results for trait evolution

## BEL = Body Elongation

# Generate a BAMMdata object from BAMM output
BAMMdata_BEL <- BAMMtools::getEventData(
  phy = Fish_tree, type = "trait",
  eventdata = paste0(morpho_folder, "BEL/BAMM_outputs/event_data.txt"),
  burnin = 0.2) 

head(BAMMdata_BEL$eventData)

# Save BAMMdata_BEL
saveRDS(object = BAMMdata_BEL, file = "./outputs/BAMM/BAMMdata_BEL.rds")

## MBL = Body Elongation

# Generate a BAMMdata object from BAMM output
BAMMdata_MBL <- BAMMtools::getEventData(
  phy = Fish_tree, type = "trait",
  eventdata = paste0(morpho_folder, "MBL/BAMM_outputs/event_data.txt"),
  burnin = 0.2) 

head(BAMMdata_MBL$eventData)

# Save BAMMdata_MBL
saveRDS(object = BAMMdata_MBL, file = "./outputs/BAMM/BAMMdata_MBL.rds")


##### 2/ Convert BAMM trait data into contMap ####

# Can I convert this to a contMap with mean trait values along branches ???


##### 3/ Categorize trait data #####

# It trait data are categorized, need to rerun the 

##### 4/ Extract abiotic variables at species level #####

library(sf)

Basin_shp <- read_sf(dsn = "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Shapefiles/Hydroatlas/", layer = "BasinATLAS_v10_lev05")

plot(Basin_shp[, "MAIN_BAS"])

### 4.1/ Get binary table of subbasin occurrence ####

Basin_occurence_binary_table_df <- read.csv("./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Occurrences/PresAbs_SubBasin_subset_2638spp.csv")

### 4.2/ Get subbasin level data ####

SubBasins_data_df <- read.csv("./input_data/Cerezer_2023_Neotropical_freshwater_fishes/All_data_combined/SubBasins_estimates.csv")

# bio1 = Mean temperature
# bio12 = Annual total Precipitation
# aet = Evapotranspiration (mm / syr)
# Elevation

### 4.3/ Aggregate at species level

# Keep only shared basins
SubBasins_data_df_subset <- SubBasins_data_df %>%
  filter(HYBAS_ID %in% Basin_occurence_binary_table_df$HYBAS_ID) %>% 
  select(HYBAS_ID, bio1, bio12, aet, Elevation)
Basin_occurence_binary_table_df_subset <- Basin_occurence_binary_table_df %>% 
  filter(HYBAS_ID %in% SubBasins_data_df$HYBAS_ID)

# Ensure basins are similarly ordered
Basin_occurence_binary_table_df_subset$HYBAS_ID == SubBasins_data_df_subset$HYBAS_ID

# Remove basin ID for computation
row.names(Basin_occurence_binary_table_df_subset) <- Basin_occurence_binary_table_df_subset$HYBAS_ID
Basin_occurence_binary_table_df_subset <- Basin_occurence_binary_table_df_subset %>% 
  select(-HYBAS_ID)
row.names(SubBasins_data_df_subset) <- SubBasins_data_df_subset$HYBAS_ID
SubBasins_data_df_subset<- SubBasins_data_df_subset %>% 
  select(-HYBAS_ID)

# Use matrix product to compute the sum of variable values where each species is found
abiotic_data_total <- t(as.matrix(Basin_occurence_binary_table_df_subset)) %*% as.matrix(SubBasins_data_df_subset)
# Divide by the number of basins to get average per species
basin_prevalence <- colSums(Basin_occurence_binary_table_df_subset)
abiotic_data_mean <- abiotic_data_total / basin_prevalence
  
head(abiotic_data_mean)

# Convert to df

abiotic_data_df <- as.data.frame(abiotic_data_mean)
abiotic_data_df$Taxa <- row.names(abiotic_data_df)

abiotic_data_df <- abiotic_data_df %>% 
  rename(Temp = bio1,
         Prec = bio12,
         Evapotranspiration = aet) %>%
  select(Taxa, Elevation, Temp, Prec, Evapotranspiration)

head(abiotic_data_df)

# Save abiotic data per species
saveRDS(object = abiotic_data_df, file = "./outputs/abiotic_data_df.rds")


##### 5/ Extract regions based on subbasin occurrences #####

# Use definition of regions in the Neotropics as in Morrone, 2022? (Terrestrial and floristic...)
# https://doi.org/10.1590/0001-3765202220211167 
# Morrone, 2022
#  - Antillean, Brazilian and Chacoan subregions
#  - Mexican and South American transition zones

library(sf)

Morrone_shp <- read_sf(dsn = "./input_data/Morrone_2022/NeotropicMap_Geo/", layer = "NeotropicMap_Geo")

plot(Morrone_shp[, "Subregion"])

## Make a union between the two shp files and attribute to each basin the Subregion with the most overlap in terms of area

## Then assign Subregion to species based on Sub-basins


# May still be too many regions for proper deepSTRAPP testing. Need to aggregate? Check patterns in the initial paper results

# Results from Boris's paper at lower level than Region? Does not seem to be available. Check his paper resutls to see if he explains if the clustering stopped at level 2

Leroy_output <- readRDS(file = "./input_data/Leroy_2019/metrics_per_site_wide.RDS")

#### 6/ Extract additional ecological traits ####

# See https://www.nature.com/articles/s41597-025-04674-w


