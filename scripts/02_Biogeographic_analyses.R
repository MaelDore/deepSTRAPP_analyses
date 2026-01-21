
##### Script 02: Biogeographic analyses #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Run deepSTRAPP on biogeographic range data

###

### Inputs

# Time-calibrated phylogeny
# Binary table of species occurrence in subregions

###

### Outputs

# densityMaps for ancestral ranges
# Plot p-values STRAPP tests throuh time
# Plot rates through time

###


# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### 1.1/ Load packages ####

library(deepSTRAPP)
library(ape)
library(phytools)

### 1.2/ Load the time-calibrated phylogeny

Fish_tree <- read.tree(file = "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Phylogeny/NFF_subset_2638spp.tre")
length(Fish_tree$tip.label)

nodeHeights(Fish_tree)

Fish_tree <- force.ultrametric(Fish_tree)

### 1.3/ Load Biogeographic range data ####

Subregions_occurrence_binary_table <- readRDS(file = "./outputs/Biogeo/Subregions_occurrence_binary_table.rds")

table(Fish_tree$tip.label == row.names(Subregions_occurrence_binary_table))

# Reorder the table
Biogeo_tip_data_table <- as.data.frame(Subregions_occurrence_binary_table[Fish_tree$tip.label, ])

# Convert it to BioGeoBEARS format
names(Biogeo_tip_data_table) <- c("B", "C", "S")

### 1.4/ Load BAMM results for diversification

# Load BAMMdata_object_diversification
BAMMdata_object_diversification <- readRDS(file = "./outputs/BAMM/BAMMdata_object_diversification.rds")
length(BAMMdata_object_diversification$eventData)


##### 2/ Prepare data = run BioGEOBEARS and get densityMaps #####

colors_per_ranges <- c("limegreen", "gold", "dodgerblue3")
names(colors_per_ranges) <- c("B", "C", "S")

## Run evolutionary models
Fish_biogeo_data <- prepare_trait_data(
  tip_data = Biogeo_tip_data_table,
  trait_data_type = "biogeographic",
  phylo = Fish_tree,
  # Default = "DEC" for biogeographic
  evolutionary_models = c("BAYAREALIKE", "DIVALIKE", "DEC",
                          "BAYAREALIKE+J", "DIVALIKE+J", "DEC+J"),
  prefix_for_files = "Fish",
  max_range_size = 3,
  split_multi_area_ranges = TRUE, # Set to TRUE to display the two outputs
  # Reduce the number of Stochastic Mapping simulations to save time (Default = '1000')
  nb_simulations = 100,
  colors_per_levels = colors_per_ranges,
  return_simmaps = FALSE,
  return_best_model_fit = TRUE,
  return_model_selection_df = TRUE,
  verbose = TRUE) 

# Save output
saveRDS(object = Fish_biogeo_data, file = "./outputs/Biogeo/Fish_biogeo_data.rds")

## Explore output
str(Fish_biogeo_data, 1)

# Summary of model selection
Fish_biogeo_data$model_selection_df 
# Parameter estimates and optimization summary of the best model
# (Here, the best model is DEC+J)
Fish_biogeo_data$best_model_fit$optim_result

# Posterior probabilities of each state (= ACE) at internal nodes
Fish_biogeo_data$ace # Only with unique areas
Fish_biogeo_data$ace_all_ranges # Including multi-area ranges (Here, BC, BS, CS, BCS)

## Plot densityMaps
# densityMap for range n°1 ("B")
plot(Fish_biogeo_data$densityMaps[[1]])
# densityMaps with all unique areas overlaid
plot_densityMaps_overlay(Fish_biogeo_data$densityMaps)
# densityMaps with all ranges (including multi-area ranges) overlaid
plot_densityMaps_overlay(Fish_biogeo_data$densityMaps_all_ranges)

## Export summary of model selection
write.csv(x = Fish_biogeo_data$model_selection_df, file = "./outputs/Biogeo/model_selection_df.csv")

##### 3/ Run deepSTRAPP analyses #####

### Load modeled biogeographic data
Fish_biogeo_data <- readRDS(file = "./outputs/Biogeo/Fish_biogeo_data.rds")

## Convert biogeo table to vector of character strings
Biogeo_tip_data_table_df <- as.data.frame(Biogeo_tip_data_table)
Biogeo_tip_data <- setNames(object = rep(NA, times = nrow(Biogeo_tip_data_table_df)), nm = row.names(Biogeo_tip_data_table_df))
for (i in 1:nrow(Biogeo_tip_data_table_df))
{
  # i <- 1
  
  # Extract unique areas
  unique_areas_i <- names(Biogeo_tip_data_table_df)[Biogeo_tip_data_table_df[i, ] == 1]
  unique_areas_i <- unique_areas_i[order(unique_areas_i)]
  
  # Collapse into range
  range_i <- paste(unique_areas_i, collapse = "")
  # Store range
  Biogeo_tip_data[i] <- range_i
}
Biogeo_tip_data

## Set for five time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 100 Mya.
time_step_duration <- 5
time_range <- c(0, 100)

# Run deepSTRAPP on net diversification rates
## This step is time-consuming. You can skip it and load directly the result if needed
Fish_deepSTRAPP_biogeo_0_100 <- run_deepSTRAPP_over_time(
  densityMaps = Fish_biogeo_data$densityMaps,
  ace = Fish_biogeo_data$ace,
  tip_data = Biogeo_tip_data,
  trait_data_type = "biogeographic",
  BAMM_object = BAMMdata_object_diversification,
  time_range = time_range,
  time_step_duration = time_step_duration,
  seed = 1234, # Set seed for reproducibility
  alpha = 0.05,
  # Run post hoc tests too
  posthoc_pairwise_tests = TRUE,
  # Needed to obtain STRAPP stats and plot evaluation histograms (See 4.2)
  return_perm_data = TRUE, 
  # Needed to get trait data and plot rates through time (See 4.3)
  extract_trait_data_melted_df = TRUE,
  # Needed to get diversification data and plot rates through time (See 4.3)
  extract_diversification_data_melted_df = TRUE, 
  # Needed to obtain STRAPP stats and plot evaluation histograms (See 4.2)
  return_STRAPP_results = TRUE, 
  # Needed to plot updated densityMaps (See 4.4)
  return_updated_trait_data_with_Map = TRUE, 
  # Needed to map diversification rates on updated phylogenies (See 4.5)
  return_updated_BAMM_object = TRUE, 
  verbose = TRUE,
  verbose_extended = TRUE)


## Explore output
str(Fish_deepSTRAPP_biogeo_0_100, max.level = 1)

# Display test summary
# Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
# showing the evolution of the test results across time
Fish_deepSTRAPP_biogeo_0_100$pvalues_summary_df

# Access STRAPP test results
# Can be passed down to [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] to generate plot
# showing the null distribution of the test statistics
str(Fish_deepSTRAPP_biogeo_0_100$STRAPP_results, max.level = 2)

# Save results
# saveRDS(Fish_deepSTRAPP_biogeo_0_50, file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0_50.rds")
# saveRDS(Fish_deepSTRAPP_biogeo_0_50, file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0_50_backup.rds")
saveRDS(Fish_deepSTRAPP_biogeo_0_100, file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0_100.rds")


##### 4/ Plot results #####

# Load deepSTRAPP results
Fish_deepSTRAPP_biogeo_0_100 <- readRDS(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0_100.rds")

### 4.1/ Plot evolution of STRAPP tests p-values through time ####

## Overall Kruskal-Wallis tests

pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0_100_pvalues_over_time_ggplot.pdf", width = 8, height = 6)
deepSTRAPP::plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  # time_range = c(0, 40),
  time_range = c(0, 80),
  alpha = 0.05)
dev.off()

## Post-hoc Dunn's tests

pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0_100_posthoc_pvalues_over_time_ggplot.pdf", width = 8, height = 6)
deepSTRAPP::plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  # time_range = c(0, 40),
  time_range = c(0, 80),
  plot_posthoc_tests = TRUE) # To plot results of post hoc pairwise tests instead
dev.off()

### 4.2/ Plot histogram of STRAPP test stats ####

## Overall Kruskal-Wallis tests

pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0My_overall_stats_histo_ggplot.pdf", width = 8, height = 6)
# Plot the histogram of overall Kruskal-Wallis stats for time-step n°1 = 0 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 0)
dev.off()

# Plot the histogram of overall Kruskal-Wallis stats for time-step n°5 = 20 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_20My_overall_stats_histo_ggplot.pdf", width = 8, height = 6)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 20)
dev.off()

# Plot the histogram of overall Kruskal-Wallis stats for time-step n°9 = 40 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_40My_overall_stats_histo_ggplot.pdf", width = 8, height = 6)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 40)
dev.off()

# Plot the histogram of overall Kruskal-Wallis stats for time-step n°13 = 60 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_60My_overall_stats_histo_ggplot.pdf", width = 8, height = 6)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 60)
dev.off()

## Post-hoc Dunn's tests

# Plot the histogram of overall Kruskal-Wallis stats for time-step n°1 = 0 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0My_posthoc_stats_histo_ggplot.pdf", width = 10, height = 7.5)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  plot_posthoc_tests = TRUE, # To plot results of post hoc pairwise tests instead
  focal_time = 0)
dev.off()

# Plot the histogram of overall Kruskal-Wallis stats for time-step n°5 = 20 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_20My_posthoc_stats_histo_ggplot.pdf", width = 10, height = 7.5)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  plot_posthoc_tests = TRUE, # To plot results of post hoc pairwise tests instead
  focal_time = 20)
dev.off()

# Plot the histogram of overall Kruskal-Wallis stats for time-step n°9 = 40 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_40My_posthoc_stats_histo_ggplot.pdf", width = 10, height = 7.5)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  plot_posthoc_tests = TRUE, # To plot results of post hoc pairwise tests instead
  focal_time = 40)
dev.off()

# Plot the histogram of overall Kruskal-Wallis stats for time-step n°13 = 60 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_60My_posthoc_stats_histo_ggplot.pdf", width = 10, height = 7.5)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  plot_posthoc_tests = TRUE, # To plot results of post hoc pairwise tests instead
  focal_time = 60)
dev.off()

### 4.3/ Plot evolution of rates though time in relation to trait values ####

colors_per_ranges <- c("limegreen", "gold", "dodgerblue3")
names(colors_per_ranges) <- c("B", "C", "S")

pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0_100_RTT_ggplot.pdf", width = 10, height = 6)
plot_rates_through_time(deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100, 
                        colors_per_levels = colors_per_ranges,
                        # time_range = c(0, 40),
                        time_range = c(0, 80),
                        plot_CI = TRUE)
dev.off()

### 4.4/ Plot rates vs. states across branches for a given 'focal_time' ####

# Generate ggplot for time = 0 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_0My_rates_vs_ranges_ggplot.pdf", width = 8, height = 6)
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 0,
  colors_per_levels = colors_per_ranges)
dev.off()

# Generate ggplot for time = 20 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_20My_rates_vs_ranges_ggplot.pdf", width = 8, height = 6)
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 20,
  colors_per_levels = colors_per_ranges)
dev.off()

# Generate ggplot for time = 40 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_40My_rates_vs_ranges_ggplot.pdf", width = 8, height = 6)
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 40,
  colors_per_levels = colors_per_ranges)
dev.off()

# Generate ggplot for time = 60 My
pdf(file = "./outputs/Biogeo/Fish_deepSTRAPP_biogeo_60My_rates_vs_ranges_ggplot.pdf", width = 8, height = 6)
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 60,
  colors_per_levels = colors_per_ranges)
dev.off()

### 4.5/ Plot updated densityMaps mapping trait evolution for a given 'focal_time' ####

# Root age
max(phytools::nodeHeights(tree = Fish_tree))

# Plot initial densityMaps (t = 0)
densityMaps_0My <- Fish_deepSTRAPP_biogeo_0_100$updated_trait_data_with_Map_over_time[[1]]
plot_densityMaps_overlay(densityMaps_0My$densityMaps,
                         colors_per_levels = colors_per_ranges,
                         fsize = 0.1) # Reduce tip label size
title(main = "Trait evolution for 250-0 My")

# Plot updated densityMaps for time-step n°3 = 10 My
densityMaps_20My <- Fish_deepSTRAPP_biogeo_0_100$updated_trait_data_with_Map_over_time[[3]]
plot_densityMaps_overlay(densityMaps_20My$densityMaps,
                         colors_per_levels = colors_per_ranges,
                         fsize = 0.1) # Reduce tip label size
title(main = "Trait evolution for 250-20 My")

# Plot updated densityMaps for time-step n°11 = 50 My
densityMaps_50My <- Fish_deepSTRAPP_biogeo_0_100$updated_trait_data_with_Map_over_time[[11]]
plot_densityMaps_overlay(densityMaps_50My$densityMaps,
                         colors_per_levels = colors_per_ranges,
                         fsize = 0.1) # Reduce tip label size
title(main = "Trait evolution for 250-50 My")


### 4.6/ Plot updated diversification rates and regimes for a given 'focal_time' ####

# Extract root age
root_age <- max(phytools::nodeHeights(tree = Fish_tree))

# Plot diversification rates on initial phylogeny (t = 0)
BAMM_map_0My <- Fish_deepSTRAPP_biogeo_0_100$updated_BAMM_objects_over_time[[1]]
plot_BAMM_rates(BAMM_map_0My, labels = FALSE, par.reset = FALSE)
abline(v = root_age - 10, col = "red", lty = 2) # Show where the phylogeny will be cut at 10 Mya
abline(v = root_age - 50, col = "red", lty = 2) # Show where the phylogeny will be cut at 50 Mya
title(main = "BAMM rates for 250-0 My")

# Plot diversification rates on updated phylogeny for time-step n°3 = 10 My
BAMM_map_10My <- Fish_deepSTRAPP_biogeo_0_100$updated_BAMM_objects_over_time[[3]]
plot_BAMM_rates(BAMM_map_10My, labels = FALSE,
                colorbreaks = BAMM_map_10My$initial_colorbreaks$net_diversification)
title(main = "BAMM rates for 250-10 My")

# Plot diversification rates on updated phylogeny for time-step n°11 = 50 My
BAMM_map_50My <- Fish_deepSTRAPP_biogeo_0_100$updated_BAMM_objects_over_time[[11]]
plot_BAMM_rates(BAMM_map_50My, labels = FALSE,
                colorbreaks = BAMM_map_50My$initial_colorbreaks$net_diversification)
title(main = "BAMM rates for 250-50 My")


### 4.7/ Combine both mapped phylogenies with trait evolution (4.5) and diversification rates and regimes (4.6) ####

# Plot both mapped phylogenies in the present (t = 0)
plot_traits_vs_rates_on_phylogeny_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 0,
  ftype = "off", lwd = 0.7,
  colors_per_levels = colors_per_ranges,
  labels = FALSE, legend = FALSE,
  par.reset = FALSE)

# Plot both mapped phylogenies for time-step n°3 = 10 My
plot_traits_vs_rates_on_phylogeny_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 10,
  ftype = "off", lwd = 0.7,
  colors_per_levels = colors_per_ranges,
  labels = FALSE, legend = FALSE,
  par.reset = FALSE)

# Plot both mapped phylogenies for time-step n°11 = 50 My
plot_traits_vs_rates_on_phylogeny_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_biogeo_0_100,
  focal_time = 50,
  ftype = "off", lwd = 0.7,
  colors_per_levels = colors_per_ranges,
  labels = FALSE, legend = FALSE,
  par.reset = FALSE)

