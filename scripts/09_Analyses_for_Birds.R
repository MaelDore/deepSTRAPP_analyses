
##### Script 09: Run analyses on Rabosky & Huang, 2016 (STRAPP Paper) Bird data #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Curate taxonomic information from species list to prepare for analyses
# Extract species-level information for morphological data
# Load BAMM results
# Run deepSTRAPP on bird dichromatism data

###

### Inputs

# Time-calibrated phylogeny
# BAMM objects
# Plumage dichromatism data

###

### Sources 

## Phylogeny and Ecological data (Secondary sources)

# Rabosky, D. L., & Huang, H. (2016). A robust semi-parametric test for detecting trait-dependent diversification. Systematic biology, 65(2), 181-193.
# Huang, H., & Rabosky, D. L. (2014). Sexual selection and diversification: reexamining the correlation between dichromatism and speciation rate in birds. The American Naturalist, 184(5), E101-E114.

## Primary sources:

# For the trait data:
# Armenta, J. K., P. O. Dunn, and L. A. Whittingham. 2008. Quantifying avian sexual dichromatism: a comparison of methods. Journal of Experimental Biology 211:2423–2430.

# For the phylogeny:
# Jetz, W., G. H. Thomas, J. B. Joy, K. Hartmann, and A. O. Mooers. 2012. The global diversity of birds in space and time. Nature 491: 444–448.

###

### Outputs

# Pruned phylogeny
# Matching for ecological data at species level
# BAMM object with 1000 posterior samples
# deepSTRAPP run on dichromatism variables (3 continuous variables)

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff #####

### 1.1/ Load packages ####

library(tidyverse)
library(readxl)
library(xlsx)      # Need the Java Development Kit (JDK) installed
library(openxlsx)  # Use Rccp. No need of Java
library(ape)
library(BAMMtools)
library(deepSTRAPP)
library(Rphylopars)
library(deeptime)
library(ggtree)

# BiocManager::install("ggtree")

# # Install non-CRAN dependencies
# remotes::install_github(repo = "nmatzke/BioGeoBEARS")
# remotes::install_github(repo = "bstaggmartin/contsimmap")

# # Install the latest deepSTRAPP version
# remotes::install_github(repo = "MaelDore/deepSTRAPP",
#                         # Time-consuming, but needed if you want to have access to the vignettes/tutorials
#                         build_vignettes = FALSE)


### 1.2/ Load Ecological trait data ####

Data_df <-  readr::read_table("./input_data/STRAPP_paper/Armenta_data.txt")

summary(Data_df)

## Clean variable names

# PCA = Euclidean distances in PCA on reflectance curves
# segclass = Euclidean distances in HSV space extracted from reflectance curves
# colordiscrim = Differences in cone stimulations from UV-tuned avian vision model of male/female images on deciduous forest background

Bird_data_df <- Data_df %>% 
  dplyr::rename(Taxa_labels = `"species"`,
                Reflectance_dist = `"pca"`,
                HSV_dist = `"segclass"`,
                Stimulation_diff = `"colordiscrim"`) %>%
  dplyr::select(Taxa_labels, Reflectance_dist, HSV_dist, Stimulation_diff)
Bird_data_df$Taxa_labels <- str_remove_all(string = Bird_data_df$Taxa_labels, pattern = '\"')

# Check corrrelations

plot(x = Bird_data_df$Reflectance_dist, y = Bird_data_df$HSV_dist)
plot(x = Bird_data_df$Reflectance_dist, y = Bird_data_df$Stimulation_diff)
plot(x = Bird_data_df$HSV_dist, y = Bird_data_df$Stimulation_diff)

# Check missing data
table(is.na(Bird_data_df$Reflectance_dist)) # Missing 16
table(is.na(Bird_data_df$HSV_dist)) # Missing 2
table(is.na(Bird_data_df$Stimulation_diff)) # Missing 24

table(is.na(Bird_data_df$Reflectance_dist) | is.na(Bird_data_df$HSV_dist) | is.na(Bird_data_df$Stimulation_diff)) # 26 species with missing data in total

## Save object
saveRDS(object = Bird_data_df, file = "./input_data/STRAPP_paper/Bird_data_df.rds")


### 1.3/ Load phylogeny

Bird_tree_full <- ape::read.tree(file = "./input_data/STRAPP_paper/birds_hackett.tre")

# Inspect the three options
plot(Bird_tree_full) # Visually ultrametric
length(Bird_tree_full$tip.label) # 6670 species
ape::is.ultrametric(Bird_tree_full) # Is ultrametric
range(Bird_tree_full$edge.length) # All positive branch lengths

## Save phylo
saveRDS(object = Bird_tree_full, file = "./input_data/STRAPP_paper/Bird_tree_full.rds")


##### 2/ Match ecological data with phylogeny #####

Bird_data_df <- readRDS(file = "./input_data/STRAPP_paper/Bird_data_df.rds")
Bird_tree_full <- readRDS(file = "./input_data/STRAPP_paper/Bird_tree_full.rds")

### 2.1/ Prune the phylogeny and trait_df ####

# 979 species with stimulation data available
table(is.na(Bird_data_df$Stimulation_diff))
summary(Bird_data_df$Stimulation_diff)
sd(Bird_data_df$Stimulation_diff, na.rm = T)

# 977 species with all morpho data available
table(is.na(Bird_data_df$Reflectance_dist) | is.na(Bird_data_df$HSV_dist) | is.na(Bird_data_df$Stimulation_diff)) # 26 species with missing data in total

Morpho_labels <- Bird_data_df$Taxa_labels[!(is.na(Bird_data_df$Reflectance_dist) | is.na(Bird_data_df$HSV_dist) | is.na(Bird_data_df$Stimulation_diff))]
Phylo_labels <- Bird_tree_full$tip.label

# 895 species shared between the phylo and morpho dataset
Shared_species <- intersect(Phylo_labels, Morpho_labels)

# Prune Bird_data_df
Bird_data_df_pruned <- Bird_data_df %>%
  dplyr::filter(Taxa_labels %in% Shared_species)

# Prune the phylogeny
Bird_tree_pruned <- ape::keep.tip(phy = Bird_tree, tip = Shared_species)

# Reorder as in phylo
Bird_data_df_pruned <- Bird_data_df_pruned[match(x = Bird_tree_pruned$tip.label, table = Bird_data_df_pruned$Taxa_labels), ]

# Save objects
saveRDS(object = Bird_data_df_pruned, file = "./input_data/STRAPP_paper/Bird_data_df_pruned.rds")
saveRDS(object = Bird_tree_pruned, file = "./input_data/STRAPP_paper/Bird_tree_pruned.rds")



##### 3/ Build contMap of trait evolution #####

# Load phylo and Bird_data_df_pruned
Bird_tree_pruned <- readRDS(file = "./input_data/STRAPP_paper/Bird_tree_pruned.rds")
Bird_data_df_pruned <- readRDS(file = "./input_data/STRAPP_paper/Bird_data_df_pruned.rds")

### 3.1/ For Reflectance_dist ####

## Prepare trait data

hist(Bird_data_df_pruned$Reflectance_dist)
range(Bird_data_df_pruned$Reflectance_dist)
hist(log1p(Bird_data_df_pruned$Reflectance_dist))

# Almost centered. But we don't know the initial mean to rebuild raw data...
summary(Bird_data_df_pruned$Reflectance_dist) 

# Use log+1 transformation. Not sure how distance could have been negative in the first place...

# Extract continuous trait data as a named vector
Reflectance_dist_tip_data <- setNames(
   # object = Bird_data_df_pruned$Reflectance_dist,
   object = log1p(Bird_data_df_pruned$Reflectance_dist),
   nm = Bird_tree_pruned$tip.label)
head(Reflectance_dist_tip_data)

# Select a color scheme from lowest to highest values
color_scale = c("darkgreen", "limegreen", "orange", "red")

## Map trait evolution as ML estimates
Reflectance_dist_mapped_data <- prepare_trait_data(
  tip_data = Reflectance_dist_tip_data,
  trait_data_type = "continuous",
  phylo = Bird_tree_pruned,
  seed = 1234,
  evolutionary_models = "BM",
  plot_map = FALSE,
  verbose = TRUE)

## Plot contMap = ML estimates of continuous trait evolution
plot_contMap(contMap = Reflectance_dist_mapped_data$contMap,
             color_scale = color_scale,
             fsize = c(0, 1)) # Remove tip labels

## Save output
saveRDS(object = Reflectance_dist_mapped_data, file = "./outputs/Birds/Reflectance_dist_mapped_data.rds")


### 3.2/ For HSV_dist ####

## Prepare trait data

hist(Bird_data_df_pruned$HSV_dist)
range(Bird_data_df_pruned$HSV_dist)
hist(log1p(Bird_data_df_pruned$HSV_dist))

# Almost centered. But we don't know the initial mean to rebuild raw data...
summary(Bird_data_df_pruned$HSV_dist) 

# Use log+1 transformation. Not sure how distance could have been negative in the first place...


# Extract continuous trait data as a named vector
HSV_dist_tip_data <- setNames(
   # object = Bird_data_df_pruned$HSV_dist,
   object = log1p(Bird_data_df_pruned$HSV_dist),
   nm = Bird_tree_pruned$tip.label)
head(HSV_dist_tip_data)

# Select a color scheme from lowest to highest values
color_scale = c("darkgreen", "limegreen", "orange", "red")

## Map trait evolution as ML estimates
HSV_dist_mapped_data <- prepare_trait_data(
  tip_data = HSV_dist_tip_data,
  trait_data_type = "continuous",
  phylo = Bird_tree_pruned,
  seed = 1234,
  evolutionary_models = "BM",
  plot_map = FALSE,
  verbose = TRUE)

## Plot contMap = ML estimates of continuous trait evolution
plot_contMap(contMap = HSV_dist_mapped_data$contMap,
             color_scale = color_scale,
             fsize = c(0, 1)) # Remove tip labels

## Save output
saveRDS(object = HSV_dist_mapped_data, file = "./outputs/Birds/HSV_dist_mapped_data.rds")


### 3.3/ For Stimulation_diff ####

## Prepare trait data

hist(Bird_data_df_pruned$Stimulation_diff)
range(Bird_data_df_pruned$Stimulation_diff)
hist(log1p(Bird_data_df_pruned$Stimulation_diff))

# Almost centered. But we don't know the initial mean to rebuild raw data...
summary(Bird_data_df_pruned$Stimulation_diff) 

# Use log+1 transformation to ensure the variable is projected on R


# Extract continuous trait data as a named vector
Stimulation_diff_tip_data <- setNames(
   # object = Bird_data_df_pruned$Stimulation_diff,
   object = log1p(Bird_data_df_pruned$Stimulation_diff),
   nm = Bird_tree_pruned$tip.label)
head(Stimulation_diff_tip_data)

# Select a color scheme from lowest to highest values
color_scale = c("darkgreen", "limegreen", "orange", "red")

## Map trait evolution as ML estimates
Stimulation_diff_mapped_data <- prepare_trait_data(
  tip_data = Stimulation_diff_tip_data,
  trait_data_type = "continuous",
  phylo = Bird_tree_pruned,
  seed = 1234,
  evolutionary_models = "BM",
  plot_map = FALSE,
  verbose = TRUE)

## Plot contMap = ML estimates of continuous trait evolution
plot_contMap(contMap = Stimulation_diff_mapped_data$contMap,
             color_scale = color_scale,
             fsize = c(0, 1)) # Remove tip labels

## Save output
saveRDS(object = Stimulation_diff_mapped_data, file = "./outputs/Birds/Stimulation_diff_mapped_data.rds")


##### 4/ Build BAMM objects #####

# Load pruned phylo
Bird_tree_pruned <- readRDS(file = "./input_data/STRAPP_paper/Bird_tree_pruned.rds")

## Load the full tree used initially to run BAMM
Bird_tree_full <- readRDS(file = "./input_data/STRAPP_paper/Bird_tree_full.rds")

## Load initial BAMM object with all posterior samples
Bird_BAMM_object_full <- BAMMtools::getEventData(phy = Bird_tree_full, 
                                              eventdata = "./input_data/STRAPP_paper/hack2vr_event_data_250m.csv",
                                              burnin = 0,
                                              verbose = TRUE)

hist(Bird_BAMM_object$numberEvents)

source("../deepSTRAPP/R/prune_BAMM_object.R")

## Build BAMM object for deepSTRAPP with 1000 posterior samples
Bird_BAMM_object <- deepSTRAPP::build_BAMM_object(
  phylo = Bird_tree_full, 
  eventdata = "./input_data/STRAPP_paper/hack2vr_event_data_250m.csv",
  burn_in = 0.25, nb_posterior_samples = 1000, expectedNumberOfShifts = 34,
  verbose = TRUE)

hist(Bird_BAMM_object$numberEvents)

plot_BAMM_rates(BAMM_object = Bird_BAMM_object)

## Prune the BAMM_object to match with trait data
Bird_BAMM_object_pruned <- prune_BAMM_object(
  BAMM_object = Bird_BAMM_object, 
  tips_to_keep = Bird_tree_pruned$tip.label,
  verbose = TRUE)

hist(Bird_BAMM_object_pruned$numberEvents)

plot_BAMM_rates(BAMM_object = Bird_BAMM_object_pruned)

## Save output
saveRDS(object = Bird_BAMM_object_full, file = "./outputs/Birds/Bird_BAMM_object_full.rds")
saveRDS(object = Bird_BAMM_object_pruned, file = "./outputs/Birds/Bird_BAMM_object_pruned.rds")


##### 5/ Run deepSTRAPP on Reflectance_dist data #####

# Load output mapped trait data
Reflectance_dist_mapped_data <- readRDS(file = "./outputs/Birds/Reflectance_dist_mapped_data.rds")
# Load BAMM_object_pruned summarizing 1000 posterior samples of BAMM
Bird_BAMM_object_pruned <- readRDS(file = "./outputs/Birds/Bird_BAMM_object_pruned.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(Bird_tree_pruned)[,2])
root_age # 113 Mya

## Set for time steps of 5 My from 0 to 50 Mya
# nb_time_steps <- 10
time_step_duration <- 5
time_range <- c(0, 50)


## Run deepSTRAPP on net diversification rates
Birds_deepSTRAPP_Reflectance_dist <- run_deepSTRAPP_over_time(
  contMap = Reflectance_dist_mapped_data$contMap,
  trait_data_type = "continuous",
  BAMM_object = Bird_BAMM_object_pruned,
  seed = 1234,
  # nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  time_range = time_range,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Birds_deepSTRAPP_Reflectance_dist, max.level = 1)

# Display test summary
Birds_deepSTRAPP_Reflectance_dist$pvalues_summary_df

## Visualize test results

# Plot p-values of Spearman tests across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Reflectance_dist,
  alpha = 0.05)

# Plot evolution of mean rates through time
plot_rates_through_time(deepSTRAPP_outputs = Birds_deepSTRAPP_Reflectance_dist,
                        plot_CI = TRUE)

# Plot rates vs. trait values across branches for time step = 10 My
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Reflectance_dist,
  focal_time = 10)

# Plot histogram of Spearman test stats for time step = 10 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Reflectance_dist,
  focal_time = 10)

# Plot histograms of Spearman test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Reflectance_dist,
  display_plots = TRUE)

## Save output
saveRDS(object = Birds_deepSTRAPP_Reflectance_dist, file = "./outputs/Birds/Birds_deepSTRAPP_Reflectance_dist.rds")


##### 6/ Run deepSTRAPP on HSV_dist data #####

# Load output mapped trait data
HSV_dist_mapped_data <- readRDS(file = "./outputs/Birds/HSV_dist_mapped_data.rds")
# Load BAMM_object_pruned summarizing 1000 posterior samples of BAMM
Bird_BAMM_object_pruned <- readRDS(file = "./outputs/Birds/Bird_BAMM_object_pruned.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(Bird_tree_pruned)[,2])
root_age # 113 Mya

## Set for time steps of 5 My from 0 to 50 Mya
# nb_time_steps <- 10
time_step_duration <- 5
time_range <- c(0, 50)


## Run deepSTRAPP on net diversification rates
Birds_deepSTRAPP_HSV_dist <- run_deepSTRAPP_over_time(
  contMap = HSV_dist_mapped_data$contMap,
  trait_data_type = "continuous",
  BAMM_object = Bird_BAMM_object_pruned,
  seed = 1234,
  # nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  time_range = time_range,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Birds_deepSTRAPP_HSV_dist, max.level = 1)

# Display test summary
Birds_deepSTRAPP_HSV_dist$pvalues_summary_df

## Visualize test results

# Plot p-values of Spearman tests across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_HSV_dist,
  alpha = 0.05)

# Plot evolution of mean rates through time
plot_rates_through_time(deepSTRAPP_outputs = Birds_deepSTRAPP_HSV_dist,
                        plot_CI = TRUE)

# Plot rates vs. trait values across branches for time step = 10 My
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_HSV_dist,
  focal_time = 10)

# Plot histogram of Spearman test stats for time step = 10 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_HSV_dist,
  focal_time = 10)

# Plot histograms of Spearman test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_HSV_dist,
  display_plots = TRUE)

## Save output
saveRDS(object = Birds_deepSTRAPP_HSV_dist, file = "./outputs/Birds/Birds_deepSTRAPP_HSV_dist.rds")


##### 7/ Run deepSTRAPP on Stimulation_diff data #####

# Load output mapped trait data
Stimulation_diff_mapped_data <- readRDS(file = "./outputs/Birds/Stimulation_diff_mapped_data.rds")
# Load BAMM_object_pruned summarizing 1000 posterior samples of BAMM
Bird_BAMM_object_pruned <- readRDS(file = "./outputs/Birds/Bird_BAMM_object_pruned.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(Bird_tree_pruned)[,2])
root_age # 113 Mya

## Set for time steps of 2.5 My from 0 to 50 Mya
# nb_time_steps <- 10
time_step_duration <- 2.5
time_range <- c(0, 50)


## Run deepSTRAPP on net diversification rates
Birds_deepSTRAPP_Stimulation_diff <- run_deepSTRAPP_over_time(
  contMap = Stimulation_diff_mapped_data$contMap,
  trait_data_type = "continuous",
  BAMM_object = Bird_BAMM_object_pruned,
  seed = 1234,
  # nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  time_range = time_range,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Birds_deepSTRAPP_Stimulation_diff, max.level = 1)

# Display test summary
Birds_deepSTRAPP_Stimulation_diff$pvalues_summary_df

## Visualize test results

# Plot p-values of Spearman tests across all time-steps
Birds_Stimulation_diff_pvalues_plot <- plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff,
  alpha = 0.05)

pdf(file = "./outputs/Birds/Birds_Stimulation_diff_pvalues_plot.pdf", height = 6, width = 8)
print(Birds_Stimulation_diff_pvalues_plot)
dev.off()


# Plot evolution of mean rates through time
Birds_Stimulation_diff_RTT_plot <- plot_rates_through_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff,
  plot_CI = TRUE)$rates_TT_ggplot

# Move legend
Birds_Stimulation_diff_RTT_plot <- Birds_Stimulation_diff_RTT_plot +
  ggplot2::theme(legend.title = ggplot2::element_text(size  = 12, face = "bold"),
                 legend.text = ggplot2::element_text(size  = 9),
                 legend.key.height = unit(x = 0.8, units = "cm"),
                 legend.position = "inside",
                 legend.position.inside = c(0.20, 0.75))

# Export without geological epochs
pdf(file = "./outputs/Birds/Birds_Stimulation_diff_RTT_plot.pdf", height = 6, width = 8)
print(Birds_Stimulation_diff_RTT_plot)
dev.off()

## Add geological epochs to plots

# Define geological epochs
epochs_df <- data.frame(
  xmin = c(66, 0),         # Left boundaries (older time)
  xmax = c(143, 23),       # Right boundaries (younger time)
  ymin = -Inf, 
  ymax = Inf,
  epoch = factor(c("Cretaceous", "Quaternary-Neogene"), 
                 levels = c("Cretaceous", "Quaternary-Neogene")))

# Create the background layer
rect_layer <- geom_rect(data = epochs_df, 
                        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = epoch), 
                        fill = "grey90",
                        alpha = 1.0,
                        inherit.aes = FALSE) # Critical: isolates this layer from global aesthetics

# Inject it to the front of the existing plot's drawing list
Birds_Stimulation_diff_RTT_plot_with_epochs <- Birds_Stimulation_diff_RTT_plot
Birds_Stimulation_diff_RTT_plot_with_epochs$layers <- c(rect_layer, Birds_Stimulation_diff_RTT_plot_with_epochs$layers)

Birds_Stimulation_diff_pvalues_plot_with_epochs <- Birds_Stimulation_diff_pvalues_plot
Birds_Stimulation_diff_pvalues_plot_with_epochs$layers <- c(rect_layer, Birds_Stimulation_diff_pvalues_plot_with_epochs$layers)

# # Apply the discrete fill colors if your stored plot doesn't have a fill scale yet
# my_plot <- my_plot + 
#   scale_fill_manual(values = c("Cretaceous" = "gray90", 
#                                "Paleogene" = "gray75", 
#                                "Quaternary-Neogene" = "gray90"))

# Display the modified plot
print(Birds_Stimulation_diff_RTT_plot_with_epochs)

# Export with geological epochs
pdf(file = "./outputs/Birds/Birds_Stimulation_diff_RTT_plot_with_epochs.pdf", height = 6, width = 8)
print(Birds_Stimulation_diff_RTT_plot_with_epochs)
dev.off()

pdf(file = "./outputs/Birds/Birds_Stimulation_diff_pvalues_plot_with_epochs.pdf", height = 6, width = 8)
print(Birds_Stimulation_diff_pvalues_plot_with_epochs)
dev.off()

# Plot rates vs. trait values across branches for time step = 0 My
ggplot1 <- plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff,
  focal_time = 0)$rates_vs_trait_ggplot

# Plot rates vs. trait values across branches for time step = 7.5 My
ggplot2 <- plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff,
  focal_time = 7.5)$rates_vs_trait_ggplot

# Plot histogram of Spearman test stats for time step = 0 My
ggplot3 <- plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff,
  focal_time = 0)

# Plot histogram of Spearman test stats for time step = 7.5 My
ggplot4 <- plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff,
  focal_time = 7.5)

# Plot all plots on one faceted plot
pdf(file = "./outputs/Birds/Birds_Stimulation_diff_output_facetted_plots.pdf", height = 12, width = 14)
cowplot::plot_grid(plotlist = list(ggplot1, ggplot2, ggplot3, ggplot4))
dev.off()


# Plot histograms of Spearman test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff,
  display_plots = TRUE)


## Save output
# saveRDS(object = Birds_deepSTRAPP_Stimulation_diff, file = "./outputs/Birds/Birds_deepSTRAPP_Stimulation_diff.rds")
saveRDS(object = Birds_deepSTRAPP_Stimulation_diff, file = "./outputs/Birds/Birds_deepSTRAPP_Stimulation_diff_2.5.rds")

# ## Load output
# Birds_deepSTRAPP_Stimulation_diff <- readRDS(file = "./outputs/Birds/Birds_deepSTRAPP_Stimulation_diff.rds")
# Birds_deepSTRAPP_Stimulation_diff <- readRDS(file = "./outputs/Birds/Birds_deepSTRAPP_Stimulation_diff_2.5.rds")


##### 8/ Run one-tailed deepSTRAPP test on Stimulation_diff data #####

# Load output mapped trait data
Stimulation_diff_mapped_data <- readRDS(file = "./outputs/Birds/Stimulation_diff_mapped_data.rds")
# Load BAMM_object_pruned summarizing 1000 posterior samples of BAMM
Bird_BAMM_object_pruned <- readRDS(file = "./outputs/Birds/Bird_BAMM_object_pruned.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(ape::as.phylo(Bird_BAMM_object_pruned))[,2])
root_age # 113 Mya

## Set for time steps of 5 My from 0 to 50 Mya
# nb_time_steps <- 10
time_step_duration <- 5
time_range <- c(0, 50)


## Run deepSTRAPP on net diversification rates
Birds_deepSTRAPP_Stimulation_diff_one_tailed <- run_deepSTRAPP_over_time(
  contMap = Stimulation_diff_mapped_data$contMap,
  trait_data_type = "continuous",
  BAMM_object = Bird_BAMM_object_pruned,
  seed = 1234,
  # nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  time_range = time_range,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  two_tailed = FALSE,
  one_tailed_hypothesis = "positive",
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Birds_deepSTRAPP_Stimulation_diff_one_tailed, max.level = 1)

# Display test summary
Birds_deepSTRAPP_Stimulation_diff_one_tailed$pvalues_summary_df

## Visualize test results

# Plot p-values of Spearman tests across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff_one_tailed,
  alpha = 0.05)

# Plot evolution of mean rates through time
plot_rates_through_time(deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff_one_tailed,
                        plot_CI = TRUE)

# Plot rates vs. trait values across branches for time step = 10 My
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff_one_tailed,
  focal_time = 10)

# Plot histogram of Spearman test stats for time step = 10 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff_one_tailed,
  focal_time = 10)

# Plot histograms of Spearman test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Birds_deepSTRAPP_Stimulation_diff_one_tailed,
  display_plots = TRUE)

## Save output
saveRDS(object = Birds_deepSTRAPP_Stimulation_diff_one_tailed, file = "./outputs/Birds/Birds_deepSTRAPP_Stimulation_diff_one_tailed.rds")
