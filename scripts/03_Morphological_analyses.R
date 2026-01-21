
##### Script 03: Morphological analyses #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Run deepSTRAPP on morphological data

###

### Inputs

# Time-calibrated phylogeny
# contMap of morphological traits

###

### Outputs

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

# ### 1.3/ Load contMap of morphological rates data ####
# 
# # Evolutionary rates of Maximum Body Length (sigma² of BM)
# mean_contMap_MBL_rates <- readRDS(file = "./outputs/Morpho/mean_contMap_MBL_rates.rds")
# 
# ## Do NOT use rates! Need trait values for a STRAPP tests !

### 1.3/ Load morphological tip data

# Load trait data
morpho_tip_data_df <- readRDS(file = "./outputs/Morpho/morpho_tip_data_df.rds")

### 1.4/ Load BAMM results for diversification

# Load BAMMdata_object_diversification
BAMMdata_object_diversification <- readRDS(file = "./outputs/BAMM/BAMMdata_object_diversification.rds")
length(BAMMdata_object_diversification$eventData)

##### 2/ Run deepSTRAPP workflow on trait values #####

### 2.1/ Model ancestral trait values ####

# Extract continuous trait data as a named vector
Fish_MBL_tip_data <- setNames(object = morpho_tip_data_df$MBL,
                              nm = morpho_tip_data_df$Taxa)

# For continuous trait, a BM model is assumed by default.
Fish_MBL_ancestral_trait_object <- prepare_trait_data(
   tip_data = Fish_MBL_tip_data,
   trait_data_type = "continuous",
   evolutionary_models = c("BM","OU","EB","rate_trend","lambda","kappa","delta"),
   # evolutionary_models = c("BM","EB","rate_trend", "lambda","kappa","delta"),
   phylo = Fish_tree,
   return_best_model_fit = TRUE,
   return_model_selection_df = TRUE,
   seed = 1234) # Set seed for reproducibility

# Explore output
str(Fish_MBL_ancestral_trait_object, 1)

# Explore summary of model selection
Fish_MBL_ancestral_trait_object$model_selection_df

# Explore best model output
Fish_MBL_ancestral_trait_object$best_model_fit

# Select a color scheme from lowest to highest values
color_scale = c("darkgreen", "limegreen", "orange", "red")
color_scale = c("grey", "red", "gold", "limegreen", "forestgreen", "dodgerblue", "darkblue")

# Extract the contMap representing continuous trait evolution on the phylogeny
Fish_MBL_contMap <- Fish_MBL_ancestral_trait_object$contMap
plot_contMap(Fish_MBL_contMap, color_scale = color_scale)

# Extract the Ancestral Character Estimates (ACE) = trait values at nodes
Fish_MBL_ACE <- Fish_MBL_ancestral_trait_object$ace
head(Fish_MBL_ACE)

## Save ancestral trait values
saveRDS(object = Fish_MBL_ancestral_trait_object, file = "./outputs/Morpho/Fish_MBL_ancestral_trait_object.rds")

## Export summary of model selection
write.csv(x = Fish_MBL_ancestral_trait_object$model_selection_df, file = "./outputs/Morpho/MBL_model_selection_df.csv")

### 2.2/ Compute STRAPP tests over time ####

## Load ancestral trait values
Fish_MBL_ancestral_trait_object <- readRDS(file = "./outputs/Morpho/Fish_MBL_ancestral_trait_object.rds")

## Set for five time steps of 5 My. Will generate deepSTRAPP workflows for 0 to 100 Mya.
time_step_duration <- 5
time_range <- c(0, 100)

# Run deepSTRAPP on net diversification rates
## This step is time-consuming. You can skip it and load directly the result if needed
Fish_deepSTRAPP_MBL_0_100 <- run_deepSTRAPP_over_time(
  contMap = Fish_MBL_ancestral_trait_object$contMap,
  ace = Fish_MBL_ancestral_trait_object$ace,
  tip_data = Fish_MBL_tip_data,
  trait_data_type = "continuous",
  BAMM_object = BAMMdata_object_diversification,
  time_range = time_range,
  time_step_duration = time_step_duration,
  seed = 1234, # Set seed for reproducibility
  alpha = 0.05,
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
str(Fish_deepSTRAPP_MBL_0_100, max.level = 1)

# Display test summary
# Can be passed down to [deepSTRAPP::plot_STRAPP_pvalues_over_time()] to generate a plot
# showing the evolution of the test results across time
Fish_deepSTRAPP_MBL_0_100$pvalues_summary_df

# Access STRAPP test results
# Can be passed down to [deepSTRAPP::plot_histograms_STRAPP_tests_over_time()] to generate plot
# showing the null distribution of the test statistics
str(Fish_deepSTRAPP_MBL_0_100$STRAPP_results, max.level = 2)

# Save results
saveRDS(Fish_deepSTRAPP_MBL_0_100, file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_0_100.rds")


##### 4/ Plot results #####

# Load results
Fish_deepSTRAPP_MBL_0_100 <- readRDS(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_0_100.rds")

### 4.1/ Plot evolution of STRAPP tests p-values through time ####

## Overall Spearman tests

pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_0_80_pvalues_over_time_ggplot.pdf", width = 8, height = 6)
deepSTRAPP::plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  time_range = c(0, 80), pvalues_max = 1.0,
  alpha = 0.05)
dev.off()


### 4.2/ Plot histogram of STRAPP test stats ####

## Overall Spearman tests

pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_0My_stats_histo_ggplot.pdf", width = 8, height = 6)
# Plot the histogram of overall Spearman stats for time-step n°1 = 0 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 0)
dev.off()

# Plot the histogram of overall Spearman stats for time-step n°5 = 20 My
pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_20My_stats_histo_ggplot.pdf", width = 8, height = 6)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 20)
dev.off()

# Plot the histogram of overall Spearman stats for time-step n°9 = 40 My
pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_40My_stats_histo_ggplot.pdf", width = 8, height = 6)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 40)
dev.off()

# Plot the histogram of overall Spearman stats for time-step n°13 = 60 My
pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_60My_stats_histo_ggplot.pdf", width = 8, height = 6)
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 60)
dev.off()

### 4.3/ Plot evolution of rates though time in relation to trait values ####

pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_0_100_RTT_ggplot_Spectral.pdf", width = 11, height = 6)
plot_rates_through_time(deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
                        # color_scale = c("darkgreen", "limegreen", "orange", "red"),
                        color_scale = NULL, # Spectral
                        time_range = c(0, 80),
                        plot_CI = TRUE)
dev.off()

### 4.4/ Plot rates vs. states across branches for a given 'focal_time' ####

# # Select a color scheme from lowest to highest values
# color_scale = c("darkgreen", "limegreen", "orange", "red")
# Get 4 colors from Spectral
color_scale <- rev(RColorBrewer::brewer.pal(n = 4, name = "Spectral"))

library(showtext)
showtext_auto()  # Enable the use of Unicode with Helvetica when plotting

# Generate ggplot for time = 0 My
pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_rates_vs_traits_0My_ggplot.pdf", width = 8, height = 6)
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 0,
  color_scale = color_scale)
dev.off()

# Generate ggplot for time = 20 My
pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_rates_vs_traits_20My_ggplot.pdf", width = 8, height = 6)
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 20,
  color_scale = color_scale)
dev.off()

# Generate ggplot for time = 40 My
pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_rates_vs_traits_40My_ggplot.pdf", width = 8, height = 6)
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 40,
  color_scale = color_scale)
dev.off()

# Generate ggplot for time = 60 My
pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_rates_vs_traits_60My_ggplot.pdf", width = 8, height = 6)
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 60,
  color_scale = color_scale)
dev.off()

### 4.5/ Plot updated densityMaps mapping trait evolution for a given 'focal_time' ####

# Root age
max(phytools::nodeHeights(tree = Fish_tree))

# Plot initial contMap (t = 0)
contMap_0My <- Fish_deepSTRAPP_MBL_0_100$updated_trait_data_with_Map_over_time[[1]]
plot_contMap(contMap_0My$contMap,
             color_scale = color_scale,
             fsize = c(0.1, 1)) # Reduce tip label size
title(main = "Trait evolution for 250-0 My")

# Plot updated densityMaps for time-step n°5 = 20 My
contMap_10My <- Fish_deepSTRAPP_MBL_0_100$updated_trait_data_with_Map_over_time[[5]]
plot_contMap(contMap_10My$contMap,
             color_scale = color_scale,
             fsize = c(0.1, 1)) # Reduce tip label size
title(main = "Trait evolution for 250-20 My")

# Plot updated densityMaps for time-step n°9 = 40 My
contMap_50My <- Fish_deepSTRAPP_MBL_0_100$updated_trait_data_with_Map_over_time[[9]]
plot_contMap(contMap_50My$contMap,
             color_scale = color_scale,
             fsize = c(0.1, 1)) # Reduce tip label size
title(main = "Trait evolution for 250-40 My")


### 4.6/ Plot updated diversification rates and regimes for a given 'focal_time' ####

# Extract root age
root_age <- max(phytools::nodeHeights(tree = Fish_tree))

# Plot diversification rates on initial phylogeny (t = 0)
BAMM_map_0My <- Fish_deepSTRAPP_MBL_0_100$updated_BAMM_objects_over_time[[1]]
plot_BAMM_rates(BAMM_map_0My, labels = FALSE, par.reset = FALSE)
abline(v = root_age - 10, col = "red", lty = 2) # Show where the phylogeny will be cut at 10 Mya
abline(v = root_age - 50, col = "red", lty = 2) # Show where the phylogeny will be cut at 50 Mya
title(main = "BAMM rates for 250-0 My")

# Plot diversification rates on updated phylogeny for time-step n°3 = 10 My
BAMM_map_10My <- Fish_deepSTRAPP_MBL_0_100$updated_BAMM_objects_over_time[[3]]
plot_BAMM_rates(BAMM_map_10My, labels = FALSE,
                colorbreaks = BAMM_map_10My$initial_colorbreaks$net_diversification)
title(main = "BAMM rates for 250-10 My")

# Plot diversification rates on updated phylogeny for time-step n°11 = 50 My
BAMM_map_50My <- Fish_deepSTRAPP_MBL_0_100$updated_BAMM_objects_over_time[[11]]
plot_BAMM_rates(BAMM_map_50My, labels = FALSE,
                colorbreaks = BAMM_map_50My$initial_colorbreaks$net_diversification)
title(main = "BAMM rates for 250-50 My")


### 4.7/ Combine both mapped phylogenies with trait evolution (4.5) and diversification rates and regimes (4.6) ####

# Plot both mapped phylogenies in the present (t = 0)
plot_traits_vs_rates_on_phylogeny_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 0,
  ftype = "off", lwd = 0.7,
  color_scale = c("darkgreen", "limegreen", "orange", "red"),
  labels = FALSE, legend = FALSE,
  par.reset = FALSE)

# Plot both mapped phylogenies for time-step n°3 = 10 My
plot_traits_vs_rates_on_phylogeny_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 10,
  ftype = "off", lwd = 0.7,
  color_scale = c("darkgreen", "limegreen", "orange", "red"),
  labels = FALSE, legend = FALSE,
  par.reset = FALSE)

# Plot both mapped phylogenies for time-step n°11 = 50 My
plot_traits_vs_rates_on_phylogeny_for_focal_time(
  deepSTRAPP_outputs = Fish_deepSTRAPP_MBL_0_100,
  focal_time = 50,
  ftype = "off", lwd = 0.7,
  color_scale = c("darkgreen", "limegreen", "orange", "red"),
  labels = FALSE, legend = FALSE,
  par.reset = FALSE)


##### 5/ Plot rates and traits on phylogenies with deepSTRAPP plotting options #####

# Load results
Fish_deepSTRAPP_MBL_0_100 <- readRDS(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_0_100.rds")

### 5.1/ Plot morphological data ####

## Try with contMap in base R

contMap_0My <- Fish_deepSTRAPP_MBL_0_100$updated_trait_data_with_Map_over_time[[1]]$contMap

# # Select a color scheme from lowest to highest values
# color_scale = c("darkgreen", "limegreen", "orange", "red")
# Get 4 colors from Spectral
color_scale <- rev(RColorBrewer::brewer.pal(n = 4, name = "Spectral"))

# Plot contMap
pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_contMap_Spectral.pdf", width = 20, height = 20)
plot_contMap(contMap = contMap_0My,
             color_scale = color_scale,
             type = "fan", # part = 0.9, open.angle = 340
             lwd = 1, fsize = c(0.2, 1))
dev.off()

# Define geological intervals in Myr (from present backward)
epochs <- data.frame(
  start = c(0, 23, 66, 145),   # younger limit (Mya)
  end   = c(23, 66, 145, 252), # older limit (Mya)
  col   = c("gray95", "gray90", "gray85", "gray80"),
  label = c("Neogene", "Paleogene", "Cretaceous", "Jurassic")
)

# Re-plot with annuli in background
# (Find max distance in plot from center to tip)
usr <- par("usr")
max_r <- max(abs(usr[1:2]), abs(usr[3:4]))

# Define function to draw annuli
draw_ring <- function(r_in, r_out, col, border = NA, n = 360) {
  theta <- seq(0, 2 * pi, length.out = n)
  x_outer <- r_out * cos(theta)
  y_outer <- r_out * sin(theta)
  x_inner <- rev(r_in * cos(theta))
  y_inner <- rev(r_in * sin(theta))
  polygon(c(x_outer, x_inner), c(y_outer, y_inner), col = col, border = border)
}

# Draw rings for each epoch
for (i in seq_len(nrow(epochs))) {
  r_in <- max_r * (epochs$start[i] / max(epochs$end)) # younger boundary
  r_out <- max_r * (epochs$end[i] / max(epochs$end))  # older boundary
  draw_ring(r_in, r_out, col = epochs$col[i])
}

# Replot tree on top
plot_contMap(contMap = contMap_0My,
             color_scale = color_scale,
             type = "fan",
             lwd = 1, fsize = c(0.2, 1),
             add = TRUE)

# Optional: Add labels for epochs
text(
  x = 0,
  y = seq(-max_r * 0.85, -max_r * 0.2, length.out = nrow(epochs)),
  labels = epochs$label,
  cex = 0.8,
  col = "gray30"
)


### 5.2/ Plot diversification data ####

BAMM_object_0My <- Fish_deepSTRAPP_MBL_0_100$updated_BAMM_objects_over_time[[1]]

## Extract Marginal Shift Probability of each branch and scale branch length accordingly
MSP_tree <- BAMMtools::marginalShiftProbsTree(BAMM_object_0My)
BAMM_object_0My$MSP_tree <- MSP_tree

## Extract the Maximum A Posteriori probability (MAP) configuration = the configuration of shift location showing up the most in the posterior sample
# Ignore shifts that have an odd-ratio of marginal posterior probability / prior < 'MAP_odd_ratio_threshold' to avoid noise from non-core shifts
# Rates are then averaged across all samples with the most frequent shift configuration of core-shifts

# Detect MAP configurations
MAP_detection <- BAMMtools::credibleShiftSet(ephy = BAMM_object_0My,
                                             expectedNumberOfShifts = 1, # Default in BAMMtools::setBAMMprior()
                                             threshold = 5, # Odd-ratio threshold used to select core-shifts used to compare configurations
                                             set.limit = 0.95)
# Extract indices of MAP samples
BAMM_object_0My$MAP_indices <- MAP_detection$indices[[1]]

# Compute mean rates/regimes across MAP samples
MAP_BAMM_object <- BAMMtools::getBestShiftConfiguration(BAMM_object_0My,
                                                        expectedNumberOfShifts = 1, # Default in BAMMtools::setBAMMprior()
                                                        threshold = 5) # Odd-ratio threshold used to select core-shifts used to compare configurations

# Reorder elements to fit order in the main BAMM_object
if ("node.label" %in% names(MAP_BAMM_object))
{
  MAP_BAMM_object <- MAP_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length", "node.label",
                                       "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                       "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                       "meanTipLambda", "meanTipMu", "type")]
} else {
  MAP_BAMM_object <- MAP_BAMM_object[c("edge", "Nnode", "tip.label", "edge.length",
                                       "begin", "end", "downseq", "lastvisit", "numberEvents", "eventData",
                                       "eventVectors", "tipStates", "tipLambda", "tipMu", "eventBranchSegs",
                                       "meanTipLambda", "meanTipMu", "type")]
}
class(MAP_BAMM_object) <- "bammdata"
attr(x = MAP_BAMM_object, which = "order") <- "cladewise"

BAMM_object_0My$MAP_BAMM_object <- MAP_BAMM_object

pdf(file = "./outputs/Morpho/Fish_deepSTRAPP_BAMM_rates_on_phylo.pdf", width = 20, height = 20)
plot_BAMM_rates(BAMM_object = BAMM_object_0My,
                method = "polar", labels = TRUE,
                legend = TRUE, pal = "RdYlBu",
                breaksmethod = "jenks",
                lwd = 1, cex = 0.2)
dev.off()


##### 6/ Plot rates and traits on phylogenies using ggtree #####

# Load results
Fish_deepSTRAPP_MBL_0_100 <- readRDS(file = "./outputs/Morpho/Fish_deepSTRAPP_MBL_0_100.rds")
