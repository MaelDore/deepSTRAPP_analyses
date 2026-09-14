
##### Script 08: Run analyses on Nelsen et al., 2020 Lichen data #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Curate taxonomic information from species list to prepare for analyses
# Extract species-level information for morphological data
# Run BAMM on time-calibrated phylogeny
# Run deepSTRAPP on lichen types

###

### Inputs

# Time-calibrated phylogeny
# BAMM objects
# Fruit type, seed type and growth type data

###

### Sources 

## Phylogeny and Ecological data

# Nelsen, M. P., Lücking, R., Boyce, C. K., Lumbsch, H. T., & Ree, R. H. (2020).
# The macroevolutionary dynamics of symbiotic and phenotypic diversification in lichens. Proceedings of the National Academy of Sciences, 117(35), 21495-21503. 

###

### Outputs

# Pruned phylogeny
# Matching for ecological data at species level
# BAMM object with 1000 posterior samples
# deepSTRAPP run on Fruit type
# deepSTRAPP run on seed type
# deepSTRAPP run on growth type

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

Data_for_corHMM_df <-  readr::read_csv("./input_data/Nelsen_2020/HiSSE/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_macro_forcorhmm.csv")

Data_df <- readr::read_csv("./input_data/Nelsen_2020/Data/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_for_github.csv")

summary(Data_df)
summary(Data_for_corHMM_df)

# Thalle size data
table(Data_df$Micro_vs_Macro) # Good balance
table(Data_for_corHMM_df$Macro) # Clean version

# Interaction data
table(Data_df$Bipartite_vs_Tripartite) # Too few tripartite

# Photobionts data
table(Data_df$Primary_Grn_vs_Cyano) # Too few Cyano
table(Data_df$PrimaryCyanobacteria) # Too few Cyano
table(Data_df$PrimaryTrentepohliales) # Too few Trentophiales
table(Data_df$PrimaryTrebouxiophyceae) # Could work
table(Data_df$PrimarySapro.Moss) # Too few zeros
table(Data_df$NonLichenized) # Too few ones

Data_df$PrimaryPhotobiont <- NA
Data_df$PrimaryPhotobiont[Data_df$PrimaryTrebouxiophyceae == 1] <- "Trebouxiophyceae"
Data_df$PrimaryPhotobiont[Data_df$PrimaryTrentepohliales == 1] <- "Trentepohliales"
Data_df$PrimaryPhotobiont[Data_df$PrimaryCyanobacteria == 1] <- "Cyanobacteria"
Data_df$PrimaryPhotobiont[Data_df$NonLichenized == 1] <- "NonLichenized"

table(Data_df$PrimaryPhotobiont) # Too unbalanced?

## Make a clean df
Lichen_data_df <- Data_for_corHMM_df %>% 
  dplyr::left_join(y = Data_df[ , c("edited_name", "PrimaryPhotobiont")]) %>%
  dplyr::rename(Taxa_label = edited_name,
                Growth_form = Macro,
                Primary_Photobiont = PrimaryPhotobiont)
Lichen_data_df$Growth_form <- c("Micro", "Macro")[Lichen_data_df$Growth_form + 1]


table(is.na(Lichen_data_df$Growth_form)) # All species covered
table(is.na(Lichen_data_df$Primary_Photobiont)) # Missing 10 species

saveRDS(object = Lichen_data_df, file = "./input_data/Nelsen_2020/Lichen_data_df.rds")


### 1.3/ Load phylogeny

Lichen_tree_ML <- ape::read.tree(file = "./input_data/Nelsen_2020/Trees/Lecanoromycetes_ML_Tree_Timescaled")
Lichen_tree_ML_with_BS <- ape::read.tree(file = "./input_data/Nelsen_2020/Trees/Lecanoromycetes_ML_w_Bootstrap_Props")
Lichen_tree_for_BAMM <- ape::read.tree(file = "./input_data/Nelsen_2020/BAMM/bamm_medusa_med_families_65m/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")

# Inspect the three options
plot(Lichen_tree_ML) # Visually ultrametric
length(Lichen_tree_ML$tip.label) # 3373
ape::is.ultrametric(Lichen_tree_ML) # Not ultrametric
range(Lichen_tree_ML$edge.length) # All positive branch lengths

plot(Lichen_tree_ML_with_BS) # Visually ultrametric
length(Lichen_tree_ML_with_BS$tip.label) # 3373
ape::is.ultrametric(Lichen_tree_ML_with_BS) # Not ultrametric
range(Lichen_tree_ML_with_BS$edge.length) # All positive branch lengths

plot(Lichen_tree_for_BAMM) # Visually ultrametric
length(Lichen_tree_for_BAMM$tip.label) # 3373
ape::is.ultrametric(Lichen_tree_for_BAMM) # Not ultrametric
range(Lichen_tree_for_BAMM$edge.length) # All positive branch lengths


# Same species
table(Lichen_tree_ML_with_BS$tip.label %in% Lichen_tree_ML$tip.label)
table(Lichen_tree_for_BAMM$tip.label %in% Lichen_tree_ML$tip.label)
# All matching the trait df
table(Lichen_tree_for_BAMM$tip.label %in% Lichen_data_df$Taxa_label)

# Different root age for the _with_BS
root_age_ML <- max(phytools::nodeHeights(Lichen_tree_ML)[,2])
root_age_ML # 248 Mya

root_age_ML_with_BS <- max(phytools::nodeHeights(Lichen_tree_ML_with_BS)[,2])
root_age_ML_with_BS # 273 Mya

root_age_for_BAMM <- max(phytools::nodeHeights(Lichen_tree_for_BAMM)[,2])
root_age_for_BAMM # 248 Mya

# Different branch length _with_BS
table(Lichen_tree_ML_with_BS$edge.length[order(Lichen_tree_ML_with_BS$edge.length)] == Lichen_tree_ML$edge.length[order(Lichen_tree_ML$edge.length)])
table(Lichen_tree_ML_with_BS$edge.length[order(Lichen_tree_ML_with_BS$edge.length)] == Lichen_tree_ML$edge.length[order(Lichen_tree_ML$edge.length)])
table(Lichen_tree_for_BAMM$edge.length[order(Lichen_tree_for_BAMM$edge.length)] == Lichen_tree_ML$edge.length[order(Lichen_tree_ML$edge.length)])

# Conclusion: Use the BAMM tree!

Lichen_tree <- Lichen_tree_for_BAMM

## Save phylo
saveRDS(object = Lichen_tree, file = "./input_data/Nelsen_2020/Lichen_tree.rds")


##### 2/ Match ecological data with phylogeny #####

Lichen_data_df <- readRDS(file = "./input_data/Nelsen_2020/Lichen_data_df.rds")
Lichen_tree <- readRDS(file = "./input_data/Nelsen_2020/Lichen_tree.rds")

table(!is.na(Lichen_data_df$Growth_form))
table(!is.na(Lichen_data_df$Primary_Photobiont))

# ### 2.1/ Prune the phylogeny and trait_df ####
# 
# # Does not work because the BAMM object was modeled on the full 3373t phylogeny. Rather, we need to infer missing data!
# 
# table(is.na(Lichen_data_df$Primary_Photobiont)) # Missing 10 species
# 
# Morpho_labels <- Lichen_data_df$Taxa_label[!is.na(Lichen_data_df$Primary_Photobiont)]
# Phylo_labels <- Lichen_tree$tip.label
# 
# Shared_species <- intersect(Phylo_labels, Morpho_labels)
# # All species in the phylogeny are in the ecological dataset!
# 
# # Prune Lichen_data_df
# Lichen_data_df_pruned <- Lichen_data_df %>% 
#   dplyr::filter(Taxa_label %in% Shared_species)
# 
# # Prune the phylogeny
# Lichen_tree_pruned <- ape::keep.tip(phy = Lichen_tree, tip = Shared_species)
# 
# # Reorder as in phylo
# Lichen_data_df_pruned <- Lichen_data_df_pruned[match(x = Lichen_tree_pruned$tip.label, table = Lichen_data_df_pruned$Taxa_label), ]
# 
# # Save objects
# saveRDS(object = Lichen_data_df_pruned, file = "./input_data/Nelsen_2020/Lichen_data_df_pruned.rds")
# saveRDS(object = Lichen_tree_pruned, file = "./input_data/Nelsen_2020/Lichen_tree_pruned.rds")


## 2.2/ Infer missing data for Primary Photobiont ####

# Build named vector tip data
Primary_Photobiont_data <- setNames(object = Lichen_data_df$Primary_Photobiont, nm = Lichen_data_df$Taxa_label)

# Extract states
Primary_Photobiont_states <- sort(unique(na.omit(Primary_Photobiont_data)))

# Build prior matrix
Primary_Photobiont_prior_matrix <- phytools::to.matrix(x = Primary_Photobiont_data, seq = Primary_Photobiont_states)
# Apply uniform prior on NA tip data
Primary_Photobiont_prior_matrix[is.na(Primary_Photobiont_data), ] <- 1 / length(Primary_Photobiont_states)

# Run ACE inference + stochastic mapping
Primary_Photobiont_simmaps <- phytools::make.simmap(
  tree = Lichen_tree,
  x = Primary_Photobiont_prior_matrix, # Provide the prior matrix as tip data
  model = "ARD",
  nsim = 100)

str(Primary_Photobiont_simmaps, 1)
str(Primary_Photobiont_simmaps[[1]], 1)

## Retrieve posterior probabilities for all tips

?phytools::getStates

# Extract tip states across all simmaps
tip_states_list <- lapply(X = Primary_Photobiont_simmaps, FUN = phytools::getStates, type = "tips")
tip_states_matrix <- do.call(rbind, tip_states_list)
dim(tip_states_matrix)

# Compute posterior probabilties of states across simmaps
Primary_Photobiont_posterior_matrix <- t(apply(X = tip_states_matrix,
                                       MARGIN = 2,  
                                       FUN = function(x) { prop.table(table(factor(x, levels = Primary_Photobiont_states))) } ))
dim(Primary_Photobiont_posterior_matrix)

# Retrieve ML states for each species  
Primary_Photobiont_full <- apply(X = Primary_Photobiont_posterior_matrix,
                         MARGIN =  1,
                         FUN = function(x) names(which.max(x)))  

Lichen_data_df$Primary_Photobiont_full <- Primary_Photobiont_full

table(Lichen_data_df$Primary_Photobiont_full)
table(is.na(Lichen_data_df$Primary_Photobiont_full))

# Save objects
saveRDS(object = Primary_Photobiont_simmaps, file = "./outputs/Bromeliaceae/Primary_Photobiont_simmaps.rds")
saveRDS(object = Primary_Photobiont_posterior_matrix, file = "./outputs/Bromeliaceae/Primary_Photobiont_posterior_matrix.rds")
saveRDS(object = Lichen_data_df, file = "./input_data/Kessous_2025/Lichen_data_df.rds")



##### 3/ Build densityMaps of trait evolution #####

# Load phylo and Ecological_data_df
Lichen_tree <- readRDS(file = "./input_data/Nelsen_2020/Lichen_tree.rds")
Lichen_data_df <- readRDS(file = "./input_data/Nelsen_2020/Lichen_data_df.rds")


# # Load phylo and Ecological_data_df
# Lichen_tree_pruned <- readRDS(file = "./input_data/Nelsen_2020/Lichen_tree_pruned.rds")
# Lichen_data_df_pruned <- readRDS(file = "./input_data/Nelsen_2020/Lichen_data_df_pruned.rds")

# # Utrametrized tree
# Lichen_tree_scaled <- phytools::force.ultrametric(tree = Lichen_tree)
# ape::is.ultrametric(Lichen_tree_scaled)
# range(Lichen_tree_scaled$edge.length)

### 3.1/ For Growth_form ####

## Prepare trait data

# Extract categorical trait data as a named vector
Growth_form_tip_data <- setNames(object = Lichen_data_df$Growth_form,
                                 nm = Lichen_tree$tip.label)
head(Growth_form_tip_data)
table(Growth_form_tip_data)
table(is.na(Lichen_data_df$Primary_Photobiont))

table(names(Growth_form_tip_data) == Lichen_tree$tip.label)


# Select a color scheme per states
colors_per_states = c("Macro" = "darkgreen", "Micro" = "darkolivegreen2")

## Map trait evolution as ML estimates
Growth_form_mapped_data <- prepare_trait_data(
  tip_data = Growth_form_tip_data,
  trait_data_type = "categorical",
  phylo = Lichen_tree,
  # phylo = Lichen_tree_scaled,
  seed = 1234,
  evolutionary_models = "ARD",
  nb_simulations = 100,
  plot_map = FALSE,
  verbose = TRUE)

## Plot densityMaps = posterior frequencies of states
plot_densityMaps_overlay(densityMaps = Growth_form_mapped_data$densityMaps,
                         colors_per_levels = colors_per_states)

## Save output
saveRDS(object = Growth_form_mapped_data, file = "./outputs/Lichens/Growth_form_mapped_data.rds")


### 3.2/ For Primary_Photobiont ####

## Prepare trait data

# Extract categorical trait data as a named vector
Primary_Photobiont_tip_data <- setNames(
   # object = Lichen_data_df$Primary_Photobiont,
   object = Lichen_data_df$Primary_Photobiont_full, # Use the full version without NA
   nm = Lichen_tree$tip.label)

head(Primary_Photobiont_tip_data)
table(Primary_Photobiont_tip_data)
table(names(Primary_Photobiont_tip_data) == Lichen_tree$tip.label)

# Select a color scheme per states
colors_per_states = c("Cyanobacteria" = "deepskyblue3", "NonLichenized" = "grey70", "Trebouxiophyceae" = "purple", "Trentepohliales" = "orange")

## Map trait evolution as ML estimates
Primary_Photobiont_mapped_data <- prepare_trait_data(
  tip_data = Primary_Photobiont_tip_data,
  trait_data_type = "categorical",
  phylo = Lichen_tree,
  seed = 1234,
  evolutionary_models = "ARD",
  nb_simulations = 100,
  plot_map = FALSE,
  verbose = TRUE)

## Plot densityMaps = posterior frequencies of states
plot_densityMaps_overlay(densityMaps = Primary_Photobiont_mapped_data$densityMaps,
                         colors_per_levels = colors_per_states)

## Save output
saveRDS(object = Primary_Photobiont_mapped_data, file = "./outputs/Lichens/Primary_Photobiont_mapped_data.rds")



##### 4/ Build BAMM objects #####

## Load ultrametric phylogeny
Lichen_tree <- readRDS(file = "./input_data/Nelsen_2020/Lichen_tree.rds")

# Load initial BAMM object with all posterior samples
Lichen_BAMM_object <- BAMMtools::getEventData(phy = Lichen_tree, 
                                              eventdata = "./input_data/Nelsen_2020/BAMM/bamm_medusa_med_families_65m/bamm_medusa_med_families_65m_event_data.txt",
                                              burnin = 0,
                                              verbose = TRUE)

hist(Lichen_BAMM_object$numberEvents)

# Build BAMM object for deepSTRAPP
Lichen_BAMM_object <- deepSTRAPP::build_BAMM_object(
  phylo = Lichen_tree, 
  eventdata = "./input_data/Nelsen_2020/BAMM/bamm_medusa_med_families_65m/bamm_medusa_med_families_65m_event_data.txt",
  burn_in = 0.25, nb_posterior_samples = 1000, expectedNumberOfShifts = 34,
  verbose = TRUE)

hist(Lichen_BAMM_object$numberEvents)

plot_BAMM_rates(BAMM_object = Lichen_BAMM_object)


## Save output
saveRDS(object = Lichen_BAMM_object, file = "./outputs/Lichens/Lichen_BAMM_object.rds")


##### 5/ Run deepSTRAPP on Growth_form data #####

# Load output mapped trait data
Growth_form_mapped_data <- readRDS(file = "./outputs/Lichens/Growth_form_mapped_data.rds")
# Load BAMM_object summarizing 1000 posterior samples of BAMM
Lichen_BAMM_object <- readRDS(file = "./outputs/Lichens/Lichen_BAMM_object.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(Lichen_tree)[,2])
root_age # 248 Mya

## Set for time steps of 10 My from 0 to 150 Mya
# nb_time_steps <- 10
time_step_duration <- 10
time_range <- c(0, 150)


## Run deepSTRAPP on net diversification rates
Lichens_deepSTRAPP_Growth_form <- run_deepSTRAPP_over_time(
  densityMaps = Growth_form_mapped_data$densityMaps,
  nb_simulations = 100,
  trait_data_type = "categorical",
  BAMM_object = Lichen_BAMM_object,
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
str(Lichens_deepSTRAPP_Growth_form, max.level = 1)

# Display test summary
Lichens_deepSTRAPP_Growth_form$pvalues_summary_df

## Visualize test results

# Select a color scheme per states
colors_per_states = c("Macro" = "darkgreen", "Micro" = "darkolivegreen2")

# Plot p-values of Mann-Whitney-Wilcoxon tests across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form,
  alpha = 0.10)

# Plot evolution of mean rates through time
plot_rates_through_time(deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form,
                        colors_per_levels = colors_per_states,
                        plot_CI = TRUE)

# Plot rates vs. trait values across branches for time step = 10 My
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form,
  colors_per_levels = colors_per_states,
  focal_time = 10)

# Plot histogram of Mann-Whitney-Wilcoxon test stats for time step = 10 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form,
  focal_time = 10)

# Plot histograms of Mann-Whitney-Wilcoxon test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form,
  display_plots = TRUE)

## Save output
saveRDS(object = Lichens_deepSTRAPP_Growth_form, file = "./outputs/Lichens/Lichens_deepSTRAPP_Growth_form.rds")

## Load output
Lichens_deepSTRAPP_Growth_form <- readRDS(file = "./outputs/Lichens/Lichens_deepSTRAPP_Growth_form.rds")
    


##### 6/ Run deepSTRAPP on Growth_form data - Extented Time Range: 0 - 200 Mya #####

# Load output mapped trait data
Growth_form_mapped_data <- readRDS(file = "./outputs/Lichens/Growth_form_mapped_data.rds")
# Load BAMM_object summarizing 1000 posterior samples of BAMM
Lichen_BAMM_object <- readRDS(file = "./outputs/Lichens/Lichen_BAMM_object.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(Lichen_tree)[,2])
root_age # 248 Mya

## Set for time steps of 10 My from 0 to 200 Mya
# nb_time_steps <- 10
time_step_duration <- 10
time_range <- c(0, 200)


## Run deepSTRAPP on net diversification rates
Lichens_deepSTRAPP_Growth_form_0_200 <- run_deepSTRAPP_over_time(
  densityMaps = Growth_form_mapped_data$densityMaps,
  nb_simulations = 100,
  trait_data_type = "categorical",
  BAMM_object = Lichen_BAMM_object,
  seed = 1234,
  # nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  time_range = time_range,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  # alpha = 0.10,
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Lichens_deepSTRAPP_Growth_form_0_200, max.level = 1)

# Display test summary
Lichens_deepSTRAPP_Growth_form_0_200$pvalues_summary_df

## Visualize test results

# Select a color scheme per states
colors_per_states = c("Macro" = "darkgreen", "Micro" = "darkolivegreen2")

# Plot p-values of Mann-Whitney-Wilcoxon tests across all time-steps
Lichens_Growth_form_0_200_pvalues_plot <- plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_0_200,
  alpha = 0.10)

pdf(file = "./outputs/Lichens/Lichens_Growth_form_0_200_pvalues_plot.pdf", height = 6, width = 8)
print(Lichens_Growth_form_0_200_pvalues_plot)
dev.off()

# Plot evolution of mean rates through time
Lichens_Growth_form_0_200_RTT_plot <- plot_rates_through_time(
   deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_0_200,
   colors_per_levels = colors_per_states,
   plot_CI = TRUE)$rates_TT_ggplot

# Move legend
Lichens_Growth_form_0_200_RTT_plot <- Lichens_Growth_form_0_200_RTT_plot +
  ggplot2::theme(legend.position = "inside",
                 legend.position.inside = c(0.85, 0.25))

# Export without geological epochs
pdf(file = "./outputs/Lichens/Lichens_Growth_form_0_200_RTT_plot.pdf", height = 6, width = 8)
print(Lichens_Growth_form_0_200_RTT_plot)
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
Lichens_Growth_form_0_200_RTT_plot_with_epochs <- Lichens_Growth_form_0_200_RTT_plot
Lichens_Growth_form_0_200_RTT_plot_with_epochs$layers <- c(rect_layer, Lichens_Growth_form_0_200_RTT_plot_with_epochs$layers)

Lichens_Growth_form_0_200_pvalues_plot_with_epochs <- Lichens_Growth_form_0_200_pvalues_plot
Lichens_Growth_form_0_200_pvalues_plot_with_epochs$layers <- c(rect_layer, Lichens_Growth_form_0_200_pvalues_plot_with_epochs$layers)

# # Apply the discrete fill colors if your stored plot doesn't have a fill scale yet
# my_plot <- my_plot + 
#   scale_fill_manual(values = c("Cretaceous" = "gray90", 
#                                "Paleogene" = "gray75", 
#                                "Quaternary-Neogene" = "gray90"))

# Display the modified plot
print(Lichens_Growth_form_0_200_RTT_plot_with_epochs)

# Export with geological epochs
pdf(file = "./outputs/Lichens/Lichens_Growth_form_0_200_RTT_plot_with_epochs.pdf", height = 6, width = 8)
print(Lichens_Growth_form_0_200_RTT_plot_with_epochs)
dev.off()

pdf(file = "./outputs/Lichens/Lichens_Growth_form_0_200_pvalues_plot_with_epochs.pdf", height = 6, width = 8)
print(Lichens_Growth_form_0_200_pvalues_plot_with_epochs)
dev.off()


# Plot rates vs. trait values across branches for time step = 0 My
ggplot1 <- plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_0_200,
  colors_per_levels = colors_per_states,
  focal_time = 0)$rates_vs_trait_ggplot

# Plot rates vs. trait values across branches for time step = 60 My
ggplot2 <- plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_0_200,
  colors_per_levels = colors_per_states,
  focal_time = 60)$rates_vs_trait_ggplot

# Plot histogram of Mann-Whitney-Wilcoxon test stats for time step = 0 My
ggplot3 <- plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_0_200,
  focal_time = 0)

# Plot histogram of Mann-Whitney-Wilcoxon test stats for time step = 60 My
ggplot4 <- plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_0_200,
  focal_time = 60)

# Plot all plots on one faceted plot
pdf(file = "./outputs/Lichens/Lichens_Growth_form_output_facetted_plots.pdf", height = 12, width = 14)
cowplot::plot_grid(plotlist = list(ggplot1, ggplot2, ggplot3, ggplot4))
dev.off()



# Retrieve Q10% to adjust plot to the alpha used in RTT for T = 0 My
quantile(x = Lichens_deepSTRAPP_Growth_form_0_200$STRAPP_results_over_time[[1]]$perm_data_df$abs_delta_U, probs = 0.10)
quantile(x = Lichens_deepSTRAPP_Growth_form_0_200$STRAPP_results_over_time[[1]]$perm_data_df$abs_delta_U, probs = Lichens_deepSTRAPP_Growth_form_0_200$STRAPP_results_over_time[[1]]$p_value)

# Retrieve Q10% to adjust plot to the alpha used in RTT for T = 60 My
quantile(x = Lichens_deepSTRAPP_Growth_form_0_200$STRAPP_results_over_time[[7]]$perm_data_df$abs_delta_U, probs = 0.10)
quantile(x = Lichens_deepSTRAPP_Growth_form_0_200$STRAPP_results_over_time[[7]]$perm_data_df$abs_delta_U, probs = Lichens_deepSTRAPP_Growth_form_0_200$STRAPP_results_over_time[[7]]$p_value)


# Plot histograms of Mann-Whitney-Wilcoxon test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_0_200,
  display_plots = TRUE)

## Save output
saveRDS(object = Lichens_deepSTRAPP_Growth_form_0_200, file = "./outputs/Lichens/Lichens_deepSTRAPP_Growth_form_0_200.rds")


##### 7/ Run deepSTRAPP on Growth_form data - One-tailed test #####

# # Load faster script
# source("../deepSTRAPP/R/update_rates_and_regimes_for_focal_time.R")

# Load output mapped trait data
Growth_form_mapped_data <- readRDS(file = "./outputs/Lichens/Growth_form_mapped_data.rds")
# Load BAMM_object summarizing 1000 posterior samples of BAMM
Lichen_BAMM_object <- readRDS(file = "./outputs/Lichens/Lichen_BAMM_object.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(ape::as.phylo(Lichen_BAMM_object))[,2])
root_age # 248 Mya

## Set for time steps of 10 My from 0 to 150 Mya
# nb_time_steps <- 10
time_step_duration <- 10
time_range <- c(0, 150)


## Run deepSTRAPP on net diversification rates
Lichens_deepSTRAPP_Growth_form_one_tailed <- run_deepSTRAPP_over_time(
  densityMaps = Growth_form_mapped_data$densityMaps,
  nb_simulations = 100,
  trait_data_type = "categorical",
  BAMM_object = Lichen_BAMM_object,
  seed = 1234,
  # nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  time_range = time_range,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  # Run one-tailed test
  two_tailed = FALSE,
  one_tailed_hypothesis = c("Micro > Macro"),
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Lichens_deepSTRAPP_Growth_form_one_tailed, max.level = 1)

# Display test summary
Lichens_deepSTRAPP_Growth_form_one_tailed$pvalues_summary_df

## Visualize test results

# Select a color scheme per states
colors_per_states = c("Macro" = "darkgreen", "Micro" = "darkolivegreen2")

# Plot p-values of Mann-Whitney-Wilcoxon tests across all time-steps
pvalues_plot_Growth_form_one_tailed <- plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed,
  alpha = 0.05)

pvalues_plot_Growth_form_one_tailed <- pvalues_plot_Growth_form_one_tailed +
  ggplot2::ggtitle(label = "STRAPP tests\nDifferences in net diversification rates through time\nOne-tailed test\n")

print(pvalues_plot_Growth_form_one_tailed)

# Plot evolution of mean rates through time
RTT_plot_Growth_form_one_tailed <- plot_rates_through_time(
   deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed,
   colors_per_levels = colors_per_states,
   plot_CI = TRUE)$rates_TT_ggplot

RTT_plot_Growth_form_one_tailed <- RTT_plot_Growth_form_one_tailed +
  ggplot2::ggtitle(label = "Net diversification rates per trait states through time\nOne-tailed test\n")

print(RTT_plot_Growth_form_one_tailed)

# Plot rates vs. trait values across branches for time step = 10 My
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed,
  colors_per_levels = colors_per_states,
  focal_time = 10)

# Plot histogram of Mann-Whitney-Wilcoxon test stats for time step = 10 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed,
  focal_time = 10)

# Plot histograms of Mann-Whitney-Wilcoxon test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed,
  display_plots = TRUE)

## Save output
saveRDS(object = Lichens_deepSTRAPP_Growth_form_one_tailed, file = "./outputs/Lichens/Lichens_deepSTRAPP_Growth_form_one_tailed.rds")


##### 8/ Run deepSTRAPP on Growth_form data - One-tailed test - Extended time range #####

# # Load faster script
# source("../deepSTRAPP/R/update_rates_and_regimes_for_focal_time.R")

# Load output mapped trait data
Growth_form_mapped_data <- readRDS(file = "./outputs/Lichens/Growth_form_mapped_data.rds")
# Load BAMM_object summarizing 1000 posterior samples of BAMM
Lichen_BAMM_object <- readRDS(file = "./outputs/Lichens/Lichen_BAMM_object.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(ape::as.phylo(Lichen_BAMM_object))[,2])
root_age # 248 Mya

## Set for time steps of 10 My from 0 to 240 Mya
# nb_time_steps <- 10
time_step_duration <- 10
time_range <- c(0, 240)


## Run deepSTRAPP on net diversification rates
Lichens_deepSTRAPP_Growth_form_one_tailed_0_240 <- run_deepSTRAPP_over_time(
  densityMaps = Growth_form_mapped_data$densityMaps,
  nb_simulations = 100,
  trait_data_type = "categorical",
  BAMM_object = Lichen_BAMM_object,
  seed = 1234,
  # nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  time_range = time_range,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  # Run one-tailed test
  two_tailed = FALSE,
  one_tailed_hypothesis = c("Micro > Macro"),
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Lichens_deepSTRAPP_Growth_form_one_tailed_0_240, max.level = 1)

# Display test summary
Lichens_deepSTRAPP_Growth_form_one_tailed_0_240$pvalues_summary_df

## Visualize test results

# Select a color scheme per states
colors_per_states = c("Macro" = "darkgreen", "Micro" = "darkolivegreen2")

# Plot p-values of Mann-Whitney-Wilcoxon tests across all time-steps
pvalues_plot_Growth_form_one_tailed_0_240 <- plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed_0_240,
  alpha = 0.05)

pvalues_plot_Growth_form_one_tailed_0_240 <- pvalues_plot_Growth_form_one_tailed_0_240 +
  ggtitle(label = "STRAPP tests\nDifferences in net diversification rates through time\nOne-tailed test\n")

print(pvalues_plot_Growth_form_one_tailed_0_240)

# Plot evolution of mean rates through time
RTT_plot_Growth_form_one_tailed_0_240 <- plot_rates_through_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed_0_240,
  colors_per_levels = colors_per_states,
  plot_CI = TRUE)$rates_TT_ggplot

RTT_plot_Growth_form_one_tailed_0_240 <- RTT_plot_Growth_form_one_tailed_0_240 +
  ggtitle(label = "Net diversification rates per trait states through time\nOne-tailed test\n")

print(RTT_plot_Growth_form_one_tailed_0_240)

# Plot rates vs. trait values across branches for time step = 10 My
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed_0_240,
  colors_per_levels = colors_per_states,
  focal_time = 10)

# Plot histogram of Mann-Whitney-Wilcoxon test stats for time step = 10 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed_0_240,
  focal_time = 10)

# Plot histograms of Mann-Whitney-Wilcoxon test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Growth_form_one_tailed_0_240,
  display_plots = TRUE)

## Save output
saveRDS(object = Lichens_deepSTRAPP_Growth_form_one_tailed_0_240, file = "./outputs/Lichens/Lichens_deepSTRAPP_Growth_form_one_tailed_0_240.rds")


##### 9/ Run deepSTRAPP on Primary_Photobiont data #####

# Load output mapped trait data
Primary_Photobiont_mapped_data <- readRDS(file = "./outputs/Lichens/Primary_Photobiont_mapped_data.rds")
# Load BAMM_object summarizing 1000 posterior samples of BAMM
Lichen_BAMM_object <- readRDS(file = "./outputs/Lichens/Lichen_BAMM_object.rds")

## Identify root age
root_age <- max(phytools::nodeHeights(Lichen_tree)[,2])
root_age # 248 Mya

## Set for time steps of 10 My from 0 to 150 Mya
# nb_time_steps <- 10
time_step_duration <- 10
time_range <- c(0, 150)


## Run deepSTRAPP on net diversification rates
Lichens_deepSTRAPP_Primary_Photobiont <- run_deepSTRAPP_over_time(
  densityMaps = Primary_Photobiont_mapped_data$densityMaps,
  nb_simulations = 100,
  trait_data_type = "categorical",
  BAMM_object = Lichen_BAMM_object,
  seed = 1234,
  # nb_time_steps = nb_time_steps,
  time_step_duration = time_step_duration,
  time_range = time_range,
  uncertainty_strategy = "rates_only",
  rate_type = "net_diversification",
  posthoc_pairwise_tests = TRUE, # Run pairwise tests too
  return_perm_data = TRUE,
  extract_trait_data_melted_df = TRUE,
  extract_diversification_data_melted_df = TRUE,
  return_STRAPP_results = TRUE,
  return_updated_Maps = TRUE,
  return_updated_BAMM_objects = TRUE,
  verbose = TRUE,
  verbose_extended = TRUE)

## Explore output
str(Lichens_deepSTRAPP_Primary_Photobiont, max.level = 1)

# Display test summary for overall tests
Lichens_deepSTRAPP_Primary_Photobiont$pvalues_summary_df
# Display test summary for posthoc pairwise Dunn's tests
Lichens_deepSTRAPP_Primary_Photobiont$pvalues_summary_df_for_posthoc_pairwise_tests

## Visualize test results

# Select a color scheme per states
colors_per_states = c("Cyanobacteria" = "deepskyblue3", "NonLichenized" = "grey70", "Trebouxiophyceae" = "purple", "Trentepohliales" = "orange")

# Plot p-values of overall Kruskal-Wallis tests across all time-steps
plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Primary_Photobiont,
  alpha = 0.05)

# Plot p-values of post hoc pairwise Dunn's tests between pairs of tests across all time-steps
pvalues_plot <-plot_STRAPP_pvalues_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Primary_Photobiont,
  alpha = 0.05,
  plot_posthoc_tests = TRUE)

pvalues_plot <- pvalues_plot + ggplot2::theme(legend.key.width = unit(0.5, units = "cm"),
                                              legend.text = element_text(size = 10),
                                              legend.title = element_text(size = 12),
                                              legend.position.inside = c(0.20, 0.6))

print(pvalues_plot)

# Plot evolution of mean rates through time
plot_rates_through_time(deepSTRAPP_outputs = Lichens_deepSTRAPP_Primary_Photobiont,
                        colors_per_levels = colors_per_states,
                        plot_CI = TRUE)

# Plot rates vs. trait values across branches for time step = 10 My
plot_rates_vs_trait_data_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Primary_Photobiont,
  colors_per_levels = colors_per_states,
  focal_time = 10)

# Plot histogram of overall Kruskal-Wallis test stats for time step = 10 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Primary_Photobiont,
  focal_time = 10,
  plot_posthoc_tests = FALSE)

# Plot histogram of post hoc pairwise Dunn's tests stats for time step = 10 My
plot_histogram_STRAPP_test_for_focal_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Primary_Photobiont,
  focal_time = 10,
  plot_posthoc_tests = TRUE)

# Plot histograms of overall Kruskal-Wallis test results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Primary_Photobiont,
  display_plots = TRUE,
  plot_posthoc_tests = FALSE)

# Plot histograms of post hoc pairwise Dunn's tests results (One plot per time-step)
plot_histograms_STRAPP_tests_over_time(
  deepSTRAPP_outputs = Lichens_deepSTRAPP_Primary_Photobiont,
  display_plots = TRUE,
  plot_posthoc_tests = TRUE)


## Save output
saveRDS(object = Lichens_deepSTRAPP_Primary_Photobiont, file = "./outputs/Lichens/Lichens_deepSTRAPP_Primary_Photobiont.rds")


