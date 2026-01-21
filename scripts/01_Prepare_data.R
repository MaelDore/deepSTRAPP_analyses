
##### Script 01: Prepare data for analyses #####

####################################
#       Author: Maël Doré          #
#  Contact: mael.dore@gmail.com    #
####################################

### Goals

# Curate taxonomic information from species list to prepare for analyses
# Extract species-level information

###

### Inputs

# Time-calibrated phylogeny
# BAMM object for diversification and morphological trait evolution
# Binary table of species occurrences in subbasins
# Environmental data at subbasin level

###

### Sources 

## Phylogeny and morphological data

# Cerezer, F.O., Dambros, C.S., Coelho, M.T.P. et al. Accelerated body size evolution in upland environments is correlated with recent speciation in South American freshwater fishes.
# Nature Communications. 14, 6070 (2023). https://doi.org/10.1038/s41467-023-41812-7

# Zenodo records: https://zenodo.org/records/8301082

## Biogeographic (sub)regions

# Morrone, 2022

###

### Outputs

# ContMap for morphological traits
# Environmental (abiotic) trait data at species level
# Biogeographic binary table of species occurrence in Subregions

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

# remotes::install_github(repo = "MaelDore/deepSTRAPP",
#                         # Time-consuming, but needed if you want to have access to the vignettes/tutorials
#                         build_vignettes = FALSE)

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

saveRDS(object = morpho_tip_data_df, file = "./outputs/Morpho/morpho_tip_data_df.rds")

# Load trait data
morpho_tip_data_df <- readRDS(file = "./outputs/Morpho/morpho_tip_data_df.rds")

### 1.3/ Load the time-calibrated phylogeny ####

Fish_tree <- read.tree(file = "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Phylogeny/NFF_subset_2638spp.tre")

plot(Fish_tree)

table(Fish_tree$tip.label %in% morpho_tip_data_df$Taxa)

Fish_tree <- force.ultrametric(Fish_tree)

### 1.4/ Load BAMM results for diversification ####

# Generate an bammdata object from BAMM output
BAMMdata_object_diversification_full <- BAMMtools::getEventData(
  phy = Fish_tree, 
  eventdata = "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Speciation_rates/BAMM/event_data.txt",
  burnin = 0.2) 

head(BAMMdata_object_diversification_full$eventData)

# Save BAMMdata_object_diversification
saveRDS(object = BAMMdata_object_diversification_full, file = "./outputs/BAMM/BAMMdata_object_diversification_full.rds")

## Subsample to 1000 posterior samples to save time

set.seed(seed = 4231)
posterior_samples_ID <- sample(x = 1:length(BAMMdata_object_diversification_full$eventData), size = 1000, replace = F)

BAMMdata_object_diversification <- BAMMtools::subsetEventData(ephy = BAMMdata_object_diversification_full,
                                                              index = posterior_samples_ID)

# Save BAMMdata_object_diversification
saveRDS(object = BAMMdata_object_diversification, file = "./outputs/BAMM/BAMMdata_object_diversification.rds")


### 1.5/ Load BAMM results for trait evolution ####

## BEL = Body Elongation

# Generate a BAMMdata object from BAMM output
BAMMdata_BEL_full <- BAMMtools::getEventData(
  phy = Fish_tree, type = "trait",
  eventdata = paste0(morpho_folder, "BEL/BAMM_outputs/event_data.txt"),
  burnin = 0.2) 

length(BAMMdata_BEL_full$eventData)

# Save BAMMdata_BEL
saveRDS(object = BAMMdata_BEL_full, file = "./outputs/BAMM/BAMMdata_BEL_full.rds")

## Subsample to 1000 posterior samples to save time

set.seed(seed = 4231)
posterior_samples_ID <- sample(x = 1:length(BAMMdata_BEL_full$eventData), size = 1000, replace = F)

BAMMdata_BEL <- BAMMtools::subsetEventData(ephy = BAMMdata_BEL_full,
                                           index = posterior_samples_ID)

# Save BAMMdata_BEL
saveRDS(object = BAMMdata_BEL, file = "./outputs/BAMM/BAMMdata_BEL.rds")


## MBL = Body Elongation

# Generate a BAMMdata object from BAMM output
BAMMdata_MBL_full <- BAMMtools::getEventData(
  phy = Fish_tree, type = "trait",
  eventdata = paste0(morpho_folder, "MBL/BAMM_outputs/event_data.txt"),
  burnin = 0.2) 

head(BAMMdata_MBL_full$eventData)

# Save BAMMdata_MBL
saveRDS(object = BAMMdata_MBL_full, file = "./outputs/BAMM/BAMMdata_MBL_full.rds")

## Subsample to 1000 posterior samples to save time

set.seed(seed = 4231)
posterior_samples_ID <- sample(x = 1:length(BAMMdata_object_MBL_full$eventData), size = 1000, replace = F)

BAMMdata_MBL <- BAMMtools::subsetEventData(ephy = BAMMdata_MBL_full,
                                           index = posterior_samples_ID)

# Save BAMMdata_MBL
saveRDS(object = BAMMdata_MBL, file = "./outputs/BAMM/BAMMdata_MBL.rds")


##### 2/ Convert BAMM trait rates data into contMap ####

# Can I convert this to a contMap with mean trait rate values along branches ???

## Make a function to convert BAMM to contMap
# For a given posterior sample
# For all samples
# For mean values across all samples

BAMMdata_MBL_test <- BAMMtools::getEventData(
  phy = Fish_tree, type = "trait",
  eventdata = paste0(morpho_folder, "MBL/BAMM_outputs/event_data.txt"),
  burnin = 0.999) 

head(BAMMdata_MBL_test$eventData)

contMap_MBL_test <- BAMMdata_object_to_contMap(BAMMdata_object = BAMMdata_MBL_test, sample_ID = 1, type = "trait", display_plot = FALSE)

plot_contMap(contMap = contMap_MBL_test[[1]])
plot(contMap_MBL_test[[1]])

plot_contMap(contMap = contMap_MBL_test)

### 2.1/ Functions to convert BAMM trait rates data into contMap ####

BAMMdata_object_to_contMap <- function (BAMMdata_object, sample_ID, type,
                                        res = 100, color_scale = NULL,
                                        display_plot = TRUE, PDF_file_path = NULL,
                                        verbose = TRUE, ...)
{
  ### Check input validity
  {
    ## BAMMdata_object
    # Check that object is of class "bammdata"
    if (!any(class(BAMMdata_object) == "bammdata"))
    {
      stop(paste0("'BAMMdata_object' must be of class 'bammdata'. See BAMMtools::getEventData() to obtain this type of objects."))
    }
    
    # Check that object has $eventData slot with list of data.frame
    if (is.null(BAMMdata_object$eventData))
    {
      stop(paste0("'BAMMdata_object' must have a '$eventData' element. See BAMMtools::getEventData() to obtain this type of objects."))
    }
    
    ## sample_ID
    # Check that sample_ID is within the range of $eventData
    nb_samples <- length(BAMMdata_object$eventData)
    if (!(sample_ID %in% 1:nb_samples))
    {
      stop(paste0("'sample_ID' is outside the range of posterior samples recorded.\n",
                  "'sample_ID' = ",sample_ID,". Number of posterior samples = ",nb_samples ,".\n",
                  "Please select a 'sample_ID' that falls within the range of available posterior samples."))
    }
    
    ## type
    # Check that type is either "diversification" or "trait"
    if (!(type %in% c("diversification", "trait")))
    {
      stop(paste0("'type' must be either 'diversification' or 'trait' according to what was modeled by BAMM."))
    }
    
    ## res
    # Check that res is a positive integer
    if (!((round(res) == res) & (res > 0)))
    {
      stop(paste0("'res' describes the number of time segments used to discretize time. It must be a positive interger."))
    }
    
    ## color_scale
    # Check whether all colors are valid
    if (!is.null(color_scale))
    {
      if (!all(is_color(color_scale)))
      {
        invalid_colors <- color_scale[!is_color(color_scale)]
        stop(paste0("Some color names in 'color_scale' are not valid.\n",
                    "Invalid: ", paste(invalid_colors, collapse = ", "), "."))
      }
    }
    
    ## PDF_file_path
    # If provided, PDF_file_path must end with ".pdf"
    if (!is.null(PDF_file_path))
    {
      if (length(grep(pattern = "\\.pdf$", x = PDF_file_path)) != 1)
      {
        stop("'PDF_file_path' must end with '.pdf'")
      }
    }
  }
  
  # Extract lims agrument if provided
  add_args <- list(...)
  if ("lims" %in% names(add_args))
  {
    lims <- add_args[names(add_args) %in% "lims"]
  } else {
    lims <- NULL
  }
  
  # Extract eventData table
  eventData_i <- BAMMdata_object$eventData[[sample_ID]]
  # Extract eventBranchSegs matrix
  eventBranchSegs_i <- BAMMdata_object$eventBranchSegs[[sample_ID]]
  colnames(eventBranchSegs_i) <- c("tipward_node_ID", "rootward_height", "tipward_height", "regime_ID")
  
  # Extract tree from BAMMdata
  tree_elements <- c("edge", "edge.length", "Nnode", "tip.label", "node.label")
  tree <- BAMMdata_object[tree_elements[tree_elements %in% names(BAMMdata_object)]]
  class(tree) <- "phylo"
  
  # Set steps to discretize the 'continuous' gradient
  h <- max(phytools::nodeHeights(tree)) # Get root age
  step_times <- c(0:(res - 1)/(res - 1) * h, h + h/(res - 1)) # Compute step boundaries
  H <- phytools::nodeHeights(tree) # Get height (root - age) of edge limits
  
  # Get mean height of each segment (average of boundary heights)
  mean_step_times <- step_times[-length(step_times)] + (step_times[2] / 2)
  
  # Compute rates for each regime across mean step times
  speciation_rates_per_regimes_mat <- matrix(data = NA, nrow = nrow(eventData_i), ncol = length(mean_step_times))
  if (type == "diversification")
  {
    net_div_rates_per_regimes_mat <- extinction_rates_per_regimes_mat <- speciation_rates_per_regimes_mat
  }
  
  # Loop per regimes (j)
  for (j in 1:nrow(eventData_i))
  {
    # j <- 98
    
    # Extract regime parameters
    lambda_0_j <- eventData_i$lam1[j]
    lambda_alpha_j <- eventData_i$lam2[j]
    mu_0_j <- eventData_i$mu1[j]
    mu_alpha_j <- eventData_i$mu2[j]
    
    # Set mean relative time since regime initiation
    relative_mean_time_steps_j <- mean_step_times - eventData_i$time[j]
    # Should remove rates computed before shift, but to ensure that segments that
    # overlap with the shit timing get rate estimes, we keep them
    # relative_mean_time_steps_j[relative_mean_time_steps_j < 0] <- NA
    
    # Compute speciation rates based on regime parameters
    if (lambda_alpha_j <= 0) # If alpha <= 0 (decrease): lambda_t = lambda_0 * exp(alpha*t)
    {
      speciation_rates_j <- lambda_0_j * exp(lambda_alpha_j*relative_mean_time_steps_j)
    } else { # If alpha > 0 (increase): lambda_t = lambda_0 * (2 - exp(-alpha*t))
      speciation_rates_j <- lambda_0_j * (2 - exp(-lambda_alpha_j*relative_mean_time_steps_j))
    }
    
    # Plot to check results
    # plot(y = speciation_rates_j, x = relative_mean_time_steps_j)
    
    if (type == "diversification")
    {
      # Compute extinction rates based on regime parameters
      if (mu_alpha_j <= 0) # If alpha <= 0 (decrease): mu_t = mu_0 * exp(alpha*t)
      {
        extinction_rates_j <- mu_0_j * exp(mu_alpha_j*relative_mean_time_steps_j)
      } else { # If alpha > 0 (increase): mu_t = mu_0 * (2 - exp(-alpha*t))
        extinction_rates_j <- mu_0_j * (2 - exp(-mu_alpha_j*relative_mean_time_steps_j))
      }
      
      # Compute net div rates
      net_div_rates_j <- speciation_rates_j - extinction_rates_j
    }
    
    
    # Store results
    speciation_rates_per_regimes_mat[j, ] <- speciation_rates_j
    if (type == "diversification")
    {
      net_div_rates_per_regimes_mat[j, ] <- net_div_rates_j
      extinction_rates_per_regimes_mat[j, ] <- extinction_rates_j
    }
  }
  
  # Set default colors based on rainbow
  cols <- rainbow(1001, start = 0, end = 0.7)
  names(cols) <- 0:1000
  
  ## Set transformation vector from raw rates to 0/1000 scale
  
  ## For speciation
  
  # Extract all values to set limits
  all_speciation_rates <- as.vector(speciation_rates_per_regimes_mat)
  all_speciation_rates <- all_speciation_rates[!is.na(all_speciation_rates)]
  
  # Transform values from scale 0 to 1000
  if (is.null(lims))
  {
    speciation_lims <- c(min(all_speciation_rates), max(all_speciation_rates))
  } else {
    speciation_lims <- lims
  }
  speciation_trans <- 0:1000/1000 * (speciation_lims[2] - speciation_lims[1]) + speciation_lims[1]
  names(speciation_trans) <- 0:1000
  
  if (type == "diversification")
  {
    ## For net_div
    
    # Extract all values to set limits
    all_net_div_rates <- as.vector(net_div_rates_per_regimes_mat)
    all_net_div_rates <- all_net_div_rates[!is.na(all_net_div_rates)]
    
    # Transform values from scale 0 to 1000
    if (is.null(lims))
    {
      net_div_lims <- c(min(all_net_div_rates), max(all_net_div_rates))
    } else {
      net_div_lims <- lims
    }
    net_div_trans <- 0:1000/1000 * (net_div_lims[2] - net_div_lims[1]) + net_div_lims[1]
    names(net_div_trans) <- 0:1000
    
    ## For extinction
    
    # Extract all values to set limits
    all_extinction_rates <- as.vector(extinction_rates_per_regimes_mat)
    all_extinction_rates <- all_extinction_rates[!is.na(all_extinction_rates)]
    
    # Transform values from scale 0 to 1000
    if (is.null(lims))
    {
      extinction_lims <- c(min(all_extinction_rates), max(all_extinction_rates))
    } else {
      extinction_lims <- lims
    }
    extinction_trans <- 0:1000/1000 * (extinction_lims[2] - extinction_lims[1]) + extinction_lims[1]
    names(extinction_trans) <- 0:1000
  }

  ## Initiate "$maps"
  tree$maps <- vector(mode = "list", length = nrow(tree$edge)) # One list of trait values per edge
  tree_speciation <- tree
  if (type == "diversification")
  {
    tree_net_div <- tree_extinction <- tree
  }
  
  ## Fill "$maps"
  
  # Loop per edge (k)
  for (k in 1:nrow(tree$edge)) 
  {
    # k <- 349
    
    ## Get tipward_node_ID
    tipward_node_ID_k <- tree$edge[k, 2]
    
    ## Identify mean time steps found along the edge
    mean_time_steps_ID_k <- intersect(which(step_times > H[k, 1]), which(step_times < H[k, 2]))
    mean_time_steps_ID_k <- c(mean_time_steps_ID_k[1] - 1, mean_time_steps_ID_k)
    # mean_step_times[mean_time_steps_ID_k] # Mean step times may fall outside of the edge range if end segments are cut after/befoer their mean time
    
    ## Fix issue with edge that are within a time segment
    if (any(is.na(mean_time_steps_ID_k)))
    {
      mean_edge_time_k <- mean(H[k, 2], H[k, 1])
      mean_time_steps_ID_k <- which.min(abs(mean_step_times - mean_edge_time_k))
    }
    
    ## Identify regime ID found along segments of the edge
    
    # Extract events/regimes on edge k
    eventBranchSegs_edge_k <- as.data.frame(eventBranchSegs_i[eventBranchSegs_i[, 1] == tipward_node_ID_k, , drop = FALSE])
    # Get time boundaries of segments along the edge
    segments_time_boundaries_k <- cbind(c(H[k, 1], step_times[intersect(which(step_times > H[k, 1]), which(step_times < H[k, 2]))]),
                                        c(step_times[intersect(which(step_times >  H[k, 1]), which(step_times < H[k, 2]))], H[k, 2]))
    segments_relative_time_boundaries_k <- segments_time_boundaries_k - H[k, 1]
    # Compute time spent per each segment in each regime
    time_segment_per_regimes_k <- matrix(nrow = nrow(segments_time_boundaries_k), ncol = nrow(eventBranchSegs_edge_k))
    for (i in 1:nrow(eventBranchSegs_edge_k))
    {
      # i <- 1
      
      # Compute time overlap between segment and regime i on the edge
      overlap_start <- sapply(X = segments_time_boundaries_k[, 1], FUN = function(x) { max(x, eventBranchSegs_edge_k$rootward_height[i]) })
      overlap_end <- sapply(X = segments_time_boundaries_k[, 2], FUN = function(x) { min(x, eventBranchSegs_edge_k$tipward_height[i]) })
      overlap_i <- overlap_end - overlap_start
      overlap_i[overlap_i < 0] <- 0
      
      time_segment_per_regimes_k[, i] <- overlap_i
    }
    # Define regime identity of segments based on longest time of the segment under each regime
    segments_regimes_ID_k <- eventBranchSegs_edge_k$regime_ID[apply(X = time_segment_per_regimes_k, MARGIN = 1, FUN = which.max)]
    
    ## Get values for each segment base on regime ID and mean step times
    speciation_segments_values_k <- c()
    if (type == "diversification")
    {
      net_div_segments_values_k <- extinction_segments_values_k <- c()
    }
    
    for (i in seq_along(segments_regimes_ID_k))
    {
      # For speciation
      speciation_segments_values_k[i] <- speciation_rates_per_regimes_mat[segments_regimes_ID_k[i], mean_time_steps_ID_k[i]]
      
      if (type == "diversification")
      {
        # For net_div
        net_div_segments_values_k[i] <- net_div_rates_per_regimes_mat[segments_regimes_ID_k[i], mean_time_steps_ID_k[i]]
        # For net_div
        extinction_segments_values_k[i] <- extinction_rates_per_regimes_mat[segments_regimes_ID_k[i], mean_time_steps_ID_k[i]]
        
      }
    }
    
    # Compute length of each segment
    segments_length_k <- segments_relative_time_boundaries_k[, 2] - segments_relative_time_boundaries_k[, 1] 
    # Save them as values in the '$maps'
    tree_net_div$maps[[k]] <- tree_speciation$maps[[k]] <- tree_extinction$maps[[k]] <- segments_length_k 
      
    # Transform values into the 0 to 1000 scale
    speciation_scale_k <- sapply(speciation_segments_values_k, getState, trans = speciation_trans)
    if (type == "diversification")
    {
      net_div_scale_k <- sapply(net_div_segments_values_k, getState, trans = net_div_trans)
      extinction_scale_k <- sapply(extinction_segments_values_k, getState, trans = extinction_trans)
    }
     # Save scaled values as names (scale from 0 to 1000)
    names(tree_speciation$maps[[k]]) <- speciation_scale_k
    if (type == "diversification")
    {
      names(tree_net_div$maps[[k]]) <- net_div_scale_k 
      names(tree_extinction$maps[[k]]) <- extinction_scale_k
    }

    ## Print progress
    if (verbose & (k %% 1000 == 0))
    {
      if (type == "diversification")
      {
        cat(paste0(Sys.time(), " - Rates mapped for edge n\u00B0", k, "/", nrow(tree$edge),"\n"))
      } else {
        cat(paste0(Sys.time(), " - Trait values mapped for edge n\u00B0", k, "/", nrow(tree$edge),"\n"))
      }
    }
  }
  
  ## Build '$mapped.edge' from the '$maps'
  
  # For speciation
  tree_speciation$mapped.edge <- makeMappedEdge(edge = tree_speciation$edge, maps = tree_speciation$maps)
  tree_speciation$mapped.edge <- tree_speciation$mapped.edge[, order(as.numeric(colnames(tree_speciation$mapped.edge))), drop = FALSE]
  
  if (type == "diversification")
  {
    # For net_div
    tree_net_div$mapped.edge <- makeMappedEdge(edge = tree_net_div$edge, maps = tree_net_div$maps)
    tree_net_div$mapped.edge <- tree_net_div$mapped.edge[, order(as.numeric(colnames(tree_net_div$mapped.edge))), drop = FALSE]

    # For extinction
    tree_extinction$mapped.edge <- makeMappedEdge(edge = tree_extinction$edge, maps = tree_extinction$maps)
    tree_extinction$mapped.edge <- tree_extinction$mapped.edge[, order(as.numeric(colnames(tree_extinction$mapped.edge))), drop = FALSE]
    
  }
  
  ## Set classes and attributes of '$tree'
  
  # For speciation
  class(tree_speciation) <- c("simmap", setdiff(class(tree_speciation), "simmap"))
  attr(tree_speciation, "map.order") <- "right-to-left"
  
  if (type == "diversification")
  {
    # For net_div
    class(tree_net_div) <- c("simmap", setdiff(class(tree_net_div), "simmap"))
    attr(tree_net_div, "map.order") <- "right-to-left"

    # For extinction
    class(tree_extinction) <- c("simmap", setdiff(class(tree_extinction), "simmap"))
    attr(tree_extinction, "map.order") <- "right-to-left"
    
  }
  
  ## Build contMaps
  
  if (type == "diversification")
  {
    # For net_div
    contMap_net_div <- list(tree = tree_net_div, cols = cols, lims = net_div_lims)
    class(contMap_net_div) <- "contMap"
    
    # For speciation
    contMap_speciation <- list(tree = tree_speciation, cols = cols, lims = speciation_lims)
    class(contMap_speciation) <- "contMap"
    
    # For extinction
    contMap_extinction <- list(tree = tree_extinction, cols = cols, lims = extinction_lims)
    class(contMap_extinction) <- "contMap"
  } else {
    
    # For trait data
    contMap_trait <- list(tree = tree_speciation, cols = cols, lims = speciation_lims)
    class(contMap_trait) <- "contMap"
  }

  ## Update colors + display plot
  
  if (type == "diversification")
  {
    # For net_div
    updated_contMap_net_div <- deepSTRAPP::plot_contMap(contMap = contMap_net_div, 
                                                        color_scale = color_scale,
                                                        display_plot = display_plot,
                                                        ...)
    if (display_plot) { title(main = paste0("\n\nNet diversification rates - Sample ", sample_ID)) }
    
    # For speciation
    updated_contMap_speciation <- deepSTRAPP::plot_contMap(contMap = contMap_speciation, 
                                                           color_scale = color_scale,
                                                           display_plot = display_plot,
                                                           ...)
    if (display_plot) { title(main = paste0("\n\nSpeciation rates - Sample ", sample_ID)) }
    
    # For extinction
    updated_contMap_extinction <- deepSTRAPP::plot_contMap(contMap = contMap_extinction, 
                                                           color_scale = color_scale,
                                                           display_plot = display_plot,
                                                           ...)
    if (display_plot) { title(main = paste0("\n\nExtinction rates - Sample ", sample_ID)) }
    
  } else {
    
    # For trait data (will also print PDF)
    updated_contMap_trait <- deepSTRAPP::plot_contMap(contMap = contMap_trait, 
                                                      color_scale = color_scale,
                                                      display_plot = display_plot,
                                                      PDF_file_path,
                                                      ...)
    if (display_plot) { title(main = paste0("\n\nTrait evolution - Sample ", sample_ID)) }
  }
  
  ## Print PDF (for diversification rates)
  
  if (!is.null(PDF_file_path) & type == "diversification")
  {
    # Adjust width and height according to phylo
    nb_tips <- length(contMap$tree$tip.label)
    height <- min(nb_tips/60*10, 200) # Maximum PDF size = 200 inches
    width <- height*8/10
    
    # Open PDF
    grDevices::pdf(file = file.path(PDF_file_path),
                   width = width, height = height)
    
    ## Plot the contMaps
    phytools::plot.contMap(x = updated_contMap_net_div,
                           plot = TRUE,
                           ...)
    title(main = "\n\nNet diversification rates")
    
    phytools::plot.contMap(x = updated_contMap_speciation,
                           plot = TRUE,
                           ...)
    title(main = "\n\nSpeciation rates")
    
    phytools::plot.contMap(x = updated_contMap_extinction,
                           plot = TRUE,
                           ...)
    title(main = "\n\nExtinction rates")
    
    # Close PDF
    grDevices::dev.off()
  }

  ## Build output
  if (type == "diversification")
  {
    output_list <- list(contMap_net_div = updated_contMap_net_div,
                        contMap_speciation = updated_contMap_speciation,
                        contMap_extinction = updated_contMap_extinction)
  } else {
    output_list <- updated_contMap_trait
  }

  # Return output
  return(invisible(output_list))
}


all_contMaps_MBL <- BAMMdata_object_to_all_contMaps(BAMMdata_object = BAMMdata_MBL_test, type = "trait")
all_contMaps_MBL_div <- BAMMdata_object_to_all_contMaps(BAMMdata_object = BAMMdata_MBL_test, type = "diversification")

## Function to extract trait/rates from all posterior samples
BAMMdata_object_to_all_contMaps <- function (BAMMdata_object, type,
                                             res = 100, color_scale = NULL,
                                             verbose = TRUE, ...)
{
  # Initiate list
  all_contMaps <- list()
  
  # Loop per posterior samples
  for ( i in seq_along(BAMMdata_object$eventData))
  {
    all_contMaps[[i]] <- BAMMdata_object_to_contMap(BAMMdata_object = BAMMdata_object, sample_ID = i, type = type,
                                                    res = res, color_scale = color_scale,
                                                    display_plot = FALSE, verbose = FALSE,
                                                    ...)
    
    ## Print progress
    if (verbose & (i %% 10 == 0))
    {
      if (type == "diversification")
      {
        cat(paste0(Sys.time(), " - Rates mapped for posterior sample n\u00B0", i, "/", length(BAMMdata_object$eventData),"\n"))
      } else {
        cat(paste0(Sys.time(), " - Trait values mapped for posterior sample n\u00B0", i, "/", length(BAMMdata_object$eventData),"\n"))
      }
    }
  }
  
  # Name contMap
  names(all_contMaps) <- paste0("contMap_", seq_along(BAMMdata_object$eventData))
  
  # For diversification data: Reorganize output per type of rates
  if (type == "diversification")
  {
    all_net_div_contMaps <- lapply(X = all_contMaps, FUN = function (x) { x[["contMap_net_div"]]} )
    all_speciation_contMaps <- lapply(X = all_contMaps, FUN = function (x) { x[["contMap_speciation"]]} )
    all_extinction_contMaps <- lapply(X = all_contMaps, FUN = function (x) { x[["contMap_extinction"]]} )
    
    names(all_net_div_contMaps) <- paste0("contMap_net_div_", seq_along(BAMMdata_object$eventData))
    names(all_speciation_contMaps) <- paste0("contMap_speciation_", seq_along(BAMMdata_object$eventData))
    names(all_extinction_contMaps) <- paste0("contMap_extinction_", seq_along(BAMMdata_object$eventData))
    
    all_contMaps <- list(net_div_contMaps = all_net_div_contMaps,
                         speciation_contMaps = all_speciation_contMaps,
                         extinction_contMaps = all_extinction_contMaps)
  }
  
  # Return list(s) of contMaps
  return(all_contMaps)
}

all_contMaps <- all_contMaps_MBL

## Function to average a list of contMaps
compute_mean_contMap <- function (all_contMaps)
{
  ## Extract raw values
  all_raw_maps <- list()
  # Loop per posterior sample / contMap
  for (i in seq_along(all_contMaps))
  {
    # i <- 1
    
    # Extract '$maps'
    scaled_maps_i <- all_contMaps[[i]]$tree$maps
    # Retrieve scale
    range <- all_contMaps[[i]]$lims
    trans_i <- 0:1000/1000 * (range[2] - range[1]) + range[1]
    # Convert to raw values
    raw_maps_i <- lapply(X = scaled_maps_i, FUN = function (x) { trans_i[as.numeric(names(x))] })
    # Store
    all_scaled_maps[[i]] <- scaled_maps_i
    all_raw_maps[[i]] <- raw_maps_i
  }
  
  ## Compute mean values across posterior samples
  
  # Initiate mean maps
  mean_values <- list()
  # Loop per egde
  for (i in seq_along(all_raw_maps[[1]]))
  {
    # i <- 1
    
    # Extract all values across posterior samples for edge i 
    all_values_edge_i <- do.call(rbind, args = lapply(X = all_raw_maps, FUN = function (x) { x[[i]] }))
    # Compute mean values across samples
    mean_values_edge_i <- colMeans(all_values_edge_i)
    
    # Store as vector of values
    mean_values[[i]] <- mean_values_edge_i
  }
  
  ## Build new scale for mean values
  all_mean_values <- unlist(mean_values)
  all_mean_values_range <- range(all_mean_values)
  mean_trans <- 0:1000/1000 * (all_mean_values_range[2] - all_mean_values_range[1]) + all_mean_values_range[1]
  names(mean_trans) <- 0:1000
  
  ## Convert values to new scale
  mean_scaled_values <- lapply(X = mean_values, FUN = function (x) { sapply(X = x, FUN = getState, trans = mean_trans) } )
  
  ## Produce '$maps'
  
  # Initiate maps with segment lengths
  mean_maps <- all_contMaps[[1]]$tree$maps
  # Loop per edge to update names with scaled values
  for (i in seq_along(mean_maps))
  {
    mean_maps_i <- mean_maps[[i]]
    names(mean_maps_i) <- mean_scaled_values[[i]]
    mean_maps[[i]] <- mean_maps_i
  }
  
  ## Update '$mapped.edge'
  mean_mapped.edge <- makeMappedEdge(edge = all_contMaps[[1]]$tree$edge, maps = mean_maps)
  
  ## Produce mean '$tree'
  mean_tree <- all_contMaps[[1]]$tree
  mean_tree$maps <- mean_maps
  mean_tree$mapped.edge <- mean_mapped.edge
  
  class(mean_tree) <- c("simmap", setdiff(class(mean_tree), "simmap"))
  attr(mean_tree, "map.order") <- "right-to-left"
  
  # Produce contMap
  
  mean_contMap <- list(tree = mean_tree,
                       cols = all_contMaps[[1]]$cols,
                       lims = all_mean_values_range)
  class(mean_contMap) <- "contMap"
  
  # Return mean_contMap
  return(mean_contMap)
}


## Function to aggregate all posterior samples and plot mean rates/trait values
BAMMdata_object_to_mean_contMap <- function (BAMMdata_object, type,
                                             res = 100, color_scale = NULL,
                                             display_plot = TRUE, PDF_file_path = NULL,
                                             verbose = TRUE, ...)
{
  ## Convert all BAMM samples into contMap
  all_contMaps <- BAMMdata_object_to_all_contMaps(BAMMdata_object = BAMMdata_object, 
                                                  type = type,
                                                  res = res, color_scale = color_scale,
                                                  verbose = verbose)
  
  ## Compute mean_contMap across all samples
  if (type == "trait")
  {
    mean_contMap_trait <- compute_mean_contMap(all_contMaps)
    
    # Built output
    output_list <- mean_contMap_trait
    
  }
  
  if (type == "diversification")
  {
    # For net_div
    mean_contMap_net_div <- compute_mean_contMap(all_contMaps$net_div_contMaps)
    # For speciation
    mean_contMap_speciation <- compute_mean_contMap(all_contMaps$speciation_contMaps)
    # For extinction
    mean_contMap_extinction <- compute_mean_contMap(all_contMaps$extinction_contMaps)
    
    # Built output
    output_list <- list(mean_contMap_net_div = mean_contMap_net_div,
                        mean_contMap_speciation = mean_contMap_speciation,
                        mean_contMap_extinction = mean_contMap_extinction)
  }
  
  ## Plot if requested
  if (display_plot)
  {
    if (type == "diversification")
    {
      # For net_div
      deepSTRAPP::plot_contMap(contMap = mean_contMap_net_div, ...)
      title(main = paste0("\n\nMean net diversification rates"))
      
      # For speciation
      deepSTRAPP::plot_contMap(contMap = mean_contMap_speciation, ...)
      title(main = paste0("\n\nMean speciation rates"))
      
      # For extinction
      deepSTRAPP::plot_contMap(contMap = mean_contMap_extinction, ...)
      title(main = paste0("\n\nMean extinction rates"))
      
    } else {
      
      # For trait data (will also print PDF)
      deepSTRAPP::plot_contMap(contMap = mean_contMap_trait,
                               PDF_file_path, ...)
      title(main = paste0("\n\nMean trait evolution"))
    }
  }
  
  ## Print PDF (for diversification rates)
  
  if (!is.null(PDF_file_path) & type == "diversification")
  {
    # Adjust width and height according to phylo
    nb_tips <- length(mean_contMap_net_div$tree$tip.label)
    height <- min(nb_tips/60*10, 200) # Maximum PDF size = 200 inches
    width <- height*8/10
    
    # Open PDF
    grDevices::pdf(file = file.path(PDF_file_path),
                   width = width, height = height)
    
    ## Plot the contMaps
    phytools::plot.contMap(x = mean_contMap_net_div, plot = TRUE, ...)
    title(main = paste0("\n\nMean net diversification rates"))

    phytools::plot.contMap(x = mean_contMap_speciation, plot = TRUE, ...)
    title(main = paste0("\n\nMean speciation rates"))
    
    phytools::plot.contMap(x = mean_contMap_extinction, plot = TRUE, ...)
    title(main = paste0("\n\nMean extinction rates"))
    
    # Close PDF
    grDevices::dev.off()
  }
  
  ## Return output
  return(output_list)
}


## Define function to transform raw values to 0/1000 scale
# Internal function from phytools package. See phytools:::getState
getState <- function (x, trans) 
{
  if (x <= trans[1]) 
    state <- names(trans)[1]
  else if (x >= trans[length(trans)]) 
    state <- names(trans)[length(trans)]
  else {
    i <- 1
    while (x > trans[i]) {
      state <- names(trans)[i]
      i <- i + 1
    }
  }
  state
}

## Define function to convert list of mapping ('$maps') in matrix of state per edge ('$makeMappedEdge')
# # Internal function from phytools package. See phytools:::makeMappedEdge
makeMappedEdge <- function (edge, maps) 
{
  st <- sort(unique(unlist(sapply(maps, function(x) names(x)))))
  mapped.edge <- matrix(0, nrow(edge), length(st))
  rownames(mapped.edge) <- apply(edge, 1, function(x) paste(x, collapse = ","))
  colnames(mapped.edge) <- st
  for (i in 1:length(maps))
  {
    for (j in 1:length(maps[[i]]))
    {
      mapped.edge[i, names(maps[[i]])[j]] <- mapped.edge[i, names(maps[[i]])[j]] + maps[[i]][j]
    }
  }
  return(mapped.edge)
}


### 2.2 Convert MBL rates into contMap ####

# Load BAMM output for MBL
BAMMdata_MBL <- readRDS(file = "./outputs/BAMM/BAMMdata_MBL.rds")

# mean_contMap_MBL_rates <- BAMMdata_object_to_mean_contMap(BAMMdata_object = BAMMdata_MBL_test, type = "trait")
mean_contMap_MBL_rates <- BAMMdata_object_to_mean_contMap(BAMMdata_object = BAMMdata_MBL, type = "trait")

updated_contMap <- phytools::setMap(x = mean_contMap_MBL_rates, colors = c("grey", "red", "gold", "limegreen", "forestgreen", "dodgerblue", "darkblue"))
plot_contMap(contMap = updated_contMap, color_scale = c("grey", "red", "gold", "limegreen", "forestgreen", "dodgerblue", "darkblue"))

saveRDS(object = mean_contMap_MBL_rates, file = "./outputs/Morpho/mean_contMap_MBL_rates.rds")


##### 3/ Categorize trait data #####

# If trait data are categorized, need to rerun the trait evolution using ARE/ER/SYM models

hist(morpho_tip_data_df$BEL) # Log-normal
hist(log(morpho_tip_data_df$BEL))

hist(morpho_tip_data_df$MBL) # Normal

hist(morpho_tip_data_df$OGP) # Can be binarized!
table(morpho_tip_data_df$OGP) # 0 vs. more than 0
morpho_tip_data_df$OGP_binary <- morpho_tip_data_df$OGP > 0
table(morpho_tip_data_df$OGP_binary)

hist(morpho_tip_data_df$RES) # Normal

hist(morpho_tip_data_df$RML) # Log-normal
hist(log(morpho_tip_data_df$RML)) # Still some extreme values

morpho_tip_data_df$OGP_binary <- morpho_tip_data_df$OGP > 0
table(morpho_tip_data_df$OGP_binary)

saveRDS(object = morpho_tip_data_df, file = "./outputs/Morpho/morpho_tip_data_df.rds")


##### 4/ Extract abiotic variables at species level #####

library(sf)

Basin_shp <- read_sf(dsn = "./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Shapefiles/Hydroatlas/", layer = "BasinATLAS_v10_lev05")
Basin_shp <- st_make_valid(Basin_shp)
st_is_valid(Basin_shp)

plot(Basin_shp[, "MAIN_BAS"])

### 4.1/ Get binary table of subbasin occurrence ####

Basin_occurrence_binary_table_df <- read.csv("./input_data/Cerezer_2023_Neotropical_freshwater_fishes/Datasets/Occurrences/PresAbs_SubBasin_subset_2638spp.csv")

### 4.2/ Get subbasin level data ####

SubBasins_data_df <- read.csv("./input_data/Cerezer_2023_Neotropical_freshwater_fishes/All_data_combined/SubBasins_estimates.csv")

# bio1 = Mean temperature
# bio12 = Annual total Precipitation
# aet = Evapotranspiration (mm / syr)
# Elevation

### 4.3/ Aggregate at species level

# Keep only shared basins
SubBasins_data_df_subset <- SubBasins_data_df %>%
  filter(HYBAS_ID %in% Basin_occurrence_binary_table_df$HYBAS_ID) %>% 
  select(HYBAS_ID, bio1, bio12, aet, Elevation)
Basin_occurrence_binary_table_df_subset <- Basin_occurrence_binary_table_df %>% 
  filter(HYBAS_ID %in% SubBasins_data_df$HYBAS_ID)

# Ensure basins are similarly ordered
Basin_occurrence_binary_table_df_subset$HYBAS_ID == SubBasins_data_df_subset$HYBAS_ID

# Remove basin ID for computation
row.names(Basin_occurrence_binary_table_df_subset) <- Basin_occurrence_binary_table_df_subset$HYBAS_ID
Basin_occurrence_binary_table_df_subset <- Basin_occurrence_binary_table_df_subset %>% 
  select(-HYBAS_ID)
row.names(SubBasins_data_df_subset) <- SubBasins_data_df_subset$HYBAS_ID
SubBasins_data_df_subset<- SubBasins_data_df_subset %>% 
  select(-HYBAS_ID)

# Use matrix product to compute the sum of variable values where each species is found
abiotic_data_total <- t(as.matrix(Basin_occurrence_binary_table_df_subset)) %*% as.matrix(SubBasins_data_df_subset)
# Divide by the number of basins to get average per species
basin_prevalence <- colSums(Basin_occurrence_binary_table_df_subset)
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

abiotic_data_df

# Save abiotic data per species
saveRDS(object = abiotic_data_df, file = "./outputs/Abiotic/abiotic_data_df.rds")

## Explore data distribution

abiotic_data_df <- readRDS(file = "./outputs/Abiotic/abiotic_data_df.rds")

hist(abiotic_data_df$Elevation) # Log-normal
hist(log(abiotic_data_df$Elevation))

hist(abiotic_data_df$Temp) # Left-skewed

hist(abiotic_data_df$Prec) # Normal

hist(abiotic_data_df$Evapotranspiration) # Left-skewed


##### 5/ Extract regions based on subbasin occurrences #####

# Use definition of regions in the Neotropics as in Morrone, 2022? (Terrestrial and floristic...)
# https://doi.org/10.1590/0001-3765202220211167 
# Morrone, 2022
#  - Antillean, Brazilian and Chacoan subregions
#  - Mexican and South American transition zones

library(sf)

Morrone_shp <- read_sf(dsn = "./input_data/Morrone_2022/NeotropicMap_Geo/", layer = "NeotropicMap_Geo")
Morrone_shp <- st_make_valid(Morrone_shp)
st_is_valid(Morrone_shp)

plot(Morrone_shp[, "Subregion"])

Morrone_Subregions_shp <- Morrone_shp
# Merge N/A and South American Transition Zone
Morrone_Subregions_shp$Subregion[Morrone_Subregions_shp$Subregion == "N/A"] <- "South American Transition Zone"
# Merge Antillean into Brazilian
Morrone_Subregions_shp$Subregion[Morrone_Subregions_shp$Subregion == "Antillean"] <- "Brazilian"
# Melt to Subregion level
Morrone_Subregions_shp <- Morrone_Subregions_shp[, "Subregion"] %>% 
  group_by(Subregion) %>% 
  summarize(geometry = st_union(geometry)) %>% 
  ungroup()

plot(Morrone_Subregions_shp)

# Save Morrone subregions
saveRDS(Morrone_Subregions_shp, file = "./input_data/Morrone_2022/Morrone_Subregions_shp.rds")

## Make a union between Morrone subregions and Subbasins shp files
sf::sf_use_s2(FALSE)
Morrone_cross_basins_sf <- st_intersection(x = Morrone_Subregions_shp, y = Basin_shp[, "HYBAS_ID"])
Morrone_cross_basins_sf <- st_make_valid(Morrone_cross_basins_sf)
sf::sf_use_s2(TRUE)

plot(Morrone_cross_basins_sf)

# Attribute to each basin the Subregion with the most overlap in terms of area
sf::sf_use_s2(FALSE)
Morrone_cross_basins_sf <- Morrone_cross_basins_sf %>% 
  mutate(area = st_area(Morrone_cross_basins_sf)) %>%
  group_by(HYBAS_ID) %>%
  arrange(HYBAS_ID, area) %>%
  mutate(rank = rank(area))
sf::sf_use_s2(TRUE)

## Simpler => Use the list of intersecting areas
Morrone_cross_subbasins_list <- st_intersects(x = Morrone_Subregions_shp, y = Basin_shp[, "HYBAS_ID"])
names(Morrone_cross_subbasins_list) <- levels(as.factor(Morrone_Subregions_shp$Subregion))

# Check if all basins have a match
all_basins_with_match <- unlist(Morrone_cross_subbasins_list)
table(1:nrow(Basin_shp) %in% all_basins_with_match)

# Convert into binary table
Subregions_subbasin_binary_table <- matrix(data = 0,
                                           nrow = nrow(Basin_shp), ncol = 3)
for (i in 1:3)
{
  j <- c(1, 2, 4)[i]
  Subbasin_matched_i <- Morrone_cross_subbasins_list[[j]]
  Subregions_subbasin_binary_table[Subbasin_matched_i, i] <- 1
}
Subregions_subbasin_binary_table <- as.data.frame(Subregions_subbasin_binary_table)
names(Subregions_subbasin_binary_table) <- names(Morrone_cross_subbasins_list)[c(1, 2, 4)]
row.names(Subregions_subbasin_binary_table) <- Basin_shp$HYBAS_ID

test <- Basin_shp[, "HYBAS_ID"] %>% 
  mutate(check = HYBAS_ID == "6050040460")
plot(test[, "check"])

## Then assign Subregion to species based on Sub-basins

# Filter to match the subbasins
Basin_occurrence_binary_table_df_subset <- Basin_occurrence_binary_table_df
row.names(Basin_occurrence_binary_table_df_subset) <- Basin_occurrence_binary_table_df_subset$HYBAS_ID
Basin_occurrence_binary_table_df_subset <- Basin_occurrence_binary_table_df_subset %>% 
  select(-HYBAS_ID)

Subregions_subbasin_binary_table_subset <- Subregions_subbasin_binary_table %>% 
  filter(row.names(Subregions_subbasin_binary_table) %in% row.names(Basin_occurrence_binary_table_df_subset))


Subregions_occurrence_table <- t(as.matrix(Basin_occurrence_binary_table_df_subset)) %*% as.matrix(Subregions_subbasin_binary_table_subset) 
Subregions_occurrence_binary_table <- (Subregions_occurrence_table > 0) * 1

saveRDS(object = Subregions_occurrence_binary_table, file = "./outputs/Biogeo/Subregions_occurrence_binary_table.rds")

# Results from Boris's paper at lower level than Region? Does not seem to be available. Check his paper results to see if he explains if the clustering stopped at level 2

Leroy_output <- readRDS(file = "./input_data/Leroy_2019/metrics_per_site_wide.RDS")


#### 6/ Extract additional ecological traits ####

# See https://www.nature.com/articles/s41597-025-04674-w


