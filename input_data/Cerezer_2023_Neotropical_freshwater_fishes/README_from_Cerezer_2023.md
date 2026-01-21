# Accelerated body size evolution in upland environments is correlated with recent speciation in South American freshwater fishes
---

This repository contains scripts and documentation organized into distinct folders, each tailored to specific functions within the analysis. A brief outline of the folder structure and its contents is provided below. For detailed information on setting up and using these scripts, please refer to the README file.




## Speciation Rates Estimates (Speciation_rates)


In this folder, you'll find codes and data necessary for conducting Bayesian analysis of macroevolutionary mixtures (BAMM) and computing Diversification Rate (DR) statistics.

* **BAMM**: The R script *'BAMM_codes.R'* is provided alongside essential BAMM analysis outputs: *'chain_swap.txt'*, *'divcontrol.txt'*, *'event_data.txt'*, *'mcmc_out.txt'*, and *'myPriors.txt'*.
* **DR**: Utilize the R script *'DR_codes.R'* to execute the DR statistic calculation.
* **Geography Speciation**: The R script *'Geography_codes.R'* facilitates mapping recent speciation rates (tip rates) across geographical areas, specifically South American sub-basins, with averaged rates based on co-occurring species lists.

## Biotic Variables (Biotic_factors)

This folder offers codes and data for assessing species diversity and rates of morphological evolution.

* **Species Diversity**: Quantify species richness in each sub-basin using the R script *'Species_diversity.R'*.
* **Morphological Evolution**: Estimate trait evolution (e.g., body elongation, maximum body length, oral gape position, relative eye size, and relative maxillary length) using the R script *'TraitEvol_codes.R'*. Additionally, BAMM outputs (*'chain_swap.txt'*, *'divcontrol.txt', 'event_data.txt', 'mcmc_out.txt', and 'myPriors.txt'*) are provided.

## Abiotic Variables (Abiotic_factors)

This folder provides codes and data for extracting climate and habitat variables.

* **Climate**: Extract climatic variables within each sub-basin using the R script *'Climatic_codes.R'*.
* **Habitat**: Extract habitat variables within each sub-basin using the R script *'Habitat_codes.R'*.

## All Extracted Data Combined (All_data_combined)

Here, you'll find the combined estimates for all sub-basins in the file named *'SubBasins_estimates.csv'*.

## Additional Source of Data (Datasets)

This folder contains supplementary files used in the study, including morphological traits, species occurrences, phylogenetic data, and shapefiles.

* **Morphology**: A dataset containing five morphological traits for 2,638 species, sourced from the fishmorph database (*'Morphology_subset_2638spp.csv'*).
* **Occurrences**: A presence/absence matrix for 2,638 species across South American sub-basins (*'PresAbs_SubBasin_subset_2638.csv'*).
* **Phylogeny**: The phylogenetic relationship among 2,638 freshwater fish species (*'NFF_subset_2638spp.tre'*).
* **Shapefiles**: This folder comprises polygons delineating South American boundaries (*'Continent_border'*), hydroatlas database (*'Hydroatlas'*), and sub-basins (*SubBasins*).

## Statistical Analyses (Statistics)

Execute various statistical analyses such as multiple linear regression, hierarchical partitioning, and variance partitioning using the R script *'Statistics_codes.R'*.

## Sensitivity Analyses (Sensitivity_test)

Perform sensitivity analyses, including the removal of biological outliers and correlation assessments between alternative metrics of speciation rates, using the R script *'Sensitivity_codes.R'*.
