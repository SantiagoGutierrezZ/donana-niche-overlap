# ==============================================================================
# Script: 01_spatial_analysis.R
# ==============================================================================
# Description: Analysis of spatial overlap between the mesocarnivore mammal
# community and leporid species in Doñana National Park.
# This script calculates the spatial overlap (Pianka's index) based on
# detections across different sites.
# ==============================================================================

# 1. Loading packages ----------------------------------------------------------
library(dplyr)
library(tidyr)

# ==============================================================================
# 2. Loading and Preparing Data ------------------------------------------------
# ==============================================================================

# load("RData/01_spatial_analysis.RData")

# Read data
data <- read.csv("data/dataNicheOverlap.csv", header = TRUE) %>%
    select(-date_time)

# Create frequency table (Site x Species matrix)
# Rows: Sites, Columns: Species, Values: Number of detections of each species
# at each site
detection_matrix <- data %>%
    count(site, sp) %>%
    pivot_wider(names_from = sp, values_from = n, values_fill = 0)

# Define carnivores and leporids species
carnivores <- c(
    "Vulpes vulpes", "Herpestes ichneumon", "Genetta genetta",
    "Lynx pardinus", "Meles meles"
)
leporids <- c("Lepus granatensis", "Oryctolagus cuniculus")

# Convert to data frame and set row names to site names
species_df <- as.data.frame(detection_matrix)
rownames(species_df) <- species_df$site
species_df <- species_df %>% select(-site) # Remove the 'site' column

# Convert counts to proportions by species (column-wise)
# This represents the proportion of a species' total activity that occurs
# at each site.
prop_matrix <- apply(species_df, 2, function(x) {
    if (sum(x) == 0) {
        return(x)
    } # Handle species with 0 detections to avoid division by zero
    return(x / sum(x))
})

# ==============================================================================
# 3. Define Pianka's Index Function --------------------------------------------
# ==============================================================================

# Function to calculate Pianka's Index for spatial overlap
# Input: Two vectors of proportions (p and q) representing spatial usage
pianka_index_filtered <- function(p, q) {
    # Filter to include only sites where at least one species is present
    # (Note: Mathematically, sites with 0 for both don't affect the sum,
    # but this can be useful for data handling)
    valid_idx <- (p > 0) | (q > 0)
    p <- p[valid_idx]
    q <- q[valid_idx]

    numerator <- sum(p * q)
    denominator <- sqrt(sum(p^2) * sum(q^2))

    # Avoid division by zero if denominator is 0 (e.g., one species absent)
    if (denominator == 0) {
        return(0)
    }

    return(numerator / denominator)
}

# ==============================================================================
# 4. Calculate Overlap ---------------------------------------------------------
# ==============================================================================

# Calculate the index for all Predator (Carnivore) - Prey (Leporid) pairs
results <- expand.grid(
    predator_sp = carnivores,
    prey_sp = leporids,
    stringsAsFactors = FALSE
) %>%
    rowwise() %>%
    mutate(
        pianka_index = round(pianka_index_filtered(
            prop_matrix[, predator_sp],
            prop_matrix[, prey_sp]
        ), 3)
    ) %>%
    ungroup()

# View results
print(results)

# Save results to CSV
write.csv(results, "results/Pianka_index_results.csv", row.names = FALSE)
# save RData
save.image("RData/01_spatial_analysis.RData")
