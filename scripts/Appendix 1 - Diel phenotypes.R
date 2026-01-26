# Clear the workspace environment to ensure a clean start
rm(list = ls())

# ==============================================================================
# Script: Appendix 1 - Diel Phenotypes Classification
# ==============================================================================
# Purpose: This script classifies species into diel phenotypes (activity
#          patterns across the 24-hour cycle) using a Bayesian model-based
#          hypothesis framework proposed by Gerber et al. (2024).
#
# Method:  The analysis uses the "General" hypothesis set from Diel.Niche
#          package, which compares seven alternative models representing
#          different diel activity patterns:
#          - Unimodal: diurnal (day-active) or nocturnal (night-active)
#          - Bimodal: crepuscular (twilight-active)
#          - Trimodal: cathemeral (active throughout the diel cycle)
#
# Output:  The script produces posterior probability estimates for each
#          species' activity across three diel periods (twilight, day, night)
#          and visualizes these probabilities.
# ==============================================================================

# 1. Loading Required Packages ------------------------------------------------
# load("RData/Appendix 1 - Diel phenotypes.RData")
library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# Specialized packages for diel niche analysis
library(Diel.Niche)
library(bayesplot) # Visualization of Bayesian model outputs

# ==============================================================================
# 2. Loading and Formatting the Data ------------------------
# ==============================================================================

# Load data
detection_data <- fread("data/dataNicheOverlap.csv") %>%
    mutate(
        # Extract two-letter species codes from species names
        sp = substr(sp, 1, 2),

        # Convert date-time strings to POSIXct
        date_time = as.POSIXct(
            date_time,
            format = "%Y-%m-%d %H:%M:%S",
            tz = "UTC"
        ),

        # Assign study site coordinates (Doñana National Park, Spain)
        # These coordinates are required for calculating solar positions
        lat = 36.966,
        lon = -6.467
    )

head(detection_data)

# ==============================================================================
# 3. Classify Detections into Diel Time Periods -------------------------------
# ==============================================================================

# Create a bin list defining the time periods for diel classification
# This uses default parameters from the Diel.Niche package
diel_bins <- make.diel.bin.list(
    plot.bins = FALSE
)

# Classify each detection into diel time periods (twilight, day, night)
# based on solar position at the time and location of detection
#
# The bin.diel.times() function calculates solar elevation angles and
# assigns each observation to:
#   - Twilight: civil twilight periods (sun 0-6° below horizon)
#   - Day: sun above horizon
#   - Night: sun more than 6° below horizon
detections_classified <- bin.diel.times(
    data = detection_data,
    datetime.column = "date_time",
    lat.column = "lat",
    lon.column = "lon",
    bin.type.list = diel_bins
)

# ==============================================================================
# 4. Aggregate Detection Counts by Species and Diel Period --------------------
# ==============================================================================

# Count the number of detections for each species in each diel period
# This creates a long-format table with species, diel period, and counts
detection_counts_long <- detections_classified %>%
    group_by(sp, dielBin) %>%
    summarise(count = n())

# Reshape data to wide format for Bayesian analysis
# Each row = one species, columns = detection counts in each diel period
detection_matrix <- detection_counts_long %>%
    pivot_wider(
        names_from = dielBin,
        values_from = count,
        values_fill = 0
    ) %>%
    column_to_rownames(var = "sp") %>% # Set species as row names
    select(twilight, day, night)

# ==============================================================================
# 5. Bayesian Model Comparison: Fit Alternative Diel Phenotype Hypotheses -----
# ==============================================================================
# The "General" hypothesis set compares seven alternative models:
#   H1: Diurnal (primarily day-active)
#   H2: Nocturnal (primarily night-active)
#   H3: Diurnal-Crepuscular (day + twilight active)
#   H4: Nocturnal-Crepuscular (night + twilight active)
#   H5: Twilight-only (crepuscular)
#   H6: Cathemeral (active across all periods)
#   H7: Bimodal (day-night avoiding twilight)

# Define a function to fit diel phenotype models for a single species
# This function:
#   1. Fits all seven hypotheses to the detection data
#   2. Calculates Bayes factors (BF) to compare model support
#   3. Returns the most supported model and posterior probabilities
fit_diel_models <- function(y_numeric) {
    out <- diel.fit(
        t(as.matrix(y_numeric)), # Transpose count vector to required format
        hyp.set = hyp.sets("General"), # Use the "General" hypothesis set
        post.fit = FALSE, # Don't extract posterior samples yet
        n.chains = 3, # Run 3 MCMC chains for convergence check
        burnin = 5000, # Discard first 5000 iterations per chain
        n.mcmc = 10000, # Keep 10000 post-burnin iterations
        prints = FALSE # Suppress progress messages
    )

    # Return most supported model and Bayes factor table
    list(ms.model = out$ms.model, prob = out$bf.table)
}

# Apply the fitting function to each species (each row of detection_matrix)
# This performs model comparison separately for each species
model_fits <- apply(
    detection_matrix,
    1,
    fit_diel_models
)

# ==============================================================================
# 6. Extract Model Comparison Results ------------------------------------------
# ==============================================================================

# Extract the posterior probabilities for all seven hypotheses for each species
# Posterior probabilities sum to 1 for each species and indicate relative
# support for each hypothesis given the data
posterior_list <- sapply(
    model_fits,
    function(x) x[2]
)
model_posteriors <- matrix(
    unlist(
        lapply(
            posterior_list,
            "[", , "Posterior"
        )
    ),
    ncol = 7, # Seven hypotheses in "General" set
    byrow = TRUE
)
rownames(model_posteriors) <- rownames(detection_matrix) # Label rows with species codes
colnames(model_posteriors) <- rownames(posterior_list[[1]]) # Label columns with hypotheses

print(round(model_posteriors, digits = 2))

# Identify the most supported hypothesis for each species
# This is the hypothesis with the highest posterior probability
best_models <- unlist(
    lapply(model_fits, "[", 1)
)

print(best_models)

# Extract the posterior probability of the most supported hypothesis
# This indicates the strength of evidence for the best-supported model
best_model_probs <- unlist(
    lapply(
        lapply(model_fits, "[", 2),
        FUN = function(x) {
            max(x$prob[, 2])
        }
    )
)

print(round(best_model_probs, digits = 2))

# Combine most supported hypothesis and its probability into a data frame
best_models <- data.frame(
    best_models, best_model_probs
)

print(best_models)

# Create a data frame combining detection counts with most supported hypothesis
# This will be used for extracting posterior samples in the next step
species_counts_with_model <- data.frame(
    detection_matrix,
    hyp = best_models[, 1]
)

print(species_counts_with_model)

# ==============================================================================
# 7. Extract Posterior Samples for Activity Probabilities ---------------------
# ==============================================================================

# Define a function to extract posterior samples from the most supported
# hypothesis for each species. These samples represent the uncertainty in
# the estimated activity probabilities for each diel period.
extract_posterior_samples <- function(species_row) {
    out <- diel.fit(
        t(as.integer(species_row[-4])), # Detection counts (exclude hypothesis column)
        hyp.set = species_row[4], # Use the most supported hypothesis only
        post.fit = TRUE, # Extract full posterior samples
        prints = FALSE,
        n.chains = 3, # Three MCMC chains
        n.mcmc = 10000, # 10000 iterations per chain
        burnin = 1000 # Discard first 1000 iterations
    )

    # Return posterior samples and Gelman-Rubin convergence diagnostic
    list(post.samp = out$post.samp[[1]], gelman.diag = out$gelm.diag)
}

# Apply the posterior sampling function to each species
posterior_fits <- apply(
    species_counts_with_model,
    1,
    extract_posterior_samples
)

# Check MCMC convergence using Gelman-Rubin diagnostic
# Values close to 1.0 indicate good convergence
sapply(
    posterior_fits,
    function(x) x[2]
)

# Combine posterior samples from all chains for each species
# This creates a single matrix of posterior samples per species
posterior_samples <- lapply(
    sapply(
        posterior_fits,
        function(x) x[1]
    ),
    FUN = function(x) {
        do.call("rbind", x)
    }
)

# Calculate posterior quantiles (95% credible intervals and median)
# for activity probabilities in each diel period
posterior_quantiles <- lapply(
    posterior_samples,
    FUN = function(x) {
        apply(
            x,
            2,
            quantile,
            probs = c(
                0.025,
                0.5,
                0.975
            )
        )
    }
)


# ==============================================================================
# 8. Prepare Data for Visualization -------------------------------------------
# ==============================================================================

# Extract posterior medians for each species and diel period
# These represent the best point estimates of activity probabilities
posterior_medians <- matrix(
    unlist(
        lapply(
            posterior_quantiles,
            FUN = function(x) {
                x[2, ] # Extract the median (2nd row of quantiles: 0.025, 0.5, 0.975)
            }
        )
    ),
    ncol = 3,
    byrow = TRUE
)
rownames(posterior_medians) <- names(posterior_samples) # Label with species codes
colnames(posterior_medians) <- colnames(posterior_quantiles$Ge.post.samp) # Diel periods
print(round(posterior_medians, digits = 2))

# Prepare posterior samples for interval plots
plot_data <- posterior_samples
num_species <- length(plot_data)

# Rename columns of posterior samples to be more descriptive
plot_data <- lapply(
    plot_data,
    FUN = function(x) {
        colnames(x) <- c("P(twilight)", "P(daytime)", "P(nighttime)")
        x
    }
)

# Convert posterior samples to interval data format for ggplot
# This extracts median (prob = 0.5) and 95% credible intervals (prob_outer = 0.95)
plot_data <- do.call(
    "rbind",
    lapply(
        plot_data,
        FUN = function(x) {
            mcmc_intervals_data(
                x,
                prob = 0.5,
                prob_outer = 0.95
            )
        }
    )
)
plot_data$Species <- rep(
    rownames(detection_counts_long),
    each = 3,
    length.out = nrow(plot_data)
)

# Convert numeric species codes to full species names
plot_data <- plot_data %>%
    mutate(Species = case_when(
        Species == 1 ~ "Genet",
        Species == 2 ~ "Mongoose",
        Species == 3 ~ "Hare",
        Species == 4 ~ "Lynx",
        Species == 5 ~ "Badger",
        Species == 6 ~ "Rabbit",
        Species == 7 ~ "Fox"
    ))

# ==============================================================================
# 9. Create Visualization of Diel Activity Patterns ---------------------------
# ==============================================================================

# Add random vertical jitter to species points to prevent overlap
# Each species gets a small random offset, applied to all three diel periods
vertical_jitter <- rep(
    rnorm(nrow(plot_data) / 3, 0, 0.07),
    each = 3
)
point_positions <- position_nudge(
    y = vertical_jitter,
    x = rep(0, nrow(plot_data))
)

# Create a stacked bar chart showing the relative activity probabilities
# for each species across daytime, twilight, and nighttime periods.

plot_data$Species <- factor(
    plot_data$Species,
    levels = c(
        "Lynx", "Genet", "Mongoose",
        "Badger", "Fox", "Hare", "Rabbit"
    )
)

activity_bar_plot <- ggplot(
    plot_data,
    aes(x = m, y = Species, fill = parameter)
) +
    geom_bar(
        stat = "identity",
        position = "fill"
    ) +
    scale_fill_manual(values = c(
        "P(daytime)" = "#F2D027",
        "P(twilight)" = "#238C82",
        "P(nighttime)" = "#430E59"
    )) +
    labs(
        x = "Activity probabilities",
        y = "Species",
        fill = "Diel Period"
    ) +
    theme(
        text = element_text(
            size = 20,
            family = "sans"
        ),
        panel.background = element_rect(
            fill = "white",
            colour = "black"
        ),
        legend.key = element_rect(
            fill = "white",
            colour = "white"
        )
    )


print(activity_bar_plot)

# Save the plot
ggsave(
    filename = "results/Appendix 1 - Diel phenotypes.png",
    plot = activity_bar_plot,
    width = 8,
    height = 6,
    dpi = 300
)
