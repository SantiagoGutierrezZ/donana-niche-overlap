# ==============================================================================
# Script: 02_temporal_analysis.R
# ==============================================================================
# Description: Analysis of temporal activity patterns and overlap between
# predator and prey species in Donana National Park.
# This script calculates average anchored times to account for seasonal variation
# in daylight, fits circular kernel density estimates, and estimates the
# coefficient of overlap (Dhat) with bootstrap confidence intervals.
# ==============================================================================

# 1. Loading packages -----------------------------------------------------
library(dplyr)
library(activity)
library(overlap)

# ==============================================================================
# 2. Loading data -----------------------------------------------
# ==============================================================================
# load("RData/02_temporal_analysis.RData")
data <- read.csv(
    "data/dataNicheOverlap.csv",
    header = TRUE
)

head(data)
# ==============================================================================
# 3. Calculate average anchored times -------------------------------------
# Since the daylight regime in the study area varies across the year (ranging
# from 9.6 hours in winter to 14.7 hours in summer), we use the `solartime`
# function to calculate the average sunrise and sunset times. These are used
# as anchor points to rescale the temporal data to a common scale where
# sunrise and sunset distinct events are fixed.
# See details in Vázquez et al. (2019).
# ==============================================================================

# Convert date_time to POSIXct.
# Cameras were set to UTC time throughout the study period.
data$date_time <- as.POSIXct(
    data$date_time,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
)

data$anchored_time <- solartime(
    data$date_time,
    lat = 36.9669, # Doñana Latitude
    lon = -6.4671, # Doñana Longitude
    tz = 0 # Offset for UTC data
)$solar # The 'solar' column contains the anchored times

# ==============================================================================
# 4. Calculate circular kernel density estimates --------------------------
# ==============================================================================

vulpes <- fitact(
    data$anchored_time[data$sp == "Vulpes vulpes"],
    sample = "data"
)
meles <- fitact(
    data$anchored_time[data$sp == "Meles meles"],
    sample = "data"
)
genetta <- fitact(
    data$anchored_time[data$sp == "Genetta genetta"],
    sample = "data"
)
herpestes <- fitact(
    data$anchored_time[data$sp == "Herpestes ichneumon"],
    sample = "data"
)
lynx <- fitact(
    data$anchored_time[data$sp == "Lynx pardinus"],
    sample = "data"
)
oryctolagus <- fitact(
    data$anchored_time[data$sp == "Oryctolagus cuniculus"],
    sample = "data"
)
lepus <- fitact(
    data$anchored_time[data$sp == "Lepus granatensis"],
    sample = "data"
)

# ==============================================================================
# 5. Calculate temporal overlap between predators and prey.
# We calculate the coefficient of overlap (Delta). Use Dhat4 if both species
# have > 50 detections, otherwise use Dhat1 (Ridout & Linkie 2009).
# ==============================================================================

# Initialize an empty data frame to store results
overlap_results <- data.frame(
    pair = character(),
    overlap = numeric(),
    type = character(),
    bootCI = character(),
    result = character(),
    stringsAsFactors = FALSE
)

predators <- c(
    "Vulpes vulpes",
    "Meles meles",
    "Genetta genetta",
    "Herpestes ichneumon",
    "Lynx pardinus"
)

prey <- c(
    "Oryctolagus cuniculus",
    "Lepus granatensis"
)

# WARNING!!!: Number of bootstrap repetitions
n_reps <- 10000

for (i in seq_along(predators)) {
    for (j in seq_along(prey)) {
        pair_name <- paste(predators[i], prey[j], sep = "_")

        # Calculate number of detections per species
        n_pred <- sum(data$sp == predators[i])
        n_prey <- sum(data$sp == prey[j])

        # Determine estimator type based on sample size
        # (use Dhat4 if both > 50, otherwise Dhat1)
        type <- ifelse(n_pred > 50 & n_prey > 50, "Dhat4", "Dhat1")

        # Calculate overlap
        est <- overlapEst(
            data$anchored_time[data$sp == prey[j]],
            data$anchored_time[data$sp == predators[i]],
            type = type
        )

        # Bootstrap for Confidence Intervals
        bs <- bootstrap(
            data$anchored_time[data$sp == prey[j]],
            data$anchored_time[data$sp == predators[i]],
            n_reps,
            smooth = TRUE,
            type = type
        )

        # Calculate CI
        ci <- bootCIlogit(est, bs, conf = 0.99)

        lower_ci <- ci["basic", 1]
        upper_ci <- ci["basic", 2]

        # Format CI string
        # Accessing the "norm" row of the CI matrix
        bootCI_string <- paste0(
            round(lower_ci, 2),
            "-",
            round(upper_ci, 2)
        )
        # Format text for figure
        result_text <- paste0(
            type, ": ",
            round(as.numeric(est), 2), " ",
            "(", bootCI_string, ")"
        )

        # Add a new row to the results
        overlap_results <- rbind(overlap_results, data.frame(
            pair = pair_name,
            overlap = round(as.numeric(est), 2),
            type = type,
            bootCI = bootCI_string,
            result = result_text,
            stringsAsFactors = FALSE
        ))
    }
}

print(overlap_results)

# ==============================================================================
# 6. Plot activity patterns ----------------------------------------------------
# We define a custom function ("activity_pattern") to plot the kernel density
# estimates, axis ticks, and anchored sunrise/sunset times.
# All plots are centered on midnight.
# ==============================================================================

# Define the activity pattern plot function
activity_pattern <- function(model, add, lty) {
    plot(model,
        data = "none",
        xunit = "radians", yunit = "density",
        ylim = c(0, 0.8),
        centre = "night",
        ylab = "", xlab = "",
        tline = list(col = "black", lwd = 3, lty = lty),
        dline = list(col = "transparent"),
        cline = list(col = "transparent"),
        xaxis = list(xaxt = "n"),
        add = add
    )
}

# Define tick marks in radians for a nocturnal x-axis (centered on midnight)
hours_in_radians_night <- c(-3, -2, -1, 0, 1, 2, 3)
draw_night_xaxis <- function() {
    axis(
        1,
        at = hours_in_radians_night,
        labels = c("3", "4", "5", "0", "1", "2", "3"),
        tick = TRUE,
        las = 1,
        padj = 0.5
    )
}

# Define average sunrise and sunset anchored times.
# We use as.POSIXlt to avoid "invalid 'tz' value" error in get_suntimes.
# Raw solar times are calculated in hours.
solar_times_raw <- get_suntimes(
    as.POSIXlt(data$date_time, tz = "UTC"),
    lat = 36.9669,
    lon = -6.4671,
    tz = 0,
    offset = 0
)

# Convert solar times to radians (x * pi / 12)
data$sunrise <- solar_times_raw[, 1] * pi / 12
data$sunset <- solar_times_raw[, 2] * pi / 12

draw_average_sunrise <- function() {
    lines(
        rep(cmean(data$sunrise), 2),
        c(0, 2500),
        col = "red",
        lwd = 2
    )
}

draw_average_sunset <- function() {
    lines(
        rep(cmean(data$sunset), 2) - 2 * pi,
        c(0, 2500),
        col = "red",
        lwd = 2
    )
}

draw_text <- function(text) {
    text(
        x = 0,
        y = 0.74,
        labels = text,
        cex = 4.5,
        font = 2,
        pos = 3
    )
}

# ==============================================================================
# 7. Generate Figures ---------------------------------------------------------
# ==============================================================================

png(
    "results/figure 3 - Temporal overlap for each predator-prey pair.png",
    width = 2000,
    height = 2800
)
par(
    family = "serif",
    cex.lab = 5,
    cex.axis = 5,
    mar = c(4.5, 6, 3, 3), mgp = c(6, 2, 0)
)
layout(
    matrix(
        1:10,
        ncol = 2,
        nrow = 5
    )
)

# Vulpes vulpes vs Oryctolagus cuniculus
activity_pattern(
    model = vulpes,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = oryctolagus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[1, 5])

# Meles meles vs Oryctolagus cuniculus
activity_pattern(
    model = meles,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = oryctolagus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[3, 5])

# Herpestes ichneumon vs Oryctolagus cuniculus
activity_pattern(
    model = herpestes,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = oryctolagus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[7, 5])

# Genetta genetta vs Oryctolagus cuniculus
activity_pattern(
    model = genetta,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = oryctolagus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[5, 5])

# Lynx pardinus vs Oryctolagus cuniculus
activity_pattern(
    model = lynx,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = oryctolagus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[9, 5])

# Vulpes vulpes vs Lepus granatensis
activity_pattern(
    model = vulpes,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = lepus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[2, 5])

# Meles meles vs Lepus granatensis
activity_pattern(
    model = meles,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = lepus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[4, 5])

# Herpestes ichneumon vs Lepus granatensis
activity_pattern(
    model = herpestes,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = lepus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[8, 5])

# Genetta genetta vs Lepus granatensis
activity_pattern(
    model = genetta,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = lepus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[6, 5])

# Lynx pardinus vs Lepus granatensis
activity_pattern(
    model = lynx,
    add = FALSE,
    lty = 1
)
activity_pattern(
    model = lepus,
    add = TRUE,
    lty = 2
)
draw_night_xaxis()
draw_average_sunset()
draw_average_sunrise()
draw_text(overlap_results[10, 5])

layout(1)
dev.off()


save.image("RData/02_temporal_analysis.RData")
