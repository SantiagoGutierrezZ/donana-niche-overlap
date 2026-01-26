rm(list = ls())

# ==============================================================================
# Script: 03_spatiotemporal_analysis.R
# ==============================================================================
# Description: Spatiotemporal analysis of predator-prey interactions in Doñana
# National Park. We applied the recurrent event analysis framework for camera
# trap data proposed by Ferry et al. (2024). This approach allows testing the
# effect of a primary species on the detection of a secondary species,
# assessing potential attraction or avoidance patterns.
# ==============================================================================

# 1. Loading packages -----------------------------------------------------
library("data.table")
library("lubridate")
library("dplyr")
library("ctrecurrent")
library("pammtools")
library("ggplot2")
library("patchwork")

# ==============================================================================
# 2. Loading and formatting data -----------------------------------------------
# ==============================================================================

# load("RData/03_spatiotemporal_analysis.RData")

data <- fread("data/dataNicheOverlap.csv")
data$date_time <- ymd_hms(data$date_time, tz = "UTC")
data <- data %>% rename(
    DateTime = date_time,
    Site = site,
    Species = sp
)

# We include human detections (researchers) as a tertiary species to account for
# potential confounding effects related to field maintenance. Since site visits
# for SD card replacement occur every two months and may influence species
# behavior, incorporating researcher presence in the recurrent event analysis
# allows us to control for these disturbances. Additionally, these detections
# help confirm survey end times in cases of camera malfunction.


# Load the start time of each camera deployment, corresponding to field visits
# for SD card replacement.
researchers <- fread("data/detectionsResearchers.csv")
researchers$DateTime <- dmy_hm(researchers$DateTime, tz = "UTC")
data <- rbind(data, researchers)

# ==============================================================================
# 3. Obtain recurrent events of the secondary species --------------------------
# Predators and prey were considered as primary species to assess bidirectional
# effects.
# ==============================================================================

#
# Predators as primary species
primary_predators <- c(
    "Vulpes vulpes",
    "Meles meles",
    "Genetta genetta",
    "Herpestes ichneumon",
    "Lynx pardinus"
)

secondary_prey <- c(
    "Lepus granatensis",
    "Oryctolagus cuniculus"
)

survey_duration <- 30

# Get recurrent events for each predator prey combination
predator_prey_recu <- list()
for (p in primary_predators) {
    for (s in secondary_prey) {
        data.recu <- ct_to_recurrent(
            data = data,
            primary = p,
            secondary = s,
            tertiary = NULL, # Treat all other species as tertiary to account for background activity
            survey_end_date = max(data$DateTime),
            survey_duration = survey_duration,
            species_var = "Species",
            datetime_var = "DateTime"
        )

        predator_prey_recu[[paste(p, s, sep = " - ")]] <- data.recu
    }
}

#
# Prey as primary species
primary_prey <- c(
    "Lepus granatensis",
    "Oryctolagus cuniculus"
)
secondary_predators <- c(
    "Vulpes vulpes",
    "Meles meles",
    "Genetta genetta",
    "Herpestes ichneumon",
    "Lynx pardinus"
)

# Generate recurrent event data for prey-predator combinations
prey_predator_recu <- list()

for (p in primary_prey) {
    for (s in secondary_predators) {
        data.recu <- ct_to_recurrent(
            data = data,
            primary = p,
            secondary = s,
            tertiary = NULL, # Treat all other species as tertiary to account for background activity
            survey_end_date = max(data$DateTime),
            survey_duration = survey_duration,
            species_var = "Species",
            datetime_var = "DateTime"
        )

        prey_predator_recu[[paste(p, s, sep = " - ")]] <- data.recu
    }
}

# ==============================================================================
# 4. Convert recurrent events to PED format ------------------------------------
# ==============================================================================

#
# PED format for predators as primary species
predator_prey_ped <- list()
for (i in seq_along(predator_prey_recu)) {
    # Process each recurrent event dataset
    df_tmp <- predator_prey_recu[[i]]
    name_tmp <- names(predator_prey_recu)[i]

    # Skip if the dataframe is empty or has no events
    if (nrow(df_tmp) == 0) {
        message("Skipping ", name_tmp, " (no events or empty dataset)")
        next
    }
    ped <- df_tmp %>%
        as_ped(
            formula    = Surv(t.start, t.stop, event) ~ Site,
            id         = "survey_id",
            transition = "enum",
            timescale  = "calendar"
        ) %>%
        mutate(
            Site = as.factor(Site),
            survey_id = as.factor(survey_id)
        )
    predator_prey_ped[[name_tmp]] <- ped
}

#
# PED format for prey as primary species
prey_predator_ped <- list()
for (i in seq_along(prey_predator_recu)) {
    # Process each recurrent event dataset
    df_tmp <- prey_predator_recu[[i]]
    name_tmp <- names(prey_predator_recu)[i]

    # Skip if the dataframe is empty or has no events
    if (nrow(df_tmp) == 0) {
        message("Skipping ", name_tmp, " (no events or empty dataset)")
        next
    }
    ped <- df_tmp %>%
        as_ped(
            formula    = Surv(t.start, t.stop, event) ~ Site,
            id         = "survey_id",
            transition = "enum",
            timescale  = "calendar"
        ) %>%
        mutate(
            Site = as.factor(Site),
            survey_id = as.factor(survey_id)
        )
    prey_predator_ped[[name_tmp]] <- ped
}


# ==============================================================================
# 5. Model fitting -------------------------------------------------------------
# Fit Piece-wise Additive Mixed Models (PAMMs) to recurrent event data
# ==============================================================================

#
# Fit the models for predators as primary species
predator_prey_pamm <- list()

for (i in seq_along(predator_prey_ped)) {
    name_pamm <- names(predator_prey_ped)[i]
    dat <- predator_prey_ped[[i]]

    message("➡️ Processing: ", name_pamm)

    # Helper function to attempt model fitting and handle errors
    fit_attempt <- function(formula, data, offset) {
        tryCatch(
            {
                fit <- pamm(
                    formula = formula,
                    data = data,
                    offset = offset,
                    engine = "bam",
                    method = "fREML",
                    discrete = TRUE
                )
                return(fit) # Return the fitted model object if successful
            },
            error = function(e) {
                message("❌ (Error) ", substring(e$message, 1, 200))
                return(NULL) # Return NULL if fitting fails
            }
        )
    }

    # 1) First attempt
    fit1 <- fit_attempt(ped_status ~ s(tend) + s(Site, bs = "re"), dat, offset)

    if (!is.null(fit1)) {
        predator_prey_pamm[[name_pamm]] <- fit1
        message("✅ Model fitted correctly: ", name_pamm)
        next
    }

    # 2) If the first attempt fails, check the error message to decide the fallback
    # Trying with k = 3
    message("⚠️ First attempt failed for ", name_pamm, " → Reattempting with k = 3")
    fit2 <- fit_attempt(ped_status ~ s(tend, k = 3) + s(Site, bs = "re"), dat, offset)
    if (!is.null(fit2)) {
        predator_prey_pamm[[name_pamm]] <- fit2
        message("✅ Model fitted with k = 3: ", name_pamm)
        next
    }

    # 3) If it still fails, try without the Site term
    message("⚠️ Failed with k=3 for ", name_pamm, " → Reattempting without Site")
    fit3 <- fit_attempt(ped_status ~ s(tend), dat, offset)
    if (!is.null(fit3)) {
        predator_prey_pamm[[name_pamm]] <- fit3
        message("✅ Model fitted without Site: ", name_pamm)
        next
    }

    # 4) If all fails, save NULL and continue
    predator_prey_pamm[[name_pamm]] <- NULL
    message("❌ All attempts failed for ", name_pamm)
}


#
# Fit the models for prey as primary species
prey_predator_pamm <- list()

for (i in seq_along(prey_predator_ped)) {
    name_pamm <- names(prey_predator_ped)[i]
    dat <- prey_predator_ped[[i]]

    message("➡️ Processing: ", name_pamm)

    # Function to attempt fitting and return fit or NULL
    fit_attempt <- function(formula, data, offset) {
        tryCatch(
            {
                fit <- pamm(
                    formula = formula,
                    data = data,
                    offset = offset,
                    engine = "bam",
                    method = "fREML",
                    discrete = TRUE
                )
                return(fit) # Return the fitted model object if successful
            },
            error = function(e) {
                message("❌ (Error) ", substring(e$message, 1, 200))
                return(NULL) # Return NULL if fitting fails
            }
        )
    }

    # 1) First attempt
    fit1 <- fit_attempt(ped_status ~ s(tend) + s(Site, bs = "re"), dat, offset)

    if (!is.null(fit1)) {
        prey_predator_pamm[[name_pamm]] <- fit1
        message("✅ Model fitted correctly: ", name_pamm)
        next
    }

    # 2) If the first attempt fails, check the error message to decide the fallback
    # Trying with k = 3
    message("⚠️ First attempt failed for ", name_pamm, " → Reattempting with k = 3")
    fit2 <- fit_attempt(ped_status ~ s(tend, k = 3) + s(Site, bs = "re"), dat, offset)
    if (!is.null(fit2)) {
        prey_predator_pamm[[name_pamm]] <- fit2
        message("✅ Model fitted with k = 3: ", name_pamm)
        next
    }

    # 3) If it still fails, try without the Site term
    fit3 <- fit_attempt(ped_status ~ s(tend), dat, offset)
    if (!is.null(fit3)) {
        prey_predator_pamm[[name_pamm]] <- fit3
        message("✅ Model fitted without Site: ", name_pamm)
        next
    }
    # 4) If all fails, save NULL and continue
    prey_predator_pamm[[name_pamm]] <- NULL
    message("❌ All attempts failed for ", name_pamm)
}

# ==============================================================================
# 6. Model summary -------------------------------------------------------------
# Extract model summaries and store them in a list. We retain the effective
# degrees of freedom (edf), reference degrees of freedom (Ref.df),
# Chi-squared statistics, and p-values.
# ==============================================================================

results_pred_prey <- vector(
    "list",
    length(predator_prey_pamm)
)

names(results_pred_prey) <- names(predator_prey_pamm)

for (i in seq_along(predator_prey_pamm)) {
    model_i <- summary(
        predator_prey_pamm[[i]]
    )
    # Check if the model exists and contains an s.table
    if (!is.null(model_i) && !is.null(model_i$s.table)) {
        s_table <- as.data.frame(model_i$s.table)
        s_table$pairs <- names(predator_prey_pamm)[i]
        s_table <- s_table %>%
            mutate(
                edf = round(edf, 2),
                Ref.df = round(Ref.df, 2),
                Chi.sq = round(Chi.sq, 2),
                `p-value` = round(`p-value`, 2)
            )
        results_pred_prey[[i]] <- s_table
    } else {
        results_pred_prey[[i]] <- NULL
        warning(paste("No s.table for element", i, names(predator_prey_pamm)[i]))
    }
}

results_pred_prey <- bind_rows(results_pred_prey)

results_prey_predator <- vector("list", length(prey_predator_pamm))
names(results_prey_predator) <- names(prey_predator_pamm)

for (i in seq_along(prey_predator_pamm)) {
    model_i <- summary(prey_predator_pamm[[i]])
    # Check if the model exists and contains an s.table
    if (!is.null(model_i) && !is.null(model_i$s.table)) {
        s_table <- as.data.frame(model_i$s.table)
        s_table$pairs <- names(prey_predator_pamm)[i]
        s_table <- s_table %>%
            mutate(
                edf = round(edf, 2),
                Ref.df = round(Ref.df, 2),
                Chi.sq = round(Chi.sq, 2),
                `p-value` = round(`p-value`, 2)
            )
        results_prey_predator[[i]] <- s_table
    } else {
        results_prey_predator[[i]] <- NULL
        warning(paste("No s.table for element", i, names(prey_predator_pamm)[i]))
    }
}

results_prey_predator <- bind_rows(results_prey_predator)

print(results_prey_predator)

# ==============================================================================
# 7. Generate Predictions and Plots --------------------------------------------
# Create a matrix of plots for each predator-prey pair. Each plot displays the
# predicted log hazard probabilities of species presence over time. Note that
# some pairs (e.g., Lynx - Hare, Genet - Rabbit) may lack sufficient data for
# prediction, resulting in empty plots.
# ==============================================================================


# Define species
lep <- c(
    "Lepus granatensis",
    "Oryctolagus cuniculus"
)
car <- c(
    "Vulpes vulpes",
    "Meles meles",
    "Genetta genetta",
    "Herpestes ichneumon",
    "Lynx pardinus"
)

# Create an empty list to store the plots
plot_matrix <- list()

# Iterate over carnivores and leporids
for (carn in car) {
    row_plots <- list() # list of plots for one row

    for (lep_sp in lep) {
        # Expected names in the lists
        prey_pred <- paste(lep_sp, "-", carn)
        pred_prey <- paste(carn, "-", lep_sp)

        # Check if they exist in the lists
        if (prey_pred %in% names(prey_predator_pamm) && pred_prey %in% names(predator_prey_pamm)) {
            # =================================================================
            # Extract data when predators are the primary species
            pred_model <- predator_prey_pamm[[pred_prey]]
            # Get the summary of the model. If the model is significant
            # we will use dashed lines, otherwise solid lines
            sig_pred_prey <- summary(pred_model)
            smooth_terms_pred_prey <- as.data.frame(sig_pred_prey$s.table)
            pval_pred_prey <- smooth_terms_pred_prey["s(tend)", "p-value"]
            lt_pred_prey <- ifelse(pval_pred_prey < 0.05, "dashed", "solid")

            # Extract data when leporids are the primary species
            prey_model <- prey_predator_pamm[[prey_pred]]
            # Get the summary of the model. If the model is significant
            # we will use dashed lines, otherwise solid lines
            sig_prey_pred <- summary(prey_model)
            smooth_terms_prey_pred <- as.data.frame(sig_prey_pred$s.table)
            pval_prey_pred <- smooth_terms_prey_pred["s(tend)", "p-value"]
            lt_prey_pred <- ifelse(pval_prey_pred < 0.05, "dashed", "solid")
            # =================================================================


            # =================================================================
            # Make predictions
            pred_ped <- predator_prey_ped[[pred_prey]]
            prey_ped <- prey_predator_ped[[prey_pred]]

            ndf_null_pred_prey <- pred_ped %>%
                make_newdata(tend = unique(tend)) %>%
                add_hazard(pred_model)

            ndf_null_prey_pred <- prey_ped %>%
                make_newdata(tend = unique(tend)) %>%
                add_hazard(prey_model)
            # =================================================================

            p <- ggplot() +
                geom_line(
                    data = ndf_null_pred_prey,
                    aes(
                        x = tend,
                        y = hazard,
                        colour = "Predator as primary sp"
                    ),
                    linewidth = 0.8,
                    linetype = lt_pred_prey
                ) +
                geom_ribbon(
                    data = ndf_null_pred_prey,
                    aes(
                        x = tend,
                        ymin = ci_lower,
                        ymax = ci_upper,
                        fill = "Predator as primary sp"
                    ),
                    alpha = .3
                ) +
                geom_line(
                    data = ndf_null_prey_pred,
                    aes(
                        x = tend,
                        y = hazard,
                        colour = "Prey as primary sp"
                    ),
                    linewidth = 0.8,
                    linetype = lt_prey_pred
                ) +
                geom_ribbon(
                    data = ndf_null_prey_pred,
                    aes(
                        x = tend,
                        ymin = ci_lower,
                        ymax = ci_upper,
                        fill = "Prey as primary sp"
                    ),
                    alpha = .3
                ) +
                scale_color_manual(
                    name = "Direction",
                    values = c(
                        "Pred as primary sp" = "steelblue",
                        "Prey as primary sp" = "firebrick"
                    )
                ) +
                scale_fill_manual(
                    name = "Primary sp",
                    values = c(
                        "Predator as primary sp" = "steelblue",
                        "Prey as primary sp" = "firebrick"
                    ),
                    labels = c("Predators", "Leporids")
                ) +
                scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
                scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
                labs(
                    x = NULL,
                    y = NULL,
                    title = NULL
                ) +
                theme_minimal(base_size = 8) +
                theme(
                    panel.grid = element_blank(),
                    axis.ticks = element_line(color = "black"),
                    axis.line = element_blank(),
                    panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        linewidth = 0.6
                    )
                )
        } else {
            # For pairs without data, create an empty plot with invisible axes
            p <- ggplot() +
                scale_x_continuous(
                    limits = c(0, 30),
                    labels = scales::number_format(accuracy = 1)
                ) +
                scale_y_continuous(
                    limits = c(0, 5),
                    labels = scales::number_format(accuracy = 0.1)
                ) +
                theme_minimal(base_size = 8) +
                theme(
                    panel.grid = element_blank(),
                    axis.text = element_text(color = "white"),
                    axis.ticks = element_line(color = "white"),
                    axis.line = element_blank(),
                    panel.border = element_rect(
                        color = "black",
                        fill = NA,
                        linewidth = 0.6
                    )
                )
        }

        row_plots[[lep_sp]] <- p
    }

    # wrap the plots
    plot_matrix[[carn]] <- wrap_plots(row_plots, ncol = length(lep))
}

plot_matrix <- plot_matrix[c(1, 2, 4, 3, 5)]

# Combine the plots
final_panel <- wrap_plots(plot_matrix, ncol = 1) +
    plot_layout(guides = "collect") &
    theme(
        text = element_text(family = "serif"),
        legend.position = "right",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        axis.text = element_text(size = 13)
    ) &
    guides(colour = "none")

print(final_panel)

# save the plot
ggsave(
    file = paste(
        "results/",
        "figure 4 - Spatio-temporal interactions.png",
        sep = ""
    ),
    plot = final_panel,
    width = 9,
    height = 12,
    dpi = 500
)

save.image("RData/03_spatiotemporal_analysis.RData")
