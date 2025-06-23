#!/usr/bin/env Rscript
#==============================================================================
# Cumulative Effect Analysis Script
#
# This script:
#  1. Loads required libraries.
#  2. Loads the time-series data and checks for necessary columns.
#     It creates a 'region_index' column if missing.
#  3. Loads MCMC model results for Model 1 and Model 3.
#  4. Generates cumulative effect plots for Temperature, Humidity, and Rainfall,
#     by Local Health District (LHD), using the extracted model coefficients.
#  5. Saves the resulting plots to a PDF file.
#
# Adjust file paths as needed.
#==============================================================================

#------------------------------------------------------------------------------
# 1. Load Required Libraries
#------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(dlnm)
library(viridis)
library(gridExtra)
library(grid)
library(scales)

#------------------------------------------------------------------------------
# 2. Load and Prepare Data
#------------------------------------------------------------------------------
data_file <- "new_output/data_timeseries.RData"
if (!file.exists(data_file)) {
  stop("Data file does not exist: ", data_file)
}
load(data_file)  # Expecting to load an object called 'data_ts'
if (!exists("data_ts")) {
  stop("The object 'data_ts' was not found in the data file.")
}

# Check for required columns in data_ts.
required_cols <- c("LHD", "Cases", "temperature_2m", "relative_humidity_2m", "rain")
missing_cols <- setdiff(required_cols, colnames(data_ts))
if (length(missing_cols) > 0) {
  stop("The following required columns are missing in data_ts: ", paste(missing_cols, collapse=", "))
}

# If 'region_index' is missing, create it from the 'LHD' column.
if (!"region_index" %in% colnames(data_ts)) {
  data_ts$region_index <- as.numeric(factor(data_ts$LHD))
  cat("Created 'region_index' from 'LHD' column.\n")
}

#------------------------------------------------------------------------------
# 3. Load Model Results
#------------------------------------------------------------------------------
# We assume that model results for Model 1 and Model 3 are stored in separate files.
models <- list()
model1_file <- "new_output_acc/model1_results.RData"
model3_file <- "new_output_acc/model3_results.RData"

if (!file.exists(model1_file)) {
  stop("Model1 results file not found: ", model1_file)
}
tmp <- load(model1_file)
# 'tmp' is a character vector with the names of the loaded objects.
models[[1]] <- get(tmp[1])

if (!file.exists(model3_file)) {
  stop("Model3 results file not found: ", model3_file)
}
tmp <- load(model3_file)
models[[3]] <- get(tmp[1])

#------------------------------------------------------------------------------
# 4. Helper Functions
#------------------------------------------------------------------------------

# Function to extract DLNM coefficients from model results that match a given pattern.
extract_dlnm_effects <- function(model_results, pattern) {
  if (is.null(model_results)) return(NULL)
  
  # Convert mcmc.list to a matrix if needed.
  samples_matrix <- if (inherits(model_results, "mcmc.list")) {
    as.matrix(model_results)
  } else {
    do.call(rbind, lapply(model_results, as.matrix))
  }
  
  # Ensure that 'pattern' is a single valid character string.
  if (is.null(pattern) || !is.character(pattern) || length(pattern) != 1) {
    cat("Invalid pattern provided:", pattern, "\n")
    return(NULL)
  }
  
  coef_cols <- grep(pattern, colnames(samples_matrix), value = TRUE)
  if (length(coef_cols) == 0) return(NULL)
  
  samples_matrix[, coef_cols, drop = FALSE]
}

# Function to create a cumulative effect plot for a given variable and LHD.
create_cumulative_plot <- function(data, model_results, variable, lhd_name, n_cases) {
  # Subset data for the specified LHD.
  lhd_data <- data[data$LHD == lhd_name, ]
  if (nrow(lhd_data) == 0) {
    cat("No data for LHD:", lhd_name, "\n")
    return(NULL)
  }
  
  # Ensure a valid region index exists and extract the first value.
  if (!"region_index" %in% colnames(lhd_data)) {
    cat("The column 'region_index' is missing for LHD:", lhd_name, "\n")
    return(NULL)
  }
  region_idx_raw <- unique(lhd_data$region_index)[1]
  region_idx <- as.numeric(as.character(region_idx_raw))
  if (length(region_idx) == 0 || is.na(region_idx)) {
    cat("Invalid region index for LHD:", lhd_name, "\n")
    return(NULL)
  }
  
  # Set a reference point for centering predictions.
  ref_point <- switch(variable,
                      "temperature_2m" = 20,
                      "relative_humidity_2m" = 60,
                      "rain" = 300)
  
  # Create the crossbasis for the variable using natural splines.
  var_knots <- quantile(lhd_data[[variable]], c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    lhd_data[[variable]],
    lag = 1,
    argvar = list(fun = "ns", knots = var_knots),
    arglag = list(fun = "ns", df = 2)
  )
  
  # Set up prediction points.
  var_range <- range(lhd_data[[variable]], na.rm = TRUE)
  var_points <- seq(var_range[1], var_range[2], length.out = 100)
  
  # Build the search pattern for coefficient names.
  # For example, if variable is temperature_2m and region_idx is 1,
  # the pattern becomes "b_temp\\[1,".
  pattern <- switch(variable,
                    "temperature_2m" = sprintf("b_temp\\[%d,", region_idx),
                    "relative_humidity_2m" = sprintf("b_hum\\[%d,", region_idx),
                    "rain" = sprintf("b_rain\\[%d,", region_idx),
                    NA_character_)
  if (is.na(pattern)) {
    cat("Variable", variable, "not recognized. Skipping.\n")
    return(NULL)
  }
  
  # Extract the DLNM effects (coefficients) from both models.
  effects1 <- extract_dlnm_effects(model_results[[1]], pattern)
  effects3 <- extract_dlnm_effects(model_results[[3]], pattern)
  
  if (is.null(effects1) || is.null(effects3)) {
    cat(sprintf("No effects found for %s in LHD %s\n", variable, lhd_name))
    return(NULL)
  }
  
  # Function to compute cumulative predictions using crosspred().
  get_predictions <- function(effects) {
    predictions <- matrix(NA, nrow = length(var_points), ncol = nrow(effects))
    for (i in seq_len(nrow(effects))) {
      pred <- crosspred(
        cb,
        coef = effects[i, ],
        vcov = matrix(0, ncol(effects), ncol(effects)),  # assume no uncertainty for predictions
        model.link = "log",
        at = var_points,
        cen = ref_point,
        cumul = TRUE
      )
      predictions[, i] <- pred$allRRfit
    }
    data.frame(
      mean = apply(predictions, 1, mean, na.rm = TRUE),
      lower = apply(predictions, 1, quantile, probs = 0.025, na.rm = TRUE),
      upper = apply(predictions, 1, quantile, probs = 0.975, na.rm = TRUE)
    )
  }
  
  # Obtain predictions for both sets of effects.
  pred1 <- get_predictions(effects1)
  pred3 <- get_predictions(effects3)
  
  if (all(is.na(pred1$mean)) || all(is.na(pred3$mean))) return(NULL)
  
  # Set labels for annotation.
  ref_label <- switch(variable,
                      "temperature_2m" = "20°C",
                      "relative_humidity_2m" = "60%",
                      "rain" = "300mm")
  var_label <- switch(variable,
                      "temperature_2m" = "Temperature (°C)",
                      "relative_humidity_2m" = "Humidity (%)",
                      "rain" = "Rainfall (mm)")
  formatted_cases <- format(n_cases, big.mark = ",")
  
  # Assemble plot data.
  plot_data <- data.frame(
    exposure = var_points,
    m1_mean = pred1$mean,
    m1_lower = pred1$lower,
    m1_upper = pred1$upper,
    m3_mean = pred3$mean,
    m3_lower = pred3$lower,
    m3_upper = pred3$upper
  )
  
  # Generate the cumulative effect plot.
  ggplot(plot_data, aes(x = exposure)) +
    geom_ribbon(aes(ymin = m1_lower, ymax = m1_upper),
                fill = "#FDB863", alpha = 0.3) +
    geom_ribbon(aes(ymin = m3_lower, ymax = m3_upper),
                fill = "#B2ABD2", alpha = 0.3) +
    geom_line(aes(y = m1_mean, color = "Independent CCO"), linewidth = 0.8) +
    geom_line(aes(y = m3_mean, color = "Spatial CCO"), linewidth = 0.8) +
    geom_vline(xintercept = ref_point, linetype = "dotted") +
    annotate("text", x = ref_point, y = 0.5, 
             label = ref_label, angle = 90, 
             vjust = -0.5, size = 2.5, hjust = 0) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_color_manual(
      values = c("Independent CCO" = "#FDB863", "Spatial CCO" = "#B2ABD2")
    ) +
    scale_y_continuous(
      trans = "log",
      breaks = c(0.5, 1, 2, 3, 4, 5, 6),
      limits = c(0.5, 6)
    ) +
    labs(
      title = sprintf("%s (%s cases)", lhd_name, formatted_cases),
      x = var_label,
      y = "Cumulative Relative Risk",
      color = ""
    ) +
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.justification = "center",
      legend.margin = margin(b = 5),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.8, "lines"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 9, face = "bold"),
      plot.margin = margin(t = 15, r = 15, b = 15, l = 15)
    )
}

#------------------------------------------------------------------------------
# 5. Create and Save Cumulative Effect Plots to PDF
#------------------------------------------------------------------------------
create_cumulative_plots <- function() {
  cat("\nCreating cumulative effect plots...\n")
  
  # Get unique Local Health Districts (LHDs) and aggregate cases by LHD.
  lhds <- unique(data_ts$LHD)
  cases_by_lhd <- aggregate(Cases ~ LHD, data = data_ts, sum)
  
  # Define variables of interest and their pretty names.
  variables <- c("temperature_2m", "relative_humidity_2m", "rain")
  var_names <- c("Temperature", "Humidity", "Rainfall")
  
  # Open a PDF device to save the plots.
  pdf("Visualisations/cumulative_effects.pdf", width = 16, height = 12)
  
  # Loop over each variable.
  for (var_idx in seq_along(variables)) {
    cat(sprintf("\nProcessing %s...\n", var_names[var_idx]))
    plots <- list()
    
    # Create a cumulative effect plot for each LHD.
    for (i in seq_along(lhds)) {
      plots[[i]] <- create_cumulative_plot(
        data = data_ts,
        model_results = models,
        variable = variables[var_idx],
        lhd_name = lhds[i],
        n_cases = cases_by_lhd$Cases[cases_by_lhd$LHD == lhds[i]]
      )
    }
    
    # Remove any NULL plots (in case extraction failed for some LHDs).
    plots <- plots[!sapply(plots, is.null)]
    
    if (length(plots) > 0) {
      # Calculate number of pages needed (6 plots per page).
      n_plots <- length(plots)
      n_pages <- ceiling(n_plots / 6)
      
      # Arrange and plot pages.
      for (page in 1:n_pages) {
        start_idx <- (page - 1) * 6 + 1
        end_idx <- min(page * 6, n_plots)
        page_plots <- plots[start_idx:end_idx]
        
        grid.arrange(
          grobs = c(
            list(textGrob(
              sprintf("Cumulative %s-Influenza Spread by LHD (Page %d/%d)", 
                      var_names[var_idx], page, n_pages),
              gp = gpar(fontsize = 16, fontface = "bold"),
              vp = viewport(height = unit(2, "lines"))
            )),
            page_plots
          ),
          layout_matrix = rbind(
            c(1, 1, 1),
            matrix(2:7, ncol = 2, byrow = TRUE)
          ),
          heights = c(0.08, rep(1, 3))
        )
      }
    }
  }
  
  dev.off()
  cat("\nVisualization complete! Check cumulative_effects.pdf\n")
}

#------------------------------------------------------------------------------
# 6. Main Execution
#------------------------------------------------------------------------------
cat("Loading data and models...\n")
# Data and models were loaded and prepared above.
create_cumulative_plots()
