# Enhanced Climate Model Analysis with Additional Features
# ------------------------------------------------------------------------------
# This script performs an advanced analysis of climate-disease associations by
# leveraging MCMC model results. It incorporates multiple features:
#
#   1) Risk analysis and mapping: Computes relative risk estimates for different
#      Local Health Districts (LHDs) and visualizes these risks on maps.
#   2) Probability mapping and uncertainty quantification: Calculates the 
#      probability that the relative risk exceeds a specified threshold.
#   3) Temporal/spatial statistics: Computes measures such as autocorrelation,
#      seasonal decomposition, and spatial clustering (e.g., Moran's I).
#   4) Enhanced error handling and validation: Checks data integrity and gracefully
#      handles errors.
#
# This script uses various R packages for spatial analysis (sf), plotting (ggplot2,
# plotly, gridExtra), MCMC diagnostics (coda), and more.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Load required libraries
# ------------------------------------------------------------------------------
library(ggplot2)      # For advanced plotting
library(dplyr)        # For data manipulation and transformation
library(tidyr)        # For reshaping data
library(patchwork)    # For combining multiple plots into panels
library(viridis)      # For color scales with perceptually uniform palettes
library(RColorBrewer) # For additional color palettes
library(lubridate)    # For handling date-time data
library(gridExtra)    # For arranging plots in grids
library(coda)         # For MCMC sample analysis and diagnostics
library(sf)           # For handling spatial data (simple features)
library(plotly)       # For interactive plotting
library(corrplot)     # For visualizing correlation matrices
library(car)          # For regression diagnostics
library(broom)        # For tidying model output
library(ape)          # For spatial autocorrelation measures, e.g., Moran's I

# ------------------------------------------------------------------------------
# Data Validation Utilities
# ------------------------------------------------------------------------------
# The following utilities ensure that data components have the expected structure
# and that the mapping between LHD names in the dataset and shapefile is standardized.

# Create a bridge table to match LHD names from the dataset with those in the shapefile.
bridge_table <- data.frame(
  data_LHD = toupper(c(
    "Central Coast", "Far West", "Hunter New England", "Illawarra Shoalhaven",
    "Mid North Coast", "Murrumbidgee", "Nepean Blue Mountains", "Northern NSW",
    "Northern Sydney", "South Eastern Sydney", "South Western Sydney",
    "Southern NSW", "Sydney", "Western NSW", "Western Sydney"
  )),
  shape_LHD = toupper(c(
    "CENTRAL COAST", "FAR WEST", "HUNTER NEW ENGLAND", "ILLAWARRA SHOALHAVEN",
    "MID NORTH COAST", "MURRUMBIDGEE", "NEPEAN BLUE MOUNTAINS", "NORTHERN NSW",
    "NORTHERN SYDNEY", "SOUTH EASTERN SYDNEY", "SOUTH WESTERN SYDNEY",
    "SOUTHERN NSW", "SYDNEY", "WESTERN NSW", "WESTERN SYDNEY"
  ))
)

# load_all_data: Loads time series data, trend/seasonality components, shapefile,
# and model results. Also standardizes LHD names by merging with the bridge table.
load_all_data <- function() {
  # Load base time series data and additional components
  load("new_output/data_timeseries.RData")
  load("new_output/trend_timeseries.RData")
  load("new_output/seasonality_timeseries.RData")
  
  # Load the spatial boundaries shapefile for NSW Local Health Districts
  nsw_lhd <- st_read("spatial_data/NSW_LHD_Boundaries.shp", quiet = TRUE)
  
  # Standardize LHD names in the time series data
  data_ts <- data_ts %>%
    mutate(
      LHD = toupper(lhd_code),
      date = as.Date(date)
    ) %>%
    left_join(bridge_table, by = c("LHD" = "data_LHD")) %>%
    mutate(LHD = shape_LHD)
  
  # Load MCMC model results from different models (model1 through model4)
  model_results <- list()
  for(i in 1:4) {
    tryCatch({
      load(sprintf("new_output_acc/model%d_results.RData", i))
      model_results[[paste0("model", i)]] <- get(sprintf("model%d_samples", i))
    }, error = function(e) {
      cat(sprintf("Error loading model %d: %s\n", i, e$message))
    })
  }
  
  # Standardize LHD names in the shapefile to uppercase for consistency
  nsw_lhd <- nsw_lhd %>%
    mutate(LHD = toupper(lhd_name))
  
  return(list(
    data_ts = data_ts,
    trend_ts = trend_ts,
    seas_ts = seas_ts,
    model_results = model_results,
    nsw_map = nsw_lhd
  ))
}

# validate_data: Checks that the essential data components exist and have the
# expected structure. Stops execution if critical components are missing.
validate_data <- function(data) {
  required_components <- c("nsw_map", "model_results", "data_ts")
  missing_components <- setdiff(required_components, names(data))
  
  if (length(missing_components) > 0) {
    stop("Missing required data components: ", 
         paste(missing_components, collapse = ", "))
  }
  
  # Ensure that the spatial map is an sf object and contains the LHD field.
  if (!inherits(data$nsw_map, "sf")) {
    stop("nsw_map must be an sf object")
  }
  if (!"LHD" %in% names(data$nsw_map)) {
    stop("nsw_map must contain LHD column")
  }
  
  # Check that model results are available as a non-empty list.
  if (!is.list(data$model_results) || length(data$model_results) == 0) {
    stop("model_results must be a non-empty list")
  }
  
  # Verify that the time series data has the required columns.
  if (!is.null(data$data_ts)) {
    required_cols <- c("temperature_2m", "relative_humidity_2m", "rain")
    missing_cols <- setdiff(required_cols, names(data$data_ts))
    if (length(missing_cols) > 0) {
      warning("Missing time series columns: ", 
              paste(missing_cols, collapse = ", "))
    }
  }
  
  TRUE
}

# ------------------------------------------------------------------------------
# 1) Spatial/Temporal Statistics
# ------------------------------------------------------------------------------
# calculate_spatial_statistics: Computes spatial statistics for risk estimates.
# Uses the centroids of polygons (LHDs) to compute distance matrices and weight
# matrices, then calculates Moran's I for spatial autocorrelation and performs 
# hierarchical clustering.
calculate_spatial_statistics <- function(risks) {
  if (!inherits(risks, "sf")) {
    stop("risks must be an sf object")
  }
  if (!"risk" %in% names(risks)) {
    stop("risks must contain a 'risk' column")
  }
  
  tryCatch({
    # Compute centroids and distances between them.
    coords <- st_coordinates(st_centroid(risks))
    dist_matrix <- as.matrix(dist(coords))
    weights <- 1 / dist_matrix
    diag(weights) <- 0  # Set self-weights to zero
    
    # Standardize weights by dividing by the row sums.
    weights <- weights / rowSums(weights)
    
    # Compute Moran's I to quantify spatial autocorrelation.
    morans_i <- ape::Moran.I(risks$risk, weights)
    
    # Hierarchical clustering of regions based on distances.
    hc <- hclust(dist(coords), method = "ward.D2")
    clusters <- cutree(hc, k = 4)
    
    list(
      morans_i = morans_i,
      clusters = clusters,
      success = TRUE
    )
  }, error = function(e) {
    warning("Error in spatial statistics calculation: ", e$message)
    list(success = FALSE, error = e$message)
  })
}

# calculate_temporal_statistics: Computes temporal statistics including
# autocorrelation functions and seasonal decomposition for key climate variables.
calculate_temporal_statistics <- function(data) {
  if (!all(c("temperature_2m", "relative_humidity_2m", "rain") %in% names(data))) {
    stop("Missing required columns in temporal data")
  }
  
  tryCatch({
    # Compute autocorrelation functions for each climate variable.
    acf_results <- list(
      temp = acf(data$temperature_2m, plot = FALSE),
      hum  = acf(data$relative_humidity_2m, plot = FALSE),
      rain = acf(data$rain, plot = FALSE)
    )
    
    # Perform seasonal decomposition using time series (frequency = 12 for monthly data).
    decomp_results <- list(
      temp = decompose(ts(data$temperature_2m, frequency = 12)),
      hum  = decompose(ts(data$relative_humidity_2m, frequency = 12)),
      rain = decompose(ts(data$rain, frequency = 12))
    )
    
    list(
      acf = acf_results,
      decomposition = decomp_results,
      success = TRUE
    )
  }, error = function(e) {
    warning("Error in temporal statistics calculation: ", e$message)
    list(success = FALSE, error = e$message)
  })
}

# ------------------------------------------------------------------------------
# 2) Probability Mapping Utilities
# ------------------------------------------------------------------------------
# calculate_exceedance_probabilities: Calculates the probability that the
# relative risk (RR) for an exposure exceeds a given threshold. The function 
# loops over regions (LHDs) and extracts coefficients from MCMC samples to compute
# mean effects and then the exceedance probability.
calculate_exceedance_probabilities <- function(model_results, variable,
                                               threshold = 1.2, n_regions = 15) {
  if (!exists("bridge_table")) {
    stop("bridge_table not found in environment")
  }
  
  tryCatch({
    # Convert MCMC samples to a matrix for easier manipulation.
    samples_matrix <- do.call(rbind, lapply(model_results, as.matrix))
    
    # Identify the appropriate parameter pattern based on the variable.
    pattern <- switch(variable,
                      "temperature" = "^b_temp\\[",
                      "humidity"    = "^b_hum\\[",
                      "rain"        = "^b_rain\\[",
                      stop("Invalid variable specified"))
    
    exceedance_prob <- numeric(n_regions)
    
    for(i in seq_len(n_regions)) {
      # Create a region-specific pattern and extract columns.
      region_pattern <- sprintf("%s%d,", pattern, i)
      region_cols <- grep(region_pattern, colnames(samples_matrix), value = TRUE)
      
      if (length(region_cols) > 0) {
        # Compute mean effects across MCMC samples and determine probability.
        region_effects <- samples_matrix[, region_cols]
        mean_effects <- rowMeans(region_effects)
        exceedance_prob[i] <- mean(exp(mean_effects) > threshold)
      } else {
        exceedance_prob[i] <- NA
      }
    }
    
    data.frame(
      LHD = bridge_table$shape_LHD,
      probability = exceedance_prob
    )
  }, error = function(e) {
    warning("Error in probability calculation: ", e$message)
    NULL
  })
}

# plot_probability_map: Plots a map showing the probability that the RR exceeds
# a specified threshold for a given variable. The function joins the probability
# data with the spatial map and visualizes it using a color scale.
plot_probability_map <- function(model_results, map_data, variable,
                                 threshold = 1.2) {
  # Calculate exceedance probabilities for the variable.
  probs <- calculate_exceedance_probabilities(model_results, variable, threshold)
  
  if (is.null(probs)) {
    warning("Unable to calculate probabilities for ", variable)
    return(ggplot() + theme_void())
  }
  
  tryCatch({
    ggplot() +
      geom_sf(
        data = map_data %>% left_join(probs, by = "LHD"),
        aes(fill = probability)
      ) +
      scale_fill_distiller(
        name = sprintf("Prob(RR > %.1f)", threshold),
        palette = "YlOrRd",
        direction = 1,
        limits = c(0, 1),
        na.value = "grey80"
      ) +
      theme_minimal() +
      labs(
        title = sprintf("%s: Probability of Exceeding RR=%.1f",
                        tools::toTitleCase(variable), threshold)
      ) +
      theme(legend.position = "right")
  }, error = function(e) {
    warning("Error in probability mapping: ", e$message)
    ggplot() + theme_void()
  })
}

# create_probability_page: Combines multiple probability maps into a single page.
# It iterates over exposures and several thresholds, then arranges the maps into a grid.
create_probability_page <- function(model_results, map_data, model_name) {
  plots <- list()
  
  for (var in c("temperature", "humidity", "rain")) {
    for (thresh in c(1.0, 1.2, 1.5)) {
      plot_key <- sprintf("%s_%.1f", var, thresh)
      plots[[plot_key]] <- plot_probability_map(
        model_results,
        map_data,
        var,
        thresh
      )
    }
  }
  
  do.call(gridExtra::grid.arrange,
          c(plots,
            ncol = 3,
            top = sprintf("%s: Probability Analysis", model_name)))
}

# ------------------------------------------------------------------------------
# 3) Risk Map Utility
# ------------------------------------------------------------------------------
# plot_time_series: Plots the standardized climate variables over time for each LHD.
# This helps to visually inspect temporal trends in temperature, humidity, and rain.
plot_time_series <- function(data) {
  # Standardize the variables for each LHD using z-scores.
  data_long <- data %>%
    select(date, LHD, temperature_2m, relative_humidity_2m, rain) %>%
    group_by(LHD) %>%
    mutate(across(c(temperature_2m, relative_humidity_2m, rain),
                  ~scale(.), .names = "scaled_{.col}")) %>%
    ungroup()
  
  ggplot(data_long, aes(x = date)) +
    geom_line(aes(y = scaled_temperature_2m, color = "Temperature")) +
    geom_line(aes(y = scaled_relative_humidity_2m, color = "Humidity")) +
    geom_line(aes(y = scaled_rain, color = "Rain")) +
    facet_wrap(~LHD) +
    theme_minimal() +
    scale_color_manual(
      values = c("Temperature" = "red", "Humidity" = "blue", "Rain" = "darkgreen")
    ) +
    labs(
      title = "Climate Variables by Local Health District",
      x = "Date",
      y = "Standardized Value",
      color = "Variable"
    ) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

# plot_lag_response: Plots the lagged effect estimates (with uncertainty) across
# different regions. This function visualizes how the effect changes with lag.
plot_lag_response <- function(effects, title) {
  ggplot(effects, aes(x = lag, y = effect, group = region)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = region), alpha = 0.2) +
    geom_line(aes(color = region)) +
    theme_minimal() +
    labs(
      title = title,
      x = "Lag (days)",
      y = "Effect"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6)
    )
}

# plot_risk_map: Visualizes the relative risk estimates on a spatial map.
# It requires that the risk data have a 'risk' (or 'relative_risk') column and joins
# the risk estimates with the spatial shapefile of LHDs.
plot_risk_map <- function(risks, map_data, title = "Risk Map") {
  # Ensure that the risks dataframe has the required columns.
  if (!"LHD" %in% names(risks)) {
    stop("risks dataframe must contain 'LHD' column")
  }
  
  # Allow both 'relative_risk' or 'risk' to be used as the risk measure.
  if ("relative_risk" %in% names(risks)) {
    risks <- risks %>% rename(risk = relative_risk)
  } else if (!"risk" %in% names(risks)) {
    stop("risks dataframe must contain either 'risk' or 'relative_risk' column")
  }
  if (!"LHD" %in% names(map_data)) {
    stop("map_data must contain 'LHD' column")
  }
  
  tryCatch({
    # Join the risk data with the spatial map data.
    plot_data <- map_data %>% 
      left_join(risks, by = "LHD")
    
    # Verify that the 'risk' column is present after the join.
    if (!"risk" %in% names(plot_data)) {
      stop("Join operation failed to preserve 'risk' column")
    }
    
    ggplot() +
      geom_sf(
        data = plot_data,
        aes(fill = risk)
      ) +
      scale_fill_distiller(
        name = "Relative Risk",
        palette = "YlOrRd",
        direction = 1,
        limits = c(0, 2),
        na.value = "grey80"
      ) +
      theme_minimal() +
      labs(title = title) +
      theme(legend.position = "right")
  }, error = function(e) {
    warning("Error in risk mapping: ", e$message)
    ggplot() + theme_void()
  })
}

# ------------------------------------------------------------------------------
# Effect Extraction Functions
# ------------------------------------------------------------------------------
# extract_relative_risks: Extracts the average relative risk (RR) from the MCMC 
# model results for each region and returns a data frame mapping LHD names to RR.
extract_relative_risks <- function(model_results, variable, n_regions = 15) {
  # Convert MCMC samples to a matrix.
  samples_matrix <- do.call(rbind, lapply(model_results, as.matrix))
  
  # Define the regex pattern for the specified variable.
  pattern <- switch(variable,
                    "temperature" = "^b_temp\\[",
                    "humidity" = "^b_hum\\[",
                    "rain" = "^b_rain\\[")
  
  # Calculate the mean effect for each region.
  rr_by_region <- numeric(n_regions)
  
  for(i in 1:n_regions) {
    region_pattern <- sprintf("%s%d,", pattern, i)
    region_cols <- grep(region_pattern, colnames(samples_matrix), value = TRUE)
    
    if(length(region_cols) > 0) {
      region_effects <- samples_matrix[, region_cols]
      rr_by_region[i] <- exp(mean(rowMeans(region_effects)))
    } else {
      rr_by_region[i] <- NA
    }
  }
  
  data.frame(
    LHD = bridge_table$shape_LHD,
    relative_risk = rr_by_region
  )
}

# extract_lag_effects: Extracts lagged effects for a specified variable from the
# MCMC model results for each region. Returns a combined data frame with effect
# estimates and corresponding confidence intervals for different lags.
extract_lag_effects <- function(model_results, variable, n_regions = 15) {
  # Convert MCMC samples to a matrix.
  samples_matrix <- do.call(rbind, lapply(model_results, as.matrix))
  
  pattern <- switch(variable,
                    "temperature" = "^b_temp\\[",
                    "humidity" = "^b_hum\\[",
                    "rain" = "^b_rain\\[")
  
  effects_list <- list()
  
  # Loop over regions to extract lag effects.
  for(i in 1:n_regions) {
    region_pattern <- sprintf("%s%d,", pattern, i)
    region_cols <- grep(region_pattern, colnames(samples_matrix), value = TRUE)
    
    if(length(region_cols) > 0) {
      effects <- samples_matrix[, region_cols]
      means <- colMeans(effects)
      cis <- apply(effects, 2, quantile, probs = c(0.025, 0.975))
      
      effects_list[[i]] <- data.frame(
        region = bridge_table$shape_LHD[i],
        lag = 0:(ncol(effects) - 1),
        effect = means,
        lower = cis[1, ],
        upper = cis[2, ]
      )
    }
  }
  
  do.call(rbind, effects_list)
}

# ------------------------------------------------------------------------------
# 4) Main Enhanced Report Generation
# ------------------------------------------------------------------------------
# create_enhanced_report: This function orchestrates the full analysis report.
# It loads the necessary data, validates inputs, creates various plots (time series,
# risk maps, probability maps), and arranges them into a PDF report.
create_enhanced_report <- function(output_dir = "Visualisations") {
  # Ensure the output directory exists; create it if not.
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Loading data...\n")
  data <- load_all_data()  # Load all required data components
  
  tryCatch({
    validate_data(data)  # Validate that all required components are present
    
    # Set the output PDF path and open the PDF device.
    pdf_path <- file.path(output_dir, "enhanced_climate_analysis.pdf")
    pdf(pdf_path, width = 15, height = 10)
    on.exit(dev.off())  # Ensure the PDF device is closed on exit or error
    
    cat("Creating overview...\n")
    # Plot time series of climate variables if available.
    if (!is.null(data$data_ts)) {
      print(plot_time_series(data$data_ts))
    }
    
    # Process each model's results to generate risk maps and probability maps.
    for (model_name in names(data$model_results)) {
      cat(sprintf("\nProcessing %s...\n", model_name))
      
      # 1) Generate Risk Maps
      plots <- list()
      for (var in c("temperature", "humidity", "rain")) {
        cat(sprintf("  Processing %s...\n", var))
        
        # Extract relative risk estimates from the model.
        risks <- extract_relative_risks(
          data$model_results[[model_name]], 
          var
        )
        
        # Print the structure of risk data for debugging.
        cat("  Risk data structure:\n")
        print(str(risks))
        
        # Create a risk map for the variable.
        plots[[var]] <- plot_risk_map(
          risks,
          data$nsw_map,
          sprintf("%s Effects", tools::toTitleCase(var))
        )
      }
      
      # Arrange the risk maps into a grid layout.
      grid.arrange(
        grobs = plots,
        ncol = 3,
        top = sprintf("%s: Risk Analysis", model_name)
      )
      
      # 2) Generate Probability Maps
      create_probability_page(
        data$model_results[[model_name]],
        data$nsw_map,
        model_name
      )
    }
    
    cat("\nAnalysis complete. Results saved in:", pdf_path, "\n")
    
  }, error = function(e) {
    cat("\nError in report generation:", e$message, "\n")
    if (dev.cur() > 1) dev.off()  # Ensure PDF is closed if an error occurs
  })
}

# ------------------------------------------------------------------------------
# 5) Run the Analysis
# ------------------------------------------------------------------------------
# The final step: call the main function to execute the enhanced report generation.
create_enhanced_report()
