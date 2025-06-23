
# -------------------------------------------------------------------------------------
# Research Objective:
# This script is designed to explore the relationship between meteorological factors 
# (temperature, relative humidity, and rain) and influenza cases. We hypothesize that 
# weather conditions may influence the timing and intensity of influenza outbreaks.
#
# Methodological Approach:
# - Data Cleaning and Preparation: Ensure that only complete, reliable data is used.
# - Distributed Lag Non-Linear Modeling (DLNM): To capture both immediate and delayed 
#   (lagged) effects of weather on influenza cases, while allowing for non-linear 
#   relationships using spline functions.
# - Time Series Analysis: Decompose the temporal trends and seasonality in the data.
# - Spatial Analysis: Incorporate spatial dependencies using Local Health District (LHD)
#   boundaries to examine regional differences.
# - Diagnostic Visualizations: Generate plots to inspect trends, seasonal patterns, and 
#   variability at both overall and LHD-specific levels.
# -------------------------------------------------------------------------------------

# ----------------------
# Load Required Libraries
# ----------------------
library(dlnm)           # For Distributed Lag Non-linear Models: Captures lagged effects.
library(splines)        # Provides spline functions to flexibly model non-linear trends.
library(survival)       # Used for survival analysis; relevant for time-to-event data.
library(lubridate)      # Simplifies date and time manipulation.
library(sf)             # Handles spatial data and shapefiles for mapping and spatial analysis.
library(ggplot2)        # Provides advanced plotting capabilities.
library(dplyr)          # Simplifies data manipulation and transformation.
library(tidyr)          # Helps in tidying data for analysis.
library(RColorBrewer)   # Supplies color palettes, especially useful for mapping.
library(spdep)          # Analyzes spatial dependence; essential for spatial modeling.

# Create directories to store output files such as plots and processed data.
dir.create("output", showWarnings = FALSE)
dir.create("output/lhd_plots", showWarnings = FALSE)

# -------------------------------------------------------------------------------------
# Data Preparation
# -------------------------------------------------------------------------------------
# Load the combined dataset containing temperature, humidity, rain, and influenza cases.
data <- read.csv("data/combined_temp_influenza_data.csv")

# Load the shapefile containing NSW Local Health District (LHD) boundaries.
# This spatial data is used to map health districts and incorporate spatial dependencies.
shapefile_nsw <- st_read("spatial_data/NSW_LHD_Boundaries.shp")

# Generate spatial neighbor information from the shapefile.
# - poly2nb: Identifies neighboring polygons (LHDs) based on contiguity.
# - nb2WB: Converts the neighbor list into a format usable for spatial models.
list_neig <- nb2WB(poly2nb(shapefile_nsw))
# Save the spatial neighbor list for future use in analyses that require spatial structure.
save(list_neig, file = "output/list_neighbours.RData")

# Data validation and cleaning:
# Define the required columns to ensure we have all necessary variables for analysis.
required_cols <- c("temperature_2m", "relative_humidity_2m", 
                   "rain", "Cases", "date",
                   "spatial_lhd", "lhd_code")  

# Filter the data to only include rows with complete information in these key columns.
data_ts <- data[complete.cases(data[, required_cols]), ]

# Convert the date column into a proper Date object and extract month and year.
# This is crucial for identifying seasonal patterns and long-term trends.
data_ts$date <- as.Date(data_ts$date)
data_ts$month <- month(data_ts$date)
data_ts$year <- year(data_ts$date)

# Order the data by date to ensure accurate time series analyses.
data_ts <- data_ts[order(data_ts$date), ]

# -------------------------------------------------------------------------------------
# Calculate Thresholds for Analysis
# -------------------------------------------------------------------------------------
# Define baseline threshold values for the environmental variables.
# These thresholds act as reference points for modeling and interpretation.
thresholds <- list(
  temperature = mean(data_ts$temperature_2m, na.rm = TRUE), # Mean temperature as a threshold
  humidity = 60, # A fixed reference value for humidity based on prior research
  rain = mean(data_ts$rain, na.rm = TRUE) # Mean rainfall as a threshold
)

# -------------------------------------------------------------------------------------
# Define DLNM (Distributed Lag Non-linear Model) Parameters
# -------------------------------------------------------------------------------------
# DLNM enables us to explore how weather variables affect influenza cases over time,
# accounting for both immediate and lagged responses and allowing non-linear relationships.
dlnm_var <- list(
  df_trend = 4,     # Degrees of freedom to flexibly capture long-term trends.
  df_seas = 3,      # Degrees of freedom for modeling seasonality (monthly patterns).
  max_lag = 1       # Maximum lag period to consider; here it is set to 1 time unit.
)

# -------------------------------------------------------------------------------------
# Create Overall Crossbasis Objects for DLNM
# -------------------------------------------------------------------------------------
# The crossbasis function constructs a two-dimensional space that captures:
# 1. The exposure-response relationship (how different levels of the variable affect cases).
# 2. The lag-response relationship (how the effects are distributed over time).
#
# Natural splines (ns) are used to flexibly model non-linear patterns.
# Knots are set at the 10th, 50th, and 90th percentiles to capture variability across the range.
create_crossbasis <- function(data, var_name, max_lag = dlnm_var$max_lag) {
  var_data <- data[[var_name]]  # Extract variable data from the dataset.
  var_knots <- quantile(var_data, c(0.10, 0.50, 0.90), na.rm = TRUE)  # Define spline knots.
  
  # Construct the crossbasis object that will be used in the DLNM.
  cb <- crossbasis(
    var_data,
    lag = max_lag, 
    argvar = list(fun = "ns", knots = var_knots),  # Exposure-response function.
    arglag = list(fun = "ns")  # Lag-response function.
  )
  return(cb)
}

# Create crossbasis objects for each environmental variable of interest.
cb_list <- list()
variables <- c("temperature_2m", "relative_humidity_2m", "rain")  
for(var in variables) {
  cb_list[[var]] <- create_crossbasis(data_ts, var)
}

# -------------------------------------------------------------------------------------
# Create Overall Time Components: Trend and Seasonality
# -------------------------------------------------------------------------------------
# To adjust for overall temporal trends and seasonal variations, we use natural splines.

# Create a time trend basis over the sequential time index.
# This helps to adjust for any long-term trends unrelated to weather effects.
time_index <- seq(nrow(data_ts))
trend_ts <- onebasis(time_index, fun = "ns", df = dlnm_var$df_trend)

# Create a seasonality basis using the month variable.
# This step is critical since influenza and weather variables exhibit strong seasonal patterns.
seas_ts <- onebasis(data_ts$month, fun = "ns", df = dlnm_var$df_seas)

# -------------------------------------------------------------------------------------
# Create Local Health District (LHD)-specific Analyses
# -------------------------------------------------------------------------------------
# Given spatial heterogeneity in epidemiology, it is important to analyze data for each LHD.
# This provides insights into regional differences and tailors public health responses.

# Calculate summary statistics for each LHD:
# - Mean and standard deviation for temperature, humidity, and rain.
# - Total and average monthly influenza cases.
lhd_stats <- data_ts %>%
  group_by(LHD) %>%
  summarise(
    mean_temp = mean(temperature_2m, na.rm = TRUE),
    sd_temp = sd(temperature_2m, na.rm = TRUE),
    mean_humidity = mean(relative_humidity_2m, na.rm = TRUE),
    sd_humidity = sd(relative_humidity_2m, na.rm = TRUE),
    mean_rain = mean(rain, na.rm = TRUE),
    sd_rain = sd(rain, na.rm = TRUE),
    total_cases = sum(Cases, na.rm = TRUE),
    avg_monthly_cases = mean(Cases, na.rm = TRUE)
  )

# -------------------------------------------------------------------------------------
# Compile Summary Information
# -------------------------------------------------------------------------------------
# Collecting summary information aids in verifying data integrity and understanding the 
# baseline conditions of our dataset. It also documents the dimensions of various components.
summary_info <- list(
  data_range = range(data_ts$date),                # Time span of the data
  n_observations = nrow(data_ts),                   # Total number of observations
  thresholds = thresholds,                          # Reference thresholds for weather variables
  variable_summaries = lapply(data_ts[variables], summary),  # Summary stats for each variable
  missing_data = colSums(is.na(data_ts[variables])),  # Count of missing values per variable
  dimensions = list(
    trend = dim(trend_ts),                          # Dimensions of the trend basis
    seasonality = dim(seas_ts),                     # Dimensions of the seasonality basis
    crossbasis = lapply(cb_list, dim)               # Dimensions for each crossbasis object
  ),
  lhd_summaries = lhd_stats                         # Summary statistics for each LHD
)

# -------------------------------------------------------------------------------------
# Create Diagnostic Plots for Data Exploration
# -------------------------------------------------------------------------------------
# Visual diagnostics are vital to confirm that the data behaves as expected and to spot 
# potential issues such as outliers, trends, or anomalies.

# Function to create a time series plot for a given variable.
# Time series plots help in visualizing overall trends and potential lags in the data.
create_ts_plot <- function(data, var_name, ylab) {
  plot(data$date, data[[var_name]], 
       type = "l", 
       main = paste(ylab, "Time Series"),
       xlab = "Date",
       ylab = ylab,
       ylim = range(data[[var_name]], na.rm = TRUE))
}

# Function to create monthly boxplots for a given variable.
# Boxplots provide insights into the variability and distribution of data within each month.
create_monthly_boxplot <- function(data, var_name, ylab) {
  boxplot(data[[var_name]] ~ data$month,
          main = paste("Monthly", ylab),
          xlab = "Month",
          ylab = ylab,
          outline = TRUE,
          range = 1.5)
}

# -------------------------------------------------------------------------------------
# Generate Overall Diagnostic Plots
# -------------------------------------------------------------------------------------
# The following code creates diagnostic plots for the entire dataset to inspect trends 
# and seasonal patterns across the variables.
pdf("output/timeseries_diagnostics.pdf", width = 12, height = 8)

# Set up a multi-panel plot (2 rows x 2 columns) for the time series plots.
par(mfrow = c(2,2), mar = c(4,4,3,1))
create_ts_plot(data_ts, "temperature_2m", "Temperature (°C)")
create_ts_plot(data_ts, "relative_humidity_2m", "Relative Humidity (%)")
create_ts_plot(data_ts, "rain", "Rain (mm)")

# Reset layout for monthly boxplots.
par(mfrow = c(2,2), mar = c(4,4,3,1))
create_monthly_boxplot(data_ts, "temperature_2m", "Temperature (°C)")
create_monthly_boxplot(data_ts, "relative_humidity_2m", "Relative Humidity (%)")
create_monthly_boxplot(data_ts, "rain", "Rain (mm)")
dev.off()

# -------------------------------------------------------------------------------------
# Generate LHD-specific Diagnostic Plots
# -------------------------------------------------------------------------------------
# In addition to overall diagnostics, it is important to visualize data for each LHD.
# This helps identify any regional differences in weather patterns and influenza cases.
for(lhd in unique(data_ts$LHD)) {
  # Subset the data to only include observations for the current LHD.
  lhd_data <- subset(data_ts, LHD == lhd)
  
  # Define a PDF file name for the current LHD, replacing spaces with underscores.
  pdf(file = paste0("output/lhd_plots/", gsub(" ", "_", lhd), "_timeseries.pdf"),
      width = 12, height = 8)
  
  # Create time series plots for temperature, humidity, and rain for the LHD.
  par(mfrow = c(2,2), mar = c(4,4,3,1))
  create_ts_plot(lhd_data, "temperature_2m", paste(lhd, "Temperature (°C)"))
  create_ts_plot(lhd_data, "relative_humidity_2m", paste(lhd, "Humidity (%)"))
  create_ts_plot(lhd_data, "rain", paste(lhd, "Rain (mm)"))
  
  # Create monthly boxplots to visualize the seasonal variability in the LHD data.
  par(mfrow = c(2,2), mar = c(4,4,3,1))
  create_monthly_boxplot(lhd_data, "temperature_2m", paste(lhd, "Temperature"))
  create_monthly_boxplot(lhd_data, "relative_humidity_2m", paste(lhd, "Humidity"))
  create_monthly_boxplot(lhd_data, "rain", paste(lhd, "Rain"))
  dev.off()
}

# -------------------------------------------------------------------------------------
# Save Processed Objects for Reproducibility and Future Analysis
# -------------------------------------------------------------------------------------
# Saving key objects ensures that the analysis is reproducible and can be easily shared 
# or extended by other researchers.
save(data_ts, file = "output/data_timeseries.RData")          # Save the cleaned time series data.
save(cb_list, file = "output/crossbasis_timeseries.RData")       # Save the crossbasis objects.
save(trend_ts, file = "output/trend_timeseries.RData")           # Save the trend component.
save(seas_ts, file = "output/seasonality_timeseries.RData")      # Save the seasonality component.
save(summary_info, file = "output/timeseries_summary.RData")     # Save overall summary information.
save(lhd_stats, file = "output/lhd_statistics.RData")            # Save LHD-specific statistics.

# -------------------------------------------------------------------------------------
# Final Messages: Informing the User of Successful Execution
# -------------------------------------------------------------------------------------
# The following messages provide a quick summary of the analysis outputs,
# including thresholds used, the data range, number of observations, and the number
# of Local Health Districts processed.
cat("\nTime series components have been created and saved to output/\n")
cat("\nSummary of thresholds:\n")
print(thresholds)
cat("\nData range:", format(range(data_ts$date), "%Y-%m-%d"), "\n")
cat("Number of observations:", nrow(data_ts), "\n")
cat("\nAnalyses completed for", length(unique(data_ts$LHD)), "Local Health Districts\n")
