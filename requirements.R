# R Package Requirements for SB-DLNM Influenza Analysis
# 
# Install all required packages by running:
# source("requirements.R")

# Check R version
if (getRversion() < "4.0.0") {
  warning("R version 4.0.0 or higher is recommended for this analysis")
}

# Function to install packages if not already installed
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  }
}

# Core statistical modeling packages
core_packages <- c(
  "dlnm",        # Distributed Lag Non-Linear Models
  "splines",     # Regression spline functions
  "survival",    # Survival analysis
  "coda"         # MCMC diagnostics
)

# Spatial analysis packages
spatial_packages <- c(
  "sf",          # Simple features for spatial data
  "spdep",       # Spatial dependence analysis
  "ape"          # Analyses of phylogenetics and evolution
)

# Data manipulation packages
data_packages <- c(
  "dplyr",       # Data manipulation
  "tidyr",       # Data tidying
  "data.table",  # Fast data manipulation
  "lubridate",   # Date/time handling
  "zoo",         # Time series infrastructure
  "readr"        # Fast file reading
)

# Visualization packages
viz_packages <- c(
  "ggplot2",     # Advanced plotting
  "plotly",      # Interactive plots
  "viridis",     # Color scales
  "RColorBrewer",# Color palettes
  "patchwork",   # Combine plots
  "gridExtra",   # Grid graphics
  "corrplot"     # Correlation plots
)

# Time series analysis packages
ts_packages <- c(
  "tseries",     # Time series analysis
  "forecast",    # Forecasting functions
  "trend"        # Trend analysis
)

# Additional statistical packages
stats_packages <- c(
  "car",         # Companion to Applied Regression
  "broom",       # Tidy model outputs
  "MASS",        # Modern Applied Statistics
  "mgcv"         # GAM models
)

# Combine all packages
all_packages <- c(core_packages, spatial_packages, data_packages, 
                  viz_packages, ts_packages, stats_packages)

# Install missing packages
cat("Checking and installing required R packages...\n")
install_if_missing(all_packages)

# Load all packages to verify installation
cat("\nLoading packages to verify installation...\n")
package_status <- sapply(all_packages, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))
})

# Report status
if (all(package_status)) {
  cat("\nAll packages successfully installed and loaded!\n")
} else {
  failed_packages <- names(package_status)[!package_status]
  cat("\nThe following packages failed to load:\n")
  print(failed_packages)
  cat("\nPlease install them manually using:\n")
  cat(paste0("install.packages(c(", paste0('"', failed_packages, '"', collapse = ", "), "))\n"))
}

# Print session info for reproducibility
cat("\n=== Session Information ===\n")
print(sessionInfo())