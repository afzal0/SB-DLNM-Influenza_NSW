# Complete Response Curves Generator
# --------------------------------------------------------------------------------
# This script generates response curves that depict the relative risk (RR) 
# associated with changes in environmental exposures (temperature, humidity, 
# and rain) using coefficients estimated from DLNM-based models.
#
# The response curves are produced by:
#   1. Creating a crossbasis for the exposure variable, which captures the 
#      non-linear and lagged relationship.
#   2. Extracting model coefficients from the MCMC samples.
#   3. Predicting the exposure–response effect over a grid of values.
#   4. Plotting the estimated RR with its 95% confidence intervals.
#
# The results from multiple models (e.g., different designs or pooling approaches)
# are combined into a single panel for visual comparison.
# --------------------------------------------------------------------------------

# Load required libraries for plotting, data manipulation, and DLNM operations.
library(ggplot2)     # For creating plots
library(dplyr)       # For data manipulation
library(dlnm)        # For creating crossbasis functions (DLNM)
library(splines)     # For spline functions in the crossbasis
library(patchwork)   # For combining multiple ggplot2 plots into panels
library(viridis)     # For color scales in plots
library(gridExtra)   # For arranging grid-based plots
library(cowplot)     # For additional plotting utilities
library(scales)      # For scaling functions in ggplot2

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------

# safe_exp: A helper function to compute the exponential safely.
# Limits the input to the range [-50, 50] to avoid numerical overflow.
# Any resulting NA values are replaced with 1 (indicating a neutral effect).
safe_exp <- function(x) {
  result <- exp(pmin(pmax(x, -50), 50))  # Bound x to prevent overflow
  result[is.na(result)] <- 1             # Replace NA's with 1 (neutral effect)
  result
}

# create_crossbasis: Constructs a crossbasis object for a given exposure variable.
# A crossbasis captures both the non-linear exposure-response and lag-response relationships.
# If no knots are provided, they are set to the 10th, 50th, and 90th percentiles of the data.
create_crossbasis <- function(var_data, var_type, knots = NULL) {
  if(is.null(knots)) {
    knots <- quantile(var_data, c(0.10, 0.50, 0.90), na.rm = TRUE)
  }
  
  # Create a crossbasis using a natural spline for both the exposure and lag dimensions.
  # Here, a lag of 1 is assumed for all variables.
  cb <- crossbasis(
    var_data,
    lag = 1,
    argvar = list(fun = "ns", knots = knots),
    arglag = list(fun = "ns")
  )
  
  return(cb)
}

#------------------------------------------------------------------------------
# Data Processing Functions
#------------------------------------------------------------------------------

# extract_coefficients: Extracts model coefficients corresponding to a specific
# exposure variable from the MCMC samples of a fitted model.
# 'samples' are the MCMC samples, 'model_name' is a label (for information),
# 'variable' is one of "temperature", "humidity", or "rain", and 'region' (default 1)
# selects the region-specific coefficients.
extract_coefficients <- function(samples, model_name, variable, region = 1) {
  # Convert MCMC samples to a matrix if they are in mcmc.list format.
  if(inherits(samples, "mcmc.list")) {
    samples_matrix <- do.call(rbind, lapply(samples, as.matrix))
  } else {
    samples_matrix <- as.matrix(samples)
  }
  
  # Define a regex pattern based on the variable and region.
  pattern <- switch(variable,
                    "temperature" = sprintf("^b_temp\\[%d,", region),
                    "humidity" = sprintf("^b_hum\\[%d,", region),
                    "rain" = sprintf("^b_rain\\[%d,", region))
  
  # Extract coefficient columns matching the pattern.
  coef_cols <- grep(pattern, colnames(samples_matrix), value = TRUE)
  if(length(coef_cols) == 0) {
    return(NULL)
  }
  
  samples_matrix[, coef_cols]
}

# calculate_predictions: Calculates predicted relative risks using the model
# coefficients and the crossbasis representation of an exposure.
# A prediction grid of 100 points is created over the range of the exposure variable.
# The function computes the dot product between each MCMC sample’s coefficients and the
# prediction basis, applies safe exponentiation to convert log-RR to RR, and calculates
# the mean effect and 95% confidence intervals.
calculate_predictions <- function(coefs, var_data, cb) {
  n_points <- 100  # Number of points in the prediction grid
  pred_range <- range(var_data, na.rm = TRUE)
  pred_values <- seq(pred_range[1], pred_range[2], length.out = n_points)
  
  # Create a prediction basis for the grid values.
  pred_basis <- create_crossbasis(pred_values, "prediction")
  
  # Initialize matrix to store predictions from each MCMC sample.
  n_samples <- nrow(coefs)
  pred_matrix <- matrix(nrow = n_samples, ncol = n_points)
  
  # For each sample, compute predicted log-relative risk using the crossbasis.
  for(i in 1:n_samples) {
    pred_matrix[i,] <- coefs[i,] %*% t(pred_basis)
  }
  
  # Compute the mean predicted relative risk and 95% confidence intervals.
  mean_effect <- colMeans(safe_exp(pred_matrix), na.rm = TRUE)
  ci <- apply(safe_exp(pred_matrix), 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  
  list(
    x = pred_values,
    mean = mean_effect,
    lower = ci[1,],
    upper = ci[2,]
  )
}

#------------------------------------------------------------------------------
# Plotting Functions
#------------------------------------------------------------------------------

# plot_response_curve: Creates a ggplot2 object for the response curve.
# It plots the predicted mean relative risk along with a ribbon showing the 95%
# confidence interval, and includes a reference horizontal line at RR = 1.
plot_response_curve <- function(predictions, variable, title = NULL) {
  ggplot(data.frame(
    x = predictions$x,
    mean = predictions$mean,
    lower = predictions$lower,
    upper = predictions$upper
  ), aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey70") +
    geom_line(aes(y = mean), linewidth = 1, color = "blue") +
    geom_hline(yintercept = 1, linetype = "dashed") +  # Reference line indicating no effect
    theme_minimal() +
    labs(
      title = ifelse(is.null(title), paste(tools::toTitleCase(variable), "Effect"), title),
      x = tools::toTitleCase(variable),
      y = "Relative Risk"
    ) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

# create_model_panel: Combines response curve plots for the three exposures into a single panel.
# For each variable (temperature, humidity, and rain), it extracts the coefficients,
# calculates the predicted effect, generates the plot, and then combines them using patchwork.
create_model_panel <- function(data_ts, model_results, model_name) {
  plots <- list()
  
  for(variable in c("temperature", "humidity", "rain")) {
    # Retrieve the corresponding exposure data from the dataset.
    var_data <- switch(variable,
                       "temperature" = data_ts$temperature_2m,
                       "humidity" = data_ts$relative_humidity_2m,
                       "rain" = data_ts$rain)
    
    # Create a crossbasis for the variable.
    cb <- create_crossbasis(var_data, variable)
    
    # Extract the variable's coefficients from the model results.
    coefs <- extract_coefficients(model_results, model_name, variable)
    
    if(!is.null(coefs)) {
      # Calculate the predicted relative risk over a grid.
      preds <- calculate_predictions(coefs, var_data, cb)
      
      # Generate the response curve plot.
      plots[[variable]] <- plot_response_curve(preds, variable)
    }
  }
  
  # Combine the individual plots into a single panel.
  if(length(plots) > 0) {
    wrap_plots(plots, ncol = 3) +
      plot_annotation(
        title = model_name,
        theme = theme(
          plot.title = element_text(size = 12, hjust = 0.5, face = "bold")
        )
      )
  } else {
    # Return an empty plot if no valid data is available.
    ggplot() + 
      theme_void() + 
      annotate("text", x = 0.5, y = 0.5, label = paste("No valid data for", model_name))
  }
}

#------------------------------------------------------------------------------
# Main Function: Generate Response Curves and Save as PDF
#------------------------------------------------------------------------------

# create_response_curves: This is the main function that orchestrates the generation
# of response curves. It loads the necessary time series data and model results,
# iterates over each model to create response curve panels, and saves them to a PDF.
create_response_curves <- function(output_file = "Visualisations/response_curves.pdf") {
  # Function to print status messages
  status <- function(msg) cat(sprintf("%s\n", msg))
  
  # Load the time series data (contains exposure and outcome variables)
  status("Loading time series data...")
  load("new_output/data_timeseries.RData")
  
  # Load model results (assumed to be saved for model1, model2, model3, and model4)
  status("Loading model results...")
  model_results <- list()
  for(i in 1:4) {
    tryCatch({
      load(sprintf("new_output_acc/model%d_results.RData", i))
      model_results[[sprintf("model%d", i)]] <- get(sprintf("model%d_samples", i))
      status(sprintf("Loaded model%d successfully", i))
    }, error = function(e) {
      warning(sprintf("Could not load model%d: %s", i, e$message))
    })
  }
  
  # Open a PDF device to save the generated response curve panels.
  pdf(output_file, width = 12, height = 4)
  
  # Generate response curves for each model.
  for(model_name in names(model_results)) {
    status(sprintf("\nProcessing %s...", model_name))
    
    tryCatch({
      # Create a panel of response curves for the current model.
      p <- create_model_panel(data_ts, model_results[[model_name]], model_name)
      print(p)  # Print the panel to the PDF device.
      status(sprintf("Created plots for %s", model_name))
    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", model_name, e$message))
    })
  }
  
  # Close the PDF device after all plots are saved.
  dev.off()
  
  status(sprintf("\nResponse curves saved to: %s", output_file))
}

#------------------------------------------------------------------------------
# Additional Statistical Summaries
#------------------------------------------------------------------------------

# calculate_statistical_summaries: Computes basic summary statistics for the model 
# coefficients related to each exposure variable. For each model and variable, it returns
# the number of parameters, the mean effect (after applying safe exponentiation),
# and the 95% confidence interval.
calculate_statistical_summaries <- function(data_ts, model_results) {
  summaries <- list()
  
  for(model_name in names(model_results)) {
    model_summary <- list()
    
    for(variable in c("temperature", "humidity", "rain")) {
      coefs <- extract_coefficients(model_results[[model_name]], model_name, variable)
      
      if(!is.null(coefs)) {
        model_summary[[variable]] <- list(
          n_parameters = ncol(coefs),
          mean_effect = mean(safe_exp(rowMeans(coefs))),
          ci = quantile(safe_exp(rowMeans(coefs)), c(0.025, 0.975))
        )
      }
    }
    
    summaries[[model_name]] <- model_summary
  }
  
  return(summaries)
}

#------------------------------------------------------------------------------
# Execute Analysis
#------------------------------------------------------------------------------

# Generate and save the response curves into a PDF file.
create_response_curves()

# Optionally, you can calculate and print statistical summaries for the model results.
# Uncomment the following lines to obtain numerical summaries:
# summaries <- calculate_statistical_summaries(data_ts, model_results)
# print(summaries)
