###############################################################################
#    LOAD MODEL RESULTS & CALCULATE DIC, BIC, QAIC FROM 'deviance' SAMPLES    #
###############################################################################
library(coda)  # Load the coda package to work with MCMC output

###############################################################################
# 1) Specify the paths to your saved model results, each containing MCMC draws
###############################################################################
model_files <- list(
  model1 = "new_output_acc/model1_results.RData",
  model2 = "new_output_acc/model2_results.RData",
  model3 = "new_output_acc/model3_results.RData",
  model4 = "new_output_acc/model4_results.RData"
)

# Define approximate parameter counts (K) and number of data points (n) for each model.
# These values are used to calculate model selection criteria (e.g., BIC and QAIC) that penalize
# model complexity relative to the sample size.
#
# Explanation of the values:
#
# Model 1 (Independent B-DLNM, Case-Crossover):
#   - K = 50:
#       This model includes stratum-specific intercepts and region-specific coefficients for each
#       climate variable (e.g., temperature, humidity, and rain). The effective parameter count
#       (including these intercepts and coefficients) is approximated as 50.
#   - n = 3405:
#       This is the number of observations used in the case-crossover design after removing missing values.
#
# Model 2 (Independent B-DLNM, Time-Series):
#   - K = 120:
#       In addition to the climate effect parameters as in Model 1, this model includes extra parameters
#       to capture overall trend and seasonal effects (such as region-specific trend coefficients and seasonal terms).
#       This increased complexity results in an approximate parameter count of 120.
#   - n = 4000:
#       The time-series design typically uses a larger dataset, here approximated as 4000 observations.
#
# Model 3 (SB-DLNM, Case-Crossover, Spatial Pooling):
#   - K = 60:
#       This sb-dlnm model still includes region-specific climate effects, but the partial pooling
#       via hyperparameters (means and scales) effectively "shrink" these estimates. This results in a
#       moderately increased effective parameter count, approximated here as 60.
#   - n = 3405:
#       The same case-crossover dataset is used as in Model 1.
#
# Model 4 (SB-DLNM, Time-Series, Spatial Pooling):
#   - K = 130:
#       This is the most complex model, as it combines the hierarchical structure for climate effects with
#       additional trend and seasonal parameters. The combined complexity leads to an approximate parameter count of 130.
#   - n = 4000:
#       As with Model 2, the time-series design uses around 4000 observations.
#
model_info <- list(
  model1 = list(K = 50,  n = 3405),  
  model2 = list(K = 120, n = 4000),
  model3 = list(K = 60,  n = 3405),
  model4 = list(K = 130, n = 4000)
)

###############################################################################
# 3) Overdispersion factor (c_hat) for QAIC calculation
#    Set c_hat = 1 if no overdispersion is suspected.
###############################################################################
c_hat <- 1

###############################################################################
# 4) Helper function to load the MCMC samples from a .RData file.
#    This function expects that the file contains an object ending with "_samples".
###############################################################################
load_coda_samples <- function(file_path) {
  # Check if the file exists; if not, print a message and return NULL.
  if (!file.exists(file_path)) {
    cat("File not found:", file_path, "\n")
    return(NULL)
  }
  
  # Create a new environment and load the RData file into it.
  env <- new.env()
  load(file_path, envir = env)
  
  # Retrieve all objects loaded in the environment.
  all_objs <- ls(env)
  
  # Identify the object that ends with "_samples" which contains the MCMC samples.
  samp_name <- all_objs[grepl("_samples$", all_objs)]
  
  # If exactly one object is found, return that object.
  if (length(samp_name) == 1) {
    return(env[[samp_name]])
  } else if (length(all_objs) == 1) {
    # If the file only has one object, assume it's the MCMC list.
    return(env[[all_objs]])
  } else {
    cat("Multiple objects found in", file_path, 
        "— unclear which is the MCMC samples.\n")
    return(NULL)
  }
}

###############################################################################
# 5) Function to compute DIC, BIC, and QAIC using the 'deviance' samples.
#
#    Arguments:
#      mcmc_list - an MCMC list (of class "mcmc.list") containing the posterior draws.
#      K         - approximate parameter count for the model.
#      n         - number of data points.
#      c_hat     - overdispersion factor for QAIC.
#
#    This function extracts the deviance samples, computes the mean deviance and 
#    variance, then uses these to calculate:
#       pD  : Effective number of parameters (approximated as var(dev)/2).
#       DIC : Deviance Information Criterion (mean deviance + pD).
#       BIC : Bayesian Information Criterion (mean deviance + K * log(n)).
#       QAIC: Quasi-Akaike Information Criterion (mean deviance + 2 * pD * c_hat).
###############################################################################
compute_model_criteria <- function(mcmc_list, K, n, c_hat = 1) {
  # Convert the mcmc.list to a matrix to work with all chains together.
  mat_draws <- as.matrix(mcmc_list)
  
  # Look for the column named "D", which holds the total deviance for each draw.
  dev_col <- which(colnames(mat_draws) == "D")
  if (length(dev_col) == 0) {
    cat("No 'deviance' column found in MCMC draws.\n")
    return(NULL)
  }
  
  # Extract the deviance samples from the identified column.
  deviance_values <- mat_draws[, dev_col]
  
  # Compute the mean and variance of the deviance samples.
  mean_dev <- mean(deviance_values)
  var_dev  <- var(deviance_values)
  
  # Effective number of parameters (pD) is approximated as half the variance.
  pD <- var_dev / 2
  
  # Calculate DIC as the sum of the mean deviance and the effective number of parameters.
  DIC <- mean_dev + pD
  
  # Approximate BIC: mean deviance plus a penalty based on parameter count and data size.
  BIC <- mean_dev + K * log(n)
  
  # QAIC adjusts for overdispersion using the factor c_hat.
  QAIC <- mean_dev + 2 * pD * c_hat
  
  # Return a list containing all the computed criteria along with mean deviance and pD.
  return(list(
    DIC      = DIC,
    BIC      = BIC,
    QAIC     = QAIC,
    mean_dev = mean_dev,
    pD       = pD
  ))
}

###############################################################################
# 6) Main script: Load each model's samples and compute DIC, BIC, QAIC.
###############################################################################
model_criteria <- list()  # List to store computed criteria for each model

# Loop over each model specified in model_files.
for (m in names(model_files)) {
  cat("\n--- Processing", m, "---\n")
  
  # 6a) Load the MCMC samples for the current model.
  mcmc_obj <- load_coda_samples(model_files[[m]])
  if (is.null(mcmc_obj)) {
    cat("No samples loaded for", m, "\n")
    next  # Skip to next model if samples cannot be loaded.
  }
  
  # 6b) Retrieve the parameter count (K) and number of data points (n) from model_info.
  if (!m %in% names(model_info)) {
    cat("No model_info specified for", m, "\n")
    next
  }
  K <- model_info[[m]]$K
  n <- model_info[[m]]$n
  
  # 6c) Compute the model criteria (DIC, BIC, QAIC) using the helper function.
  crit <- compute_model_criteria(mcmc_obj, K = K, n = n, c_hat = c_hat)
  if (is.null(crit)) {
    cat("Could not compute criteria for", m, "\n")
  } else {
    model_criteria[[m]] <- crit
    cat("DIC  =", crit$DIC,  "\n")
    cat("BIC  =", crit$BIC,  "\n")
    cat("QAIC =", crit$QAIC, "\n")
  }
}

###############################################################################
# 7) Print final model comparison results
###############################################################################
cat("\n--- Final Model Comparison Results ---\n")
print(model_criteria)
cat("\nDone.\n")
