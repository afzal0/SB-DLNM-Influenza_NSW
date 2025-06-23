###############################################################################
#                       MODEL FITTING SCRIPTS                                 #
###############################################################################
#
# OVERVIEW:
# This script fits four different models to explore the relationship between
# meteorological variables and influenza cases. We use a Distributed Lag 
# Non-Linear Modeling (DLNM) framework and account for lagged, non-linear effects.
#
# Two main study designs are implemented:
# 1. Case-Crossover: Each observation is compared with control periods within
#    the same stratum (defined by region, year, and month) to control for time-
#    invariant confounders.
# 2. Time-Series: Overall Poisson regression on the time series data with 
#    adjustments for trends and seasonality.
#
# The models are specified in JAGS and include:
# - Independent B-DLNM (Models 1 & 2)
# - Spatially/Hierarchically Pooled DLNM (SB-DLNM, Models 3 & 4)
#
# MATHEMATICAL FORMULATION:
#
# For Models 1 and 2 (Independent B-DLNM), the likelihood for each observation i is:
#
#   y_i ~ Poisson(λ_i)
#
# with the log-link defined as:
#
#   log(λ_i) = α[s_i] + ∑_{j=1}^{J} β_temp[r_i, j] * w_temp[i, j] +
#                       ∑_{j=1}^{J} β_hum[r_i, j]  * w_hum[i, j] +
#                       ∑_{j=1}^{J} β_rain[r_i, j] * w_rain[i, j] +
#                       ... [plus terms for trend/seasonality in Model 2]
#
# where:
#   - α[s_i] is the stratum-specific intercept.
#   - r_i indexes the region of observation i.
#   - w_temp, w_hum, and w_rain are the crossbasis matrices capturing the 
#     non-linear exposure-lag relationship, with J coefficients per variable.
#
# In the hierarchical models (Models 3 and 4), the region-specific coefficients 
# are modeled with partial pooling:
#
#   β[r, j] = μ[j] + σ[j] * θ[r, j]   with   θ[r, j] ~ dnorm(0, τ)
#
# This formulation allows for borrowing of strength across regions.
#
# NA removal is performed before model fitting to ensure that the crossbasis
# matrices have complete data.
#

# ----------------------
# Load required libraries
# ----------------------
library(rjags)   # Interface to JAGS for Bayesian MCMC sampling
library(coda)    # Tools for analyzing MCMC output
library(dlnm)    # For Distributed Lag Non-Linear Models (DLNM)
library(splines) # Provides functions for spline-based modeling
library(survival) # Includes survival analysis functions (often used in DLNM)
library(lubridate) # For handling date/time operations (e.g., extracting month/year)
library(dplyr)     # Data manipulation and transformation

# ----------------------
# 1) Load Data Objects
# ----------------------
# Load the pre-processed data objects. These were created in earlier steps.
load("new_output/data_timeseries.RData")         # Loads 'data_ts'
load("new_output/crossbasis_timeseries.RData")     # Loads 'cb_list'
load("new_output/trend_timeseries.RData")          # Loads 'trend_ts'
load("new_output/seasonality_timeseries.RData")    # Loads 'seas_ts'

# Optionally load spatial neighbor information if needed in spatial models.
load("output/list_neighbours.RData")               # List of spatial neighbors

# Define which models to run
execute_model <- list(
  "model1" = TRUE, # Independent B-DLNM with case-crossover design
  "model2" = TRUE, # Independent B-DLNM with time-series design
  "model3" = TRUE, # SB-DLNM with case-crossover design (partial/hierarchical)
  "model4" = TRUE  # SB-DLNM with time-series design (partial/hierarchical)
)

# Set MCMC parameters for each model
n.iter <- list(
  "model1" = 5000,
  "model2" = 5000,
  "model3" = 5000,
  "model4" = 5000
)
n.chains <- 3

# -------------------------------------------------------------------------
# 2) Prepare Data for Case-Crossover & Time-Series Designs
# -------------------------------------------------------------------------
# For the case-crossover design, we create a new dataset (data_cco) from data_ts.
# Strata are defined as combinations of region (lhd_code), year, and month.
# This controls for time-invariant confounders within each stratum.

# (A) CASE-CROSSOVER data preparation:
data_cco <- data_ts

# Create 'strata' by combining lhd_code, year, and month (formatted as two digits)
data_cco$strata <- paste(
  data_cco$lhd_code,
  year(data_cco$date),
  formatC(month(data_cco$date), width = 2, flag = "0"),
  sep = ":"
)

# (Optionally) Remove strata where the total number of cases is zero
keep <- sapply(split(data_cco$Cases, data_cco$strata), sum)
keep <- data_cco$strata %in% names(keep[keep != 0])
data_cco <- data_cco[keep, ]

# Create a numeric region index for both datasets (used in modeling)
data_cco$region_index <- as.numeric(factor(data_cco$lhd_code))
data_ts$region_index  <- as.numeric(factor(data_ts$lhd_code))
max_region <- max(data_ts$region_index)

# -------------------------------------------------------------------------
# 3) Create Crossbasis Function
# -------------------------------------------------------------------------
# The function below creates a crossbasis matrix for a given variable.
# It uses natural splines (ns) with knots at the 10th, 50th, and 90th percentiles.
#
# The crossbasis matrix (w) captures the exposure-response relationship as well as
# the lagged effects. Mathematically, for a variable X:
#
#   w[i, j] = f(X[i], lag=j)
#
# where f() is modeled using a natural spline.
#
create_crossbasis_matrix <- function(data, var_name, lag = 1) {
  var_knots <- quantile(data[[var_name]], c(0.10, 0.50, 0.90), na.rm = TRUE)
  cb <- crossbasis(
    x   = data[[var_name]],
    lag = lag,
    argvar = list(fun = "ns", knots = var_knots),  # Non-linear exposure-response
    arglag = list(fun = "ns")                       # Non-linear lag-response
  )
  as.matrix(cb)
}

###############################################################################
# MODEL 1: Independent B-DLNM, Case-Crossover (NA-removal included)
###############################################################################
if (execute_model$model1) {
  cat("\n-------------------\nFitting Model 1...\n-------------------\n")
  
  # 1) Create crossbasis matrices for temperature, humidity, and rain
  w_temp <- create_crossbasis_matrix(data_cco, "temperature_2m", lag = 1)
  w_hum  <- create_crossbasis_matrix(data_cco, "relative_humidity_2m", lag = 1)
  w_rain <- create_crossbasis_matrix(data_cco, "rain", lag = 1)
  
  # 2) Identify rows with any missing (NA) values in the crossbasis matrices
  has_na <- apply(w_temp, 1, anyNA) | apply(w_hum, 1, anyNA) | apply(w_rain, 1, anyNA)
  
  # 3) Remove rows with NA values from both the data and crossbasis matrices
  data_cco_clean <- data_cco[!has_na, ]
  w_temp_clean   <- w_temp[!has_na, , drop = FALSE]
  w_hum_clean    <- w_hum[!has_na, , drop = FALSE]
  w_rain_clean   <- w_rain[!has_na, , drop = FALSE]
  
  cat("Model 1 rows removed due to NA:", sum(has_na), "\n")
  cat("Remaining rows:", nrow(data_cco_clean), "\n")
  
  # 4) Prepare the data list for JAGS
  # This list passes all necessary data objects to the JAGS model.
  model1_data <- list(
    y            = as.numeric(data_cco_clean$Cases),
    strata       = as.numeric(factor(data_cco_clean$strata)),
    region_index = as.numeric(data_cco_clean$region_index),
    w_temp       = w_temp_clean,
    w_hum        = w_hum_clean,
    w_rain       = w_rain_clean,
    n_data       = nrow(data_cco_clean),
    n_strata     = length(unique(data_cco_clean$strata)),
    n_reg        = max(data_cco_clean$region_index),
    n_coef       = ncol(w_temp_clean)
  )
  
  # JAGS model specification for Model 1:
  # This model assumes:
  #   y_i ~ Poisson(λ_i)
  #   log(λ_i) = p.alpha[strata[i]] +
  #              inprod(b_temp[region_index[i], ], w_temp[i, ]) +
  #              inprod(b_hum[region_index[i], ],  w_hum[i, ]) +
  #              inprod(b_rain[region_index[i], ], w_rain[i, ])
  #
  # "inprod" calculates the dot product between the coefficient vector and the 
  # corresponding crossbasis vector for each observation.
  #
  # Deviance is calculated for each observation for use in model diagnostics (DIC).
  model1_jags <- "
  model {
    # Likelihood: Poisson model for count data
    for(i in 1:n_data) {
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- p.alpha[strata[i]] +
                        inprod(b_temp[region_index[i],1:n_coef], w_temp[i,1:n_coef]) +
                        inprod(b_hum[region_index[i],1:n_coef],  w_hum[i,1:n_coef]) +
                        inprod(b_rain[region_index[i],1:n_coef], w_rain[i,1:n_coef])
      
      # Calculate the log-likelihood and deviance for each observation
      log_like[i] <- y[i] * log(lambda[i]) - lambda[i] - logfact(y[i])
      dev[i] <- -2 * log_like[i]
    }
    
    # Summed deviance (D) for DIC calculations
    D <- sum(dev[1:n_data])
    pD <- 0  # Effective number of parameters (to be computed post-hoc)
    DIC <- 0 # Deviance Information Criterion (to be computed post-hoc)
    
    # Stratum-specific intercepts to account for case-crossover design
    for(s in 1:n_strata) {
      p.alpha[s] ~ dnorm(0, 0.001)
    }

    # Region-specific coefficients for each climate variable
    for(r in 1:n_reg) {
      for(j in 1:n_coef) {
        b_temp[r,j] ~ dnorm(0, 0.001)
        b_hum[r,j]  ~ dnorm(0, 0.001)
        b_rain[r,j] ~ dnorm(0, 0.001)
      }
    }
  }
  "
  
  # Initial values for the model parameters
  model1_inits <- function() {
    list(
      p.alpha = rep(0, model1_data$n_strata),
      b_temp  = matrix(0, model1_data$n_reg, model1_data$n_coef),
      b_hum   = matrix(0, model1_data$n_reg, model1_data$n_coef),
      b_rain  = matrix(0, model1_data$n_reg, model1_data$n_coef)
    )
  }
  
  # List of parameters to monitor
  model1_params <- c("p.alpha", "b_temp", "b_hum", "b_rain", "D", "log_like","dev")
  
  # 5) Fit the model in JAGS using MCMC sampling
  tryCatch({
    jags_model1 <- jags.model(
      textConnection(model1_jags),
      data  = model1_data,
      inits = model1_inits,
      n.chains = n.chains
    )
    
    # Burn-in period: discard the first half of iterations
    update(jags_model1, n.iter$model1 / 2)
    
    # Draw samples from the posterior distribution
    model1_samples <- coda.samples(
      jags_model1,
      variable.names = model1_params,
      n.iter = n.iter$model1 / 2
    )
    
    # Calculate DIC components using the posterior deviance
    D_posterior <- do.call(rbind, lapply(model1_samples, function(x) x[, "D"]))
    D_mean <- mean(D_posterior)
    D_hat <- 2 * D_mean - mean(apply(D_posterior, 1, sum))
    p_D <- D_mean - D_hat
    DIC <- D_hat + 2 * p_D
    
    model1_dic <- list(
      D_mean = D_mean,
      D_hat = D_hat,
      p_D = p_D,
      DIC = DIC
    )
    
    # Save the MCMC samples and DIC results
    save(model1_samples, model1_dic, file = "new_output_acc/model1_results.RData")
    cat("Model 1 completed and saved. DIC =", DIC, "\n")
    
  }, error = function(e) {
    cat("\nError in Model 1:", conditionMessage(e), "\n")
  })
}

###############################################################################
# MODEL 2: Independent B-DLNM, Time-Series (NA-removal included)
###############################################################################
if (execute_model$model2) {
  cat("\n-------------------\nFitting Model 2...\n-------------------\n")
  
  # 1) Create crossbasis matrices for the time-series data (data_ts)
  w_temp2 <- create_crossbasis_matrix(data_ts, "temperature_2m", lag = 1)
  w_hum2  <- create_crossbasis_matrix(data_ts, "relative_humidity_2m", lag = 1)
  w_rain2 <- create_crossbasis_matrix(data_ts, "rain", lag = 1)
  
  # 2) Identify rows with any missing values in the crossbasis matrices
  has_na2 <- apply(w_temp2, 1, anyNA) | apply(w_hum2, 1, anyNA) | apply(w_rain2, 1, anyNA)
  
  # 3) Filter out rows with NAs from the data and matrices
  data_ts_clean <- data_ts[!has_na2, ]
  w_temp2_clean <- w_temp2[!has_na2, , drop = FALSE]
  w_hum2_clean  <- w_hum2[!has_na2, , drop = FALSE]
  w_rain2_clean <- w_rain2[!has_na2, , drop = FALSE]
  
  cat("Model 2 rows removed due to NA:", sum(has_na2), "\n")
  cat("Remaining rows:", nrow(data_ts_clean), "\n")
  
  # 4) Prepare the JAGS data list for the time-series model.
  # In this model, we include additional terms for overall trend (t) and seasonality (s).
  model2_data <- list(
    y            = as.numeric(data_ts_clean$Cases),
    region_index = as.numeric(data_ts_clean$region_index),
    month        = as.numeric(data_ts_clean$month),
    year         = as.numeric(data_ts_clean$year - min(data_ts_clean$year) + 1),
    w_temp       = w_temp2_clean,
    w_hum        = w_hum2_clean,
    w_rain       = w_rain2_clean,
    t            = as.matrix(trend_ts),    # Trend component from trend_timeseries
    s            = as.matrix(seas_ts),     # Seasonality component from seasonality_timeseries
    n_data       = nrow(data_ts_clean),
    n_reg        = length(unique(data_ts_clean$region_index)),
    n_coef       = ncol(w_temp2_clean),
    n_trend      = ncol(trend_ts),
    n_seas       = ncol(seas_ts),
    n_year       = length(unique(data_ts_clean$year))
  )
  
  # JAGS model specification for Model 2:
  # The model structure is similar to Model 1 but includes additional terms:
  #   log(λ_i) = α[r, m] +
  #              inprod(b_temp[r, ], w_temp[i, ]) +
  #              inprod(b_hum[r, ], w_hum[i, ]) +
  #              inprod(b_rain[r, ], w_rain[i, ]) +
  #              inprod(gamma[r, ], t[i, ]) +
  #              inprod(delta[r, y_i, ], s[i, ])
  #
  # Here, α[r, m] represents a region- and month-specific intercept, gamma[r, ]
  # captures the trend effect, and delta[r, y, ] represents seasonal effects for each year.
  model2_jags <- "
  model {
    for(i in 1:n_data) {
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- alpha[region_index[i], month[i]] +
                        inprod(b_temp[region_index[i],1:n_coef], w_temp[i,1:n_coef]) +
                        inprod(b_hum[region_index[i],1:n_coef],  w_hum[i,1:n_coef]) +
                        inprod(b_rain[region_index[i],1:n_coef], w_rain[i,1:n_coef]) +
                        inprod(gamma[region_index[i],1:n_trend], t[i,1:n_trend]) +
                        inprod(delta[region_index[i], year[i],1:n_seas], s[i,1:n_seas])
      
      # Deviance calculation for each observation
      log_like[i] <- y[i] * log(lambda[i]) - lambda[i] - logfact(y[i])
      dev[i] <- -2 * log_like[i]
    }
    
    # Total deviance for DIC calculation
    D <- sum(dev[1:n_data])
    pD <- 0  # Effective number of parameters (to be computed)
    DIC <- 0 # Deviance Information Criterion (to be computed)

    # Priors for region-specific monthly intercepts
    for(r in 1:n_reg) {
      for(m in 1:12) {
        alpha[r, m] ~ dnorm(0, 0.001)
      }
      # Priors for the climate effects
      for(j in 1:n_coef) {
        b_temp[r, j] ~ dnorm(0, 0.001)
        b_hum[r, j]  ~ dnorm(0, 0.001)
        b_rain[r, j] ~ dnorm(0, 0.001)
      }
      # Priors for the trend effects
      for(k in 1:n_trend) {
        gamma[r, k] ~ dnorm(0, 0.001)
      }
      # Priors for the seasonal effects by year
      for(yv in 1:n_year) {
        for(ss in 1:n_seas) {
          delta[r, yv, ss] ~ dnorm(0, 0.001)
        }
      }
    }
  }
  "
  
  # Initial values for Model 2 parameters
  model2_inits <- function(){
    list(
      alpha   = matrix(0, model2_data$n_reg, 12),
      b_temp  = matrix(0, model2_data$n_reg, model2_data$n_coef),
      b_hum   = matrix(0, model2_data$n_reg, model2_data$n_coef),
      b_rain  = matrix(0, model2_data$n_reg, model2_data$n_coef),
      gamma   = matrix(0, model2_data$n_reg, model2_data$n_trend),
      delta   = array(0, dim = c(model2_data$n_reg, model2_data$n_year, model2_data$n_seas))
    )
  }
  
  model2_params <- c("alpha", "b_temp", "b_hum", "b_rain", "gamma", "delta", "D", "log_like","dev")
  
  # Fit Model 2 in JAGS
  tryCatch({
    jags_model2 <- jags.model(
      textConnection(model2_jags),
      data  = model2_data,
      inits = model2_inits,
      n.chains = n.chains
    )
    
    update(jags_model2, n.iter$model2 / 2)
    
    model2_samples <- coda.samples(
      jags_model2,
      variable.names = model2_params,
      n.iter = n.iter$model2 / 2
    )
    
    # Compute DIC components from the posterior samples
    D_posterior <- do.call(rbind, lapply(model2_samples, function(x) x[, "D"]))
    D_mean <- mean(D_posterior)
    D_hat <- 2 * D_mean - mean(apply(D_posterior, 1, sum))
    p_D <- D_mean - D_hat
    DIC <- D_hat + 2 * p_D
    
    model2_dic <- list(
      D_mean = D_mean,
      D_hat = D_hat,
      p_D = p_D,
      DIC = DIC
    )
    
    save(model2_samples, model2_dic, file = "new_output_acc/model2_results.RData")
    cat("Model 2 completed and saved. DIC =", DIC, "\n")
    
  }, error = function(e) {
    cat("\nError in Model 2:", conditionMessage(e), "\n")
  })
}

###############################################################################
# MODEL 3: SB-DLNM, Case-Crossover (Hierarchical) with NA removal
###############################################################################
if (execute_model$model3) {
  cat("\n-------------------\nFitting Model 3...\n-------------------\n")
  
  # 1) Create crossbasis matrices for the case-crossover data (data_cco)
  w_temp3 <- create_crossbasis_matrix(data_cco, "temperature_2m", lag = 1)
  w_hum3  <- create_crossbasis_matrix(data_cco, "relative_humidity_2m", lag = 1)
  w_rain3 <- create_crossbasis_matrix(data_cco, "rain", lag = 1)
  
  # 2) Identify rows with any NA in the crossbasis matrices
  has_na3 <- apply(w_temp3, 1, anyNA) | apply(w_hum3, 1, anyNA) | apply(w_rain3, 1, anyNA)
  
  # 3) Remove observations with NAs from both the dataset and the matrices
  data_cco_clean3 <- data_cco[!has_na3, ]
  w_temp3_clean   <- w_temp3[!has_na3, , drop = FALSE]
  w_hum3_clean    <- w_hum3[!has_na3, , drop = FALSE]
  w_rain3_clean   <- w_rain3[!has_na3, , drop = FALSE]
  
  cat("Model 3 rows removed due to NA:", sum(has_na3), "\n")
  cat("Remaining rows:", nrow(data_cco_clean3), "\n")
  
  # 4) Prepare the JAGS data for Model 3
  # In this hierarchical model, the region-specific climate effects are modeled with partial pooling.
  model3_data <- list(
    y            = as.numeric(data_cco_clean3$Cases),
    strata       = as.numeric(factor(data_cco_clean3$strata)),
    region_index = as.numeric(data_cco_clean3$region_index),
    w_temp       = w_temp3_clean,
    w_hum        = w_hum3_clean,
    w_rain       = w_rain3_clean,
    n_data       = nrow(data_cco_clean3),
    n_strata     = length(unique(data_cco_clean3$strata)),
    n_reg        = length(unique(data_cco_clean3$region_index)),
    n_coef       = ncol(w_temp3_clean)
  )
  
  # JAGS model specification for Model 3 (Hierarchical, case-crossover):
  # The model is similar to Model 1, but with hierarchical structure for climate effects:
  #
  #   b[r, j] = μ[j] + σ[j] * θ[r, j]   with   θ[r, j] ~ dnorm(0, τ)
  #
  # This allows for partial pooling of the coefficients across regions.
  model3_jags <- "
  model {
    # Likelihood for each observation
    for(i in 1:n_data) {
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- p.alpha[strata[i]] +
                        inprod(b_temp[region_index[i],], w_temp[i,]) +
                        inprod(b_hum[region_index[i],],  w_hum[i,]) +
                        inprod(b_rain[region_index[i],], w_rain[i,])
      
      # Deviance calculation
      log_like[i] <- y[i] * log(lambda[i]) - lambda[i] - logfact(y[i]) 
      dev[i] <- -2 * log_like[i]
    }
    
    # Total deviance for the model
    D <- sum(dev[1:n_data])
    pD <- 0  
    DIC <- 0 

    # Stratum-specific intercepts
    for(s in 1:n_strata) {
      p.alpha[s] ~ dnorm(0, 0.001)
    }

    # Hierarchical specification for region-specific climate effects
    for(r in 1:n_reg) {
      for(j in 1:n_coef) {
        b_temp[r,j] <- mu_temp[j] + sigma_temp[j] * theta_temp[r,j]
        theta_temp[r,j] ~ dnorm(0, tau_temp)

        b_hum[r,j]  <- mu_hum[j] + sigma_hum[j] * theta_hum[r,j]
        theta_hum[r,j] ~ dnorm(0, tau_hum)

        b_rain[r,j] <- mu_rain[j] + sigma_rain[j] * theta_rain[r,j]
        theta_rain[r,j] ~ dnorm(0, tau_rain)
      }
    }

    # Hyperpriors for the mean and scale of climate effects
    for(j in 1:n_coef) {
      mu_temp[j]  ~ dnorm(0, 0.001)
      mu_hum[j]   ~ dnorm(0, 0.001)
      mu_rain[j]  ~ dnorm(0, 0.001)

      sigma_temp[j]  ~ dunif(0,10)
      sigma_hum[j]   ~ dunif(0,10)
      sigma_rain[j]  ~ dunif(0,10)
    }

    # Priors for the precision of the random effects
    tau_temp ~ dgamma(0.1, 0.1)
    tau_hum  ~ dgamma(0.1, 0.1)
    tau_rain ~ dgamma(0.1, 0.1)
  }
  "
  
  # Initial values for hierarchical model parameters
  model3_inits <- function() {
    list(
      p.alpha    = rep(0, model3_data$n_strata),
      theta_temp = matrix(0, model3_data$n_reg, model3_data$n_coef),
      theta_hum  = matrix(0, model3_data$n_reg, model3_data$n_coef),
      theta_rain = matrix(0, model3_data$n_reg, model3_data$n_coef),
      mu_temp    = rep(0, model3_data$n_coef),
      mu_hum     = rep(0, model3_data$n_coef),
      mu_rain    = rep(0, model3_data$n_coef),
      sigma_temp = rep(1, model3_data$n_coef),
      sigma_hum  = rep(1, model3_data$n_coef),
      sigma_rain = rep(1, model3_data$n_coef),
      tau_temp   = 1,
      tau_hum    = 1,
      tau_rain   = 1
    )
  }
  
  model3_params <- c(
    "p.alpha",
    "b_temp", "b_hum", "b_rain",
    "mu_temp", "mu_hum", "mu_rain",
    "sigma_temp", "sigma_hum", "sigma_rain",
    "tau_temp", "tau_hum", "tau_rain",
    "D", "log_like","dev"
  )
  
  # Fit Model 3 in JAGS
  tryCatch({
    jags_model3 <- jags.model(
      textConnection(model3_jags),
      data  = model3_data,
      inits = model3_inits,
      n.chains = n.chains
    )
    
    update(jags_model3, n.iter$model3 / 2)
    model3_samples <- coda.samples(
      jags_model3,
      variable.names = model3_params,
      n.iter = n.iter$model3 / 2
    )
    
    # Calculate DIC components from the posterior samples
    D_posterior <- do.call(rbind, lapply(model3_samples, function(x) x[, "D"]))
    D_mean <- mean(D_posterior)
    D_hat <- 2 * D_mean - mean(apply(D_posterior, 1, sum))
    p_D <- D_mean - D_hat
    DIC <- D_hat + 2 * p_D
    
    model3_dic <- list(
      D_mean = D_mean,
      D_hat = D_hat,
      p_D = p_D,
      DIC = DIC
    )
    
    save(model3_samples, model3_dic, file = "new_output_acc/model3_results.RData")
    cat("Model 3 completed and saved. DIC =", DIC, "\n")
    
  }, error = function(e) {
    cat("\nError in Model 3:", conditionMessage(e), "\n")
  })
}

###############################################################################
# MODEL 4: SB-DLNM, Time-Series (Hierarchical) with NA removal
###############################################################################
if (execute_model$model4) {
  cat("\n-------------------\nFitting Model 4...\n-------------------\n")
  
  # 1) Create crossbasis matrices for the time-series data (data_ts)
  w_temp4 <- create_crossbasis_matrix(data_ts, "temperature_2m", lag = 1)
  w_hum4  <- create_crossbasis_matrix(data_ts, "relative_humidity_2m", lag = 1)
  w_rain4 <- create_crossbasis_matrix(data_ts, "rain", lag = 1)
  
  # 2) Identify rows with any missing values
  has_na4 <- apply(w_temp4, 1, anyNA) | apply(w_hum4, 1, anyNA) | apply(w_rain4, 1, anyNA)
  
  # 3) Filter out rows with NA values from both the dataset and the matrices
  data_ts_clean4 <- data_ts[!has_na4, ]
  w_temp4_clean  <- w_temp4[!has_na4, , drop = FALSE]
  w_hum4_clean   <- w_hum4[!has_na4, , drop = FALSE]
  w_rain4_clean  <- w_rain4[!has_na4, , drop = FALSE]
  
  cat("Model 4 rows removed due to NA:", sum(has_na4), "\n")
  cat("Remaining rows:", nrow(data_ts_clean4), "\n")
  
  # 4) Prepare the JAGS data list for the hierarchical time-series model.
  # This model includes both climate effects and additional terms for trend and seasonality.
  model4_data <- list(
    y            = as.numeric(data_ts_clean4$Cases),
    region_index = as.numeric(data_ts_clean4$region_index),
    month        = as.numeric(data_ts_clean4$month),
    year         = as.numeric(data_ts_clean4$year - min(data_ts_clean4$year) + 1),
    w_temp       = w_temp4_clean,
    w_hum        = w_hum4_clean,
    w_rain       = w_rain4_clean,
    t            = as.matrix(trend_ts), # Trend component
    s            = as.matrix(seas_ts),  # Seasonality component
    n_data       = nrow(data_ts_clean4),
    n_reg        = length(unique(data_ts_clean4$region_index)),
    n_coef       = ncol(w_temp4_clean),
    n_trend      = ncol(trend_ts),
    n_seas       = ncol(seas_ts),
    n_year       = length(unique(data_ts_clean4$year))
  )
  
  # JAGS model specification for Model 4 (Hierarchical, time-series):
  # The model is similar to Model 2 but introduces hierarchical (partial pooling)
  # for the climate effects as in Model 3.
  #
  #   log(λ_i) = α[r, m] +
  #              inprod(b_temp[r, ], w_temp[i, ]) +
  #              inprod(b_hum[r, ],  w_hum[i, ]) +
  #              inprod(b_rain[r, ], w_rain[i, ]) +
  #              inprod(gamma[r, ], t[i, ]) +
  #              inprod(delta[r, y_i, ], s[i,])
  #
  # with:
  #   b[r, j] = μ[j] + σ[j] * θ[r, j]  and  θ[r, j] ~ dnorm(0, τ)
  model4_jags <- "
  model {
    for(i in 1:n_data) {
      y[i] ~ dpois(lambda[i])
      log(lambda[i]) <- alpha[region_index[i], month[i]] +
                        inprod(b_temp[region_index[i],], w_temp[i,]) +
                        inprod(b_hum[region_index[i],],  w_hum[i,]) +
                        inprod(b_rain[region_index[i],], w_rain[i,]) +
                        inprod(gamma[region_index[i],], t[i,]) +   
                        inprod(delta[region_index[i], year[i],], s[i,])
      
      # Deviance calculation
      log_like[i] <- y[i] * log(lambda[i]) - lambda[i] - logfact(y[i])
      dev[i] <- -2 * log_like[i]
    }
    
    # Total deviance for DIC calculations
    D <- sum(dev[1:n_data])
    pD <- 0  
    DIC <- 0 

    for(r in 1:n_reg){
      # Monthly intercept for each region
      for(m in 1:12){
        alpha[r, m] ~ dnorm(0, 0.001)
      }
      # Hierarchical climate effects with partial pooling
      for(j in 1:n_coef){
        b_temp[r, j] <- mu_temp[j] + sigma_temp[j] * theta_temp[r, j]
        theta_temp[r, j] ~ dnorm(0, tau_temp)

        b_hum[r, j]  <- mu_hum[j] + sigma_hum[j] * theta_hum[r, j]
        theta_hum[r, j] ~ dnorm(0, tau_hum)

        b_rain[r, j] <- mu_rain[j] + sigma_rain[j] * theta_rain[r, j]
        theta_rain[r, j] ~ dnorm(0, tau_rain)
      }
      # Trend effects
      for(tt in 1:n_trend){
        gamma[r, tt] ~ dnorm(0, 0.001)
      }
      # Seasonal effects by year
      for(yy in 1:n_year){
        for(ss in 1:n_seas){
          delta[r, yy, ss] ~ dnorm(0, 0.001)
        }
      }
    }

    # Hyperpriors for the climate effects parameters
    for(j in 1:n_coef){
      mu_temp[j]  ~ dnorm(0, 0.001)
      sigma_temp[j] ~ dunif(0,10)

      mu_hum[j]   ~ dnorm(0, 0.001)
      sigma_hum[j] ~ dunif(0,10)

      mu_rain[j]  ~ dnorm(0, 0.001)
      sigma_rain[j] ~ dunif(0,10)
    }

    tau_temp ~ dgamma(0.1, 0.1)
    tau_hum  ~ dgamma(0.1, 0.1)
    tau_rain ~ dgamma(0.1, 0.1)
  }
  "
  
  # Initial values for Model 4
  model4_inits <- function() {
    list(
      alpha      = matrix(0, nrow = model4_data$n_reg, ncol = 12),
      theta_temp = matrix(0, nrow = model4_data$n_reg, ncol = model4_data$n_coef),
      theta_hum  = matrix(0, nrow = model4_data$n_reg, ncol = model4_data$n_coef),
      theta_rain = matrix(0, nrow = model4_data$n_reg, ncol = model4_data$n_coef),
      mu_temp    = rep(0, model4_data$n_coef),
      mu_hum     = rep(0, model4_data$n_coef),
      mu_rain    = rep(0, model4_data$n_coef),
      sigma_temp = rep(1, model4_data$n_coef),
      sigma_hum  = rep(1, model4_data$n_coef),
      sigma_rain = rep(1, model4_data$n_coef),
      tau_temp   = 1,
      tau_hum    = 1,
      tau_rain   = 1
    )
  }
  
  model4_params <- c(
    "alpha", "b_temp", "b_hum", "b_rain",
    "gamma", "delta",
    "mu_temp", "sigma_temp", "tau_temp",
    "mu_hum", "sigma_hum", "tau_hum",
    "mu_rain", "sigma_rain", "tau_rain",
    "D", "log_like","dev"
  )
  
  # Fit Model 4 in JAGS
  tryCatch({
    jags_model4 <- jags.model(
      textConnection(model4_jags),
      data  = model4_data,
      inits = model4_inits,
      n.chains = n.chains
    )
    
    update(jags_model4, n.iter$model4 / 2)
    model4_samples <- coda.samples(
      jags_model4,
      variable.names = model4_params,
      n.iter = n.iter$model4 / 2
    )
    
    # Calculate DIC components for Model 4
    D_posterior <- do.call(rbind, lapply(model4_samples, function(x) x[, "D"]))
    D_mean <- mean(D_posterior)
    D_hat <- 2 * D_mean - mean(apply(D_posterior, 1, sum))
    p_D <- D_mean - D_hat
    DIC <- D_hat + 2 * p_D
    
    model4_dic <- list(
      D_mean = D_mean,
      D_hat = D_hat,
      p_D = p_D,
      DIC = DIC
    )
    
    save(model4_samples, model4_dic, file = "new_output_acc/model4_results.RData")
    cat("Model 4 completed and saved. DIC =", DIC, "\n")
    
  }, error = function(e) {
    cat("\nError in Model 4:", conditionMessage(e), "\n")
  })
}
