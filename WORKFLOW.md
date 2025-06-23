# SB-DLNM Analysis Workflow Guide

This document provides a step-by-step guide for running the Spatially-informed Bayesian Distributed Lag Non-Linear Models (SB-DLNM) analysis for influenza in NSW.

## Prerequisites

1. **R Installation**: Ensure R version 4.0.0 or higher is installed
2. **Package Installation**: Run `source("requirements.R")` to install all required packages
3. **Data Files**: Ensure all data files are present in the `data/` and `spatial_data/` directories

## Analysis Workflow

### Step 1: Data Preparation (Runtime: ~5-10 minutes)

```r
source("Data Prepration.R")
```

**What this does:**
- Loads influenza case data and meteorological variables
- Processes spatial boundaries for NSW Local Health Districts
- Creates spatial neighborhood matrices
- Generates crossbasis matrices for DLNM analysis
- Performs time series decomposition
- Saves preprocessed data to `output/` directory

**Key outputs:**
- `output/list_neighbours.RData` - Spatial neighborhood relationships
- Preprocessed data objects in R workspace

### Step 2: Statistical Modeling (Runtime: ~2-4 hours)

```r
source("Data Modelling.R")
```

**What this does:**
- Fits four different Bayesian DLNM models:
  - Model 1: Independent B-DLNM with case-crossover design
  - Model 2: Independent B-DLNM with time-series design
  - Model 3: Spatially pooled DLNM with case-crossover design
  - Model 4: Spatially pooled DLNM with time-series design
- Uses MCMC sampling (10,000 iterations) for parameter estimation
- Saves model outputs for each LHD

**Key outputs:**
- `new_output/Model_*_LHD_*.RData` - Model results for each LHD
- MCMC samples and parameter estimates

**Note**: This is the most computationally intensive step. Consider running overnight or on a high-performance system.

### Step 3: Model Evaluation (Runtime: ~10-15 minutes)

```r
source("Accuracy_Assessment.R")
```

**What this does:**
- Calculates model selection criteria:
  - Deviance Information Criterion (DIC)
  - Bayesian Information Criterion (BIC)
  - Quasi-AIC (QAIC)
- Compares performance across all four models
- Generates model comparison plots

**Key outputs:**
- `new_output_acc/model_comparison_plot.pdf` - Visual comparison of models
- Console output with model selection statistics

### Step 4: Enhanced Visualization (Runtime: ~20-30 minutes)

```r
source("Visualisation.R")
```

**What this does:**
- Creates risk maps by Local Health District
- Performs uncertainty quantification
- Conducts spatial clustering analysis (Moran's I)
- Generates MCMC diagnostic plots
- Produces publication-ready figures

**Key outputs:**
- `Visualisations/risk_maps_*.pdf` - Spatial risk visualizations
- `Visualisations/uncertainty_plots_*.pdf` - Uncertainty analysis
- `new_output/lhd_plots/` - LHD-specific diagnostic plots

### Step 5: Response Curve Analysis (Runtime: ~10 minutes)

```r
source("Opt-Vis1-DLNM_Curves.R")
```

**What this does:**
- Generates lag-response curves for each meteorological variable
- Shows non-linear relationships between weather and influenza risk
- Creates 3D surface plots for visualization

**Key outputs:**
- `Visualisations/DLNM_response_curves_*.pdf` - Response curve visualizations

### Step 6: Cumulative Effects Analysis (Runtime: ~10 minutes)

```r
source("Opt-Vis2-Cumulative Effects.R")
```

**What this does:**
- Calculates cumulative effects over different lag periods
- Identifies critical temperature/humidity thresholds
- Produces summary effect plots

**Key outputs:**
- `Visualisations/cumulative_effects_*.pdf` - Cumulative effect plots

## Running the Complete Analysis

To run the entire analysis pipeline:

```r
# Set working directory to project root
setwd("/path/to/SB-DLNM-Influenza_NSW")

# Install packages
source("requirements.R")

# Run analysis pipeline
source("Data Prepration.R")
source("Data Modelling.R")
source("Accuracy_Assessment.R")
source("Visualisation.R")
source("Opt-Vis1-DLNM_Curves.R")
source("Opt-Vis2-Cumulative Effects.R")
```

## Troubleshooting

### Memory Issues
- If encountering memory errors during modeling, consider:
  - Increasing R memory limit: `memory.limit(size = 16000)` (Windows)
  - Running models for individual LHDs separately
  - Reducing MCMC iterations (though this may affect convergence)

### Missing Data Warnings
- Some LHDs may have missing data for certain time periods
- The models handle missing data through the Bayesian framework
- Check console output for specific warnings

### Convergence Issues
- Check MCMC diagnostic plots in `new_output/lhd_plots/`
- Look for Gelman-Rubin statistics < 1.1
- Consider increasing burn-in period if needed

## Customization Options

### Modifying Lag Periods
In `Data Prepration.R`, adjust:
```r
lag_max <- 30  # Maximum lag days (default: 30)
```

### Changing MCMC Settings
In `Data Modelling.R`, modify:
```r
n.iter <- 10000    # Total iterations
n.burnin <- 5000   # Burn-in period
n.thin <- 5        # Thinning interval
```

### Selecting Specific LHDs
To analyze specific LHDs only, modify the LHD list in `Data Prepration.R`.

## Expected Total Runtime

- **Minimal run** (single LHD, reduced iterations): ~1 hour
- **Standard run** (all LHDs, default settings): ~3-5 hours
- **Extended run** (increased iterations for publication): ~6-8 hours

## Contact

For questions about the implementation:
- Mohammad Afzal Khan (mkha0168@student.monash.edu)
- Code repository: https://github.com/afzal0/SB-DLNM-Influenza_NSW