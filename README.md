# Spatially-informed Bayesian Distributed Lag Non-Linear Models (SB-DLNM) for Influenza in NSW

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue.svg)](https://www.r-project.org/)


## 📋 Overview

This repository contains the complete implementation of **Spatially-informed Bayesian Distributed Lag Non-Linear Models (SB-DLNM)** for investigating the complex relationships between meteorological factors and influenza transmission across New South Wales (NSW), Australia. 

Our framework quantifies the delayed and non-linear effects of temperature, relative humidity, and rainfall on influenza notifications, analyzing over 1.2 million laboratory-confirmed cases across 15 Local Health Districts (LHDs) from 2000-2023. The model captures both non-linear and delayed effects of weather exposures while accounting for spatial dependencies between neighboring health districts, providing crucial insights for public health planning and disease surveillance.

### 🎯 Key Features

- **Bayesian DLNM Framework**: Captures non-linear and lagged effects up to 4 weeks using cross-basis splines
- **Spatial Integration**: Conditional autoregressive (CAR) priors enable neighboring LHDs to borrow strength
- **Comprehensive Analysis**: Over 1.2 million laboratory-confirmed cases (2000-2023) linked to district-level meteorological data
- **Multiple Model Comparisons**: Four designs compared (case-crossover vs. time-series, with/without spatial pooling)
- **Early Warning Capability**: Framework supports district-level surge prediction weeks in advance
- **Publication-Ready Outputs**: Automated generation of figures, tables, and diagnostic plots

## 👥 Authors

**Mohammad Afzal Khan**¹*, Oyelola Adegboye²*, Shiyang Lyu¹, Kiki Maulana Adhinugraha³, Theophilus I. Emeto⁴⁵, and David Taniar¹

¹ Faculty of Information Technology, Monash University, Melbourne, VIC 3800, Australia  
² Menzies School of Health Research, Charles Darwin University, Darwin, NT 0800, Australia  
³ School of Computing and Information Technology, La Trobe University, Melbourne, VIC 3086, Australia  
⁴ Australian Institute of Tropical Health and Medicine, James Cook University, Townsville, QLD 4811, Australia  
⁵ College of Public Health, Medical and Veterinary Sciences, James Cook University, Townsville, QLD 4811, Australia

*Corresponding authors

## 📊 Data

### Data Sources
- **Influenza Data**: NSW Health Notifiable Conditions Information Management System (NCIMS)
- **Meteorological Data**: Australian Bureau of Meteorology (BoM)
- **Spatial Boundaries**: NSW Ministry of Health

### Data Description
- **Study Period**: January 2000 - December 2023 (24 years)
- **Total Cases**: Over 1.2 million laboratory-confirmed influenza notifications
- **Geographic Coverage**: 15 Local Health Districts (LHDs) in NSW
- **Temporal Resolution**: Monthly aggregation for analysis
- **Variables**: 
  - Laboratory-confirmed influenza cases (Types A, B, and total)
  - Daily mean temperature (°C)
  - Daily mean relative humidity (%)
  - Daily total rainfall (mm)

## 🚀 Quick Start

### Prerequisites
- R ≥ 4.0.0
- Required R packages (see `requirements.R`)

### Key Dependencies
- **Statistical Modeling**: `dlnm`, `splines`, `survival`, `coda`
- **Spatial Analysis**: `sf`, `spdep`, `ape`
- **Data Processing**: `dplyr`, `tidyr`, `lubridate`, `data.table`
- **Visualization**: `ggplot2`, `plotly`, `viridis`, `patchwork`

### Installation

```bash
# Clone the repository
git clone https://github.com/afzal0/SB-DLNM-Influenza_NSW.git
cd SB-DLNM-Influenza_NSW

# Install required packages
Rscript requirements.R
```

### Running the Analysis

For a complete analysis pipeline:

```r
# Run the entire analysis
source("Data Prepration.R")
source("Data Modelling.R")
source("Accuracy_Assessment.R")
source("Visualisation.R")
source("Opt-Vis1-DLNM_Curves.R")
source("Opt-Vis2-Cumulative Effects.R")
```

For detailed workflow instructions, see [WORKFLOW.md](WORKFLOW.md).

## 📁 Repository Structure

```
SB-DLNM-Influenza_NSW/
├── data/                      # Input data files
├── spatial_data/             # Spatial boundary files
├── output/                   # Processed data outputs
├── new_output/               # Model results and diagnostics
├── new_output_acc/           # Model accuracy assessments
├── Visualisations/           # Generated figures
├── Data Prepration.R         # Data preprocessing script
├── Data Modelling.R          # Main modeling script
├── Accuracy_Assessment.R     # Model evaluation
├── Visualisation.R           # Enhanced visualizations
├── Opt-Vis1-DLNM_Curves.R   # Response curve analysis
├── Opt-Vis2-Cumulative Effects.R  # Cumulative effects
├── requirements.R            # Package installation
├── WORKFLOW.md              # Detailed workflow guide
└── README.md                # This file
```

## 📈 Results

### Key Findings
- **Temperature as dominant driver**: Cumulative relative risk (RR) peaked at 1.9 (95% CI: 1.4-2.5) near 21°C
- **Protective effects at extremes**: RR < 1.0 for temperatures <10°C or >28°C
- **Humidity effects**: Modest, location-specific effects (RR 1.3-1.6 between 55-75%)
- **Rainfall**: Only sporadically associated with influenza risk
- **Spatial patterns**: Coherent coastal hot-spots identified through spatial pooling
- **Lag structure**: Effects manifest over 0-4 weeks post-exposure

### Model Performance
| Model | Description | DIC |
|-------|-------------|-----|
| Model 1 | Case-crossover without spatial pooling | - |
| Model 2 | Time-series without spatial pooling | - |
| Model 3 | Case-crossover with spatial pooling | **153 (Best)** |
| Model 4 | Time-series with spatial pooling | - |

### Practical Implications
- Spatial pooling removed implausible extremes in data-sparse western districts
- Framework enables early-warning dashboards for district-level surge prediction
- Temperature monitoring crucial for influenza preparedness in NSW

## 🛠️ Methodology

This implementation builds upon the Spatial Bayesian Distributed Lag Non-Linear Models (SB-DLNM) framework developed by Quijal-Zamorano et al. (2024), adapting it for influenza surveillance in NSW.

### Statistical Framework
- **Distributed Lag Non-Linear Models (DLNM)**: Cross-basis splines capturing non-linear and delayed effects (0-4 weeks)
- **Bayesian Framework**: Provides uncertainty quantification through credible intervals
- **Spatial Component**: Conditional autoregressive (CAR) priors enable neighboring LHDs to borrow strength
- **MCMC Sampling**: 10,000 iterations with 5,000 burn-in
- **Model Selection**: Deviance Information Criterion (DIC) for model comparison
- **Temporal Resolution**: Monthly aggregation of daily meteorological and case data

### Model Specifications
1. **Model 1**: Independent B-DLNM with case-crossover design
2. **Model 2**: Independent B-DLNM with time-series design  
3. **Model 3**: Spatially pooled DLNM with case-crossover design
4. **Model 4**: Spatially pooled DLNM with time-series design

### Key Methodological Contributions
- Analysis of over 1.2 million laboratory-confirmed influenza cases across 24 years
- Integration of multiple meteorological factors with identification of temperature as dominant driver
- Spatial pooling to address data sparsity in western districts and identify coastal hot-spots
- Development of early-warning framework for district-level surge prediction
- Quantification of optimal temperature range (21°C) for influenza transmission in temperate Australia

## 🌟 Applications

This framework can be adapted for:
- **Public Health Surveillance**: Real-time monitoring and early warning systems
- **Resource Planning**: Hospital capacity and vaccine distribution optimization
- **Climate-Health Research**: Understanding weather-disease relationships in other regions
- **Policy Development**: Evidence-based interventions for influenza control
- **Forecasting Models**: Integration with machine learning for enhanced prediction

## 📝 Citation

If you use this code or data in your research, please cite both our work and the original SB-DLNM methodology:

### Our Implementation
```bibtex
@software{khan2024sbdlnm,
  title = {Spatially-informed Bayesian Distributed Lag Non-Linear Models for Influenza in NSW},
  author = {Khan, Mohammad Afzal and Adegboye, Oyelola and Lyu, Shiyang and 
            Adhinugraha, Kiki Maulana and Emeto, Theophilus I. and Taniar, David},
  year = {2024},
  url = {https://github.com/afzal0/SB-DLNM-Influenza_NSW},
  version = {1.0.0}
}
```

### Original SB-DLNM Methodology
```bibtex
@article{quijal2024sbdlnm,
  title = {Spatial Bayesian distributed lag non-linear models (SB-DLNM) for small-area exposure-lag-response epidemiological modelling},
  author = {Quijal-Zamorano, Marcos and Martinez-Beneito, Miguel A and Ballester, Joan and Marí-Dell'Olmo, Marc},
  journal = {International Journal of Epidemiology},
  volume = {53},
  number = {3},
  pages = {dyae061},
  year = {2024},
  doi = {10.1093/ije/dyae061}
}
```

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🤝 Contributing

We welcome contributions! Please feel free to submit a Pull Request.

## 📞 Contact

For questions or collaborations:
- **Mohammad Afzal Khan**: mkha0168@student.monash.edu
- **Repository Issues**: [GitHub Issues](https://github.com/afzal0/SB-DLNM-Influenza_NSW/issues)

## 🙏 Acknowledgments

- NSW Health for providing influenza surveillance data
- Australian Bureau of Meteorology for meteorological data
- Monash University for computational resources

### Methodological Foundation

This work builds upon the Spatial Bayesian Distributed Lag Non-Linear Models (SB-DLNM) framework:

Quijal-Zamorano, M., Martinez-Beneito, M. A., Ballester, J., & Marí-Dell'Olmo, M. (2024). Spatial Bayesian distributed lag non-linear models (SB-DLNM) for small-area exposure-lag-response epidemiological modelling. *International Journal of Epidemiology*, 53(3), dyae061. https://doi.org/10.1093/ije/dyae061

We acknowledge and thank the authors for making their methodology available, which enabled this adaptation for influenza surveillance in NSW.
