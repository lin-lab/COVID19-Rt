# COVID19-Rt

Code to calculate time-varying effective reproduction number (Rt), case rate, and death rate for COVID-19.

## Structure of Code

+ `pois_reg_estimates`: Implements our new Poisson regression-based method, as described in [our preprint](https://www.medrxiv.org/content/10.1101/2021.03.12.21253496v1).
+ `initial_estimates`: Initial Rt estimates calculated using [EpiEstim](https://cran.r-project.org/package=EpiEstim) ([Cori, A., et al., 2013](https://doi.org/10.1093/aje/kwt133))
+ `helpers`: Helper functions.
+ `preprocess_data`: Code for preprocessing the data.
+ `sens_analysis`: Sensitivity analysis

## Obtaining Rt Estimates

Within `pois_reg_estimates` and `initial_estimates` is a README with links to the latest Rt estimates.

## Misc

The `old_sensitivity` branch contains some sensitivity analysis code for Rt regression, which are depreciated.
