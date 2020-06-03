#!/usr/bin/R

#----------------------------------------------------------------
# Project: Lin Lab - Covid19
# Analysis: Country-level Rt regression
# Description: Implementation of GEE-GLM model fit
# Author: Hui Li
# Requirements: cleaned country level data
#----------------------------------------------------------------

rm(list=ls())

# data management
library(data.table)
library(dplyr)
library(tidyr)

# graphing
library(ggplot2)
library(grid)


# model fitting
library(geepack)
# library(lme4)
# library(splines)
# library(usdm)
# library(nlme)


#### ===========================
#### LOAD THE DATA
#### ===========================

# directories setup (if needed)

# load cleaned merged data
dt <- as.data.frame(fread("cleaned_country_level_2020_06_01.csv", 
                          sep=",", header = TRUE, fill = TRUE))


#### ======================
#### HELPER FUNCTIONS
#### ======================
# calculate days since an event with specified start_date
days_since <- function(start_date, dates, nz = TRUE){
  out <- as.numeric(as.Date(dates, format = "%Y-%m-%d") - as.Date(start_date, format = "%Y-%m-%d"))
  if(nz){out[out < 0] <- 0}
  return(out)
}

# Generate spline knots, fixed interval of days apart
seq_days <- function(x, delta){
  out <- seq(min(x) + delta, max(x), delta)
}

# Reverse lower case country names: for plotting aesthetics
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

# Number of ticks
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# Check how many countries + obs are left, 
# after removing certain variables from the regression
check_var_drop <- function(dt, var_remove){
  full_FE_cols <- c("mean_rt", 
                   "population", 
                   "population_density",
                   "pct_aged_70_older", 
                   # "pct_male_smokers",
                   # "pct_female_smokers",
                   # "pct_extreme_poverty",
                   "pct_diabetes", 
                   "cvd_death_rate", 
                   "Physicians_per_1000", 
                   "hospital_beds_per_100k", 
                   "Health_exp_pct_GDP_2016", 
                   "gdp_per_capita_PPP",
                   "stringencyIndex", 
                   "unique_continent",
                   colnames(dt)[colnames(dt) %like% "mobility_change"],
                   "days_s_lockdown", 
                   "days_int",
                   # "pct_handwashing_fac", 
                   # "total_tests_per_thousand",
                   "Life_expectancy_at_birth_both")
  
  sub_cols <- setdiff(full_FE_cols, var_remove)
  dt_temp <- dt[complete.cases(dt[sub_cols]), ]
  
  n_country <- length(unique(dt_temp$countryName))
  n_obs <- dim(dt_temp)[1]
  
  return(list(n_country, n_obs))
}


#### ===========================
#### PREPARE REGRESSION DATA
#### ===========================

# define continuous time variables
dt <- data.table(dt)
dt[,`:=`(
  # days since start of 2020
  days_int = days_since("2020-01-01", date),
  
  # days since lockdown order was placed
  days_s_lockdown = days_since(lockdown_date, date),

  # days since first reported case
  days_s_1p = days_since(min(date[positive > 0]), date)
), by = countryName]


# define dummy time variables
dt[,`:=`(
  lockdown = 1*(days_s_lockdown > 0),
  weekday = weekdays(as.Date(date))
),]

# knots values for B-spline terms
day_knots <- seq_days(dt$days_int, 14)

# Factorize variables for regression
dt <- as.data.frame(dt)
dt$lockdown <- as.factor(dt$lockdown)
dt$countryName <- as.factor(dt$countryName)
dt$weekday <- as.factor(dt$weekday)

# Add weekend indicator variable
dt$weekend <- ifelse(dt$weekday %in% c("Saturday", "Sunday"), 1, 0)

#### CREATE LAGGED INTERVENTION
# Option 1: push everything back by 5 days
dt$lockdown_unfac <- ifelse(dt$lockdown == 1, 1, 0)
dt$lag.lockdown_date <- as.Date(dt$lockdown_date) + 5
dt <- data.table(dt)
dt[,`:=`(
  lag.days_s_lockdown = days_since(lag.lockdown_date, date)
), by = countryName]

dt[,`:=`(
  lag.lockdown = 1*(lag.days_s_lockdown > 0)
),]

dt <- dt %>%
  group_by(countryName) %>%
  # mutate(lag.lockdown = dplyr::lag(lockdown_unfac, n = 5, default = NA)) %>%
  mutate(lag.stringencyIndex = dplyr::lag(stringencyIndex, n = 5, default = NA)) %>%
  mutate(lag.mobility = dplyr::lag(mobility_change_residential, n = 5, default = NA))

# Option 2: concatenate lockdown variable into 3 periods
dt <- data.table(dt)
dt[, stayhome_w1 := (days_s_lockdown >= 1 & days_s_lockdown <= 7),]
dt[, stayhome_w2 := (days_s_lockdown >= 8 & days_s_lockdown <= 14),]
dt[, stayhome_gt2 := (days_s_lockdown > 14),]


# Quality control of Rt
# Option 1: truncate Rt values (keep 1-99 percentile)
# dt <- (dt) %>%
#   filter(!is.na(mean_rt)) %>%
#   filter(mean_rt < quantile(mean_rt, 0.99) & mean_rt > quantile(mean_rt, 0.01))

# Option 2: Restrict country-dates to total cases > 50 <-- Avoid too much weighting on priors
dt <- dt[dt$positive >= 50, ]

# cross sectional data (use latest observation)
dt_country <- (dt) %>%
  group_by(countryName) %>%
  arrange(countryName, desc(date)) %>%
  summarise_all(first)
# 184 countries --> 156 countries


#### ===============
#### GEE functions
#### ===============

# fit the GEE model using approximate IVW for log(mean_rt)
geeglm_fit_ivw <- function(dt, model_formula, gee_id_var, corstr = "ar1",
                           weights = TRUE, sformula = "~ 1", estimate_type = "geese"){
  
  #' @param dt Data object for fitting models.
  #' @param model_formula Model formula for GEE GLM.
  #' @param gee_id_var Cluster variable name.
  #' @param corstr Correlation structure type.
  #' @param weights Whether to apply weighting or not. Weights: approximate for log(Y).
  #' @param sformula Formula for estimation of the dispersion.
  #' @param estimate_type What estimation method to use: geese is used by default. 
  
  # mold gee-id-var and sort by it
  dt <- as.data.frame(dt)
  dt$gee_id_var <- as.integer(factor(dt[[gee_id_var]], order = TRUE))
  dt <- as.data.table(dt)
  setorder(dt, gee_id_var)
  
  if (estimate_type == "glm") {
    if (weights) {
      
      # fit the model using geeglm
      fit = geeglm(data = dt,
                            formula = model_formula,
                            weights = (mean_rt/std_rt)^2,
                            id = gee_id_var,
                            corstr = corstr)
    } else {
      fit = geeglm(data = dt,
                            formula = model_formula,
                            id = gee_id_var,
                            corstr = corstr)
    }
  } else if (estimate_type == "geese") {
    
    # fit model using geese
    if (weights) {
      fit = geese(data = dt,
                           formula = model_formula,
                           weights = (mean_rt/std_rt)^2,
                           sformula = as.formula(sformula),
                           id = gee_id_var,
                           corstr = corstr)
    } else{
      fit = geese(data = dt,
                           formula = model_formula,
                           sformula = as.formula(sformula),
                           id = gee_id_var,
                           corstr = corstr)
    }
  }

  return(fit)
}

# fit the GEE model using IVW with estimated variance, iteratively
geeglm_fit_predVar_iter <- function(dt, model_formula, gee_id_var, 
                                    var_model_formula, corstr = "ar1", 
                                    iterMax = 1, tol = 1e-5){
  
  #' @param dt Data object for fitting models.
  #' @param model_formula Model formula for GEE GLM.
  #' @param gee_id_var Cluster variable name
  #' @param var_model_formula Model formula for predicting residual variances.
  #' @param corstr Correlation structure type.
  #' @param iterMax Maximum number of iterations for weighted least squares
  
  i=0
  
  # update formula for residual variances 
  var_model_formula <- as.formula(paste0("log_sq_resid ", as.character(var_model_formula)))
  
  last.coefs <- NA
  fit <- lm(model_formula, data = dt, na.action = na.exclude)
  coefs <- fit$coefficients
  
  while (is.na(last.coefs) || ( (max(coefs - last.coefs) > tol) &&  (i < iterMax))) {
    
    # create/update log squared residuals 
    dt$log_sq_resid <- log( (fit$residuals^2) )
    
    # predict the (log) residual variance 
    log_sq_resid_lm <- lm(var_model_formula, data = dt)
    
    # create/update inverse variance weights using the predicted log-residual-variance values.
    dt$ivw <- exp(-log_sq_resid_lm$fitted.values)
    
    # dt <- subset(dt, !is.na(predict(fit)))
    # print(dim(dt))
    
    # order by cluster variable
    dt <- as.data.frame(dt)
    dt$gee_id_var <- as.integer(factor(dt[[gee_id_var]], order = TRUE))
    dt <- as.data.table(dt)
    setorder(dt, gee_id_var)
    
    # model fitting
    fit <- geeglm( 
      formula = model_formula,
      id = gee_id_var, 
      corstr = corstr,
      weights = ivw, 
      data = dt, 
      family = gaussian(link = "identity"),
      control = geese.control("epsilon" = 1e-5, "maxit" = 50000) #, trace = TRUE)
    )
    
    # update iteration
    i <- i+1
    last.coefs <- coefs
    coefs <- fit$coefficients
    residuals <- fit$residuals
    
  }

  return(list(fit=fit, iterations=i))
}

# Wrapper function: Extract regression coefficients and prepare for bar plots
geeglm_coef_data <- function(fit, var_type, gee_obj, variable_name = NULL, estimate_name = "Estimate") {
  
  #' @param fit The model fit object, either geese or glm, specified by gee_obj.
  #' @param var_type Type of variance to use. 
  #' @param gee_obj What type of object is "fit". Slightly different cleaning steps for these two types.
  #' @param variable_name Labels of the variable.
  #' @param estimate_name The name of the coefficient estimates in the regression summary.
  
  if (gee_obj == "glm") {
    # extract coefficient from the fit object
    coef_table <- as.data.table(coef(summary(fit)), keep = "Variable")
    
    # filter out factor variable and time variables
    ind <- which(coef_table$Variable %like% "scale" | 
                   coef_table$Variable %like% "lockdown" | 
                   coef_table$Variable %like% "continent" | 
                   coef_table$Variable %like% "weekend")
    
    y_se <- sqrt(diag(fit$geese[[var_type]]))[ind]
    
  } else if (gee_obj == "geese") {
    # extract coefficient from the fit object
    coef_table <- as.data.table(fit$beta, keep = "Variable")
    colnames(coef_table) <- c("Variable", estimate_name)
    
    # filter out factor variable and time variables
    ind <- which(coef_table$Variable %like% "scale" | 
                   coef_table$Variable %like% "lockdown" | 
                   coef_table$Variable %like% "continent" | 
                   coef_table$Variable %like% "weekend")
    
    y_se <- sqrt(diag(fit[[var_type]]))[ind]
    
  }
 
  coef_table <- coef_table %>% filter(Variable %like% "scale" |
                                        Variable %like% "lockdown" | 
                                        Variable %like% "continent" | 
                                        Variable %like% "weekend")
  
  # clean up coefficient names and extract corresponding s.e. values 
  y <- coef_table[, estimate_name]
  
  if(is.null(variable_name)){
    coef_name <- gsub(pattern = "scale", replacement="",
                      x = coef_table$Variable)
    coef_name <- gsub(pattern="\\)", replacement="",
                      x = gsub(pattern="\\(", replacement="",
                               x = coef_name))
  }else{
    coef_name <- variable_name
  }
  
  # collect coefficients and s.e.
  plot_dt <- data.frame(name=coef_name, 
                        value=y, 
                        se=y_se)
  return(plot_dt)
}


#### ===============
#### MODEL FITTING
#### ===============

# Fixed Effects
fe_eqn = log(mean_rt) ~
  factor(weekend) +
  factor(unique_continent) + 
  scale(population) +
  scale(population_density) +
  scale(median_age) +
  scale(pct_diabetes) +
  scale(cvd_death_rate)+
  scale(Physicians_per_1000) +
  scale(hospital_beds_per_100k) +
  scale(Health_exp_pct_GDP_2016) +
  scale(gdp_per_capita_PPP) + 
  scale(Life_expectancy_at_birth_both) +
  scale(days_s_lockdown) +
  factor(lockdown) +
  scale(stringencyIndex) +
  scale(mobility_change_residential)


# Remove NA values for RHS (ONLY KEEP TERMS TO BE INCLUDED)
dt <- as.data.frame(dt)
dt_reg <- dt[complete.cases(dt[c("mean_rt", 
                                 "weekend", 
                                 "unique_continent",
                                 "population", 
                                 "population_density", 
                                 "median_age",
                                 "pct_diabetes", 
                                 "cvd_death_rate", 
                                 "Physicians_per_1000", 
                                 "hospital_beds_per_100k", 
                                 "Health_exp_pct_GDP_2016",
                                 "gdp_per_capita_PPP", 
                                 "Life_expectancy_at_birth_both",
                                 "days_s_lockdown",
                                 "lockdown",
                                 "stringencyIndex",
                                 "mobility_change_residential")]), ]
# check the date range
as.Date("2020-01-01")+min(dt_reg$days_int) #03-08
as.Date("2020-01-01")+max(dt_reg$days_int) #05-09

# set up time variable in model matrix
time_bs <- bs(dt_reg$days_int, knots = day_knots)
dt_reg[, paste0("D", seq(1,7))] <- time_bs

# clean up final equation for fitting
eqn = as.formula(paste0(fe_eqn[2], " ~ ", fe_eqn[3],
                        "+ D1 + D2 + D3 + D4 + D5 + D6 + D7"))

# 1. no weighting
fit_geelm_ar1 <- geeglm_fit_ivw(dt = dt_reg, 
                                    model_formula = eqn, 
                                    gee_id_var = "countryName",
                                    weights = FALSE,
                                    estimate_type = "geese")

# 2. use approximate weights first
fit_geelm_ar1_ivw <- geeglm_fit_ivw(dt = dt_reg, 
                            model_formula = eqn, 
                            gee_id_var = "countryName",
                            weights = TRUE,
                            estimate_type = "geese")

# 3. use approximate weights to estimate dispersion
fit_geelm_ar1_disp <- geeglm_fit_ivw(dt = dt_reg, 
                     model_formula = eqn, 
                     gee_id_var = "countryName",
                     sformula = "~ (std_rt^2)/(mean_rt^2)")

# 4. use residual variance as weighting
fit_geelm_ar1_predVar <- geeglm_fit_predVar_iter(dt = dt_reg,
                                      model_formula = eqn, 
                                      gee_id_var = "countryName",
                                      var_model_formula = "~ log((std_rt^2)/(mean_rt^2))")$fit

# 5. use iterative weighted least squares (default uses 10 iterations)
fit_geelm_ar1_predVar_iter <- geeglm_fit_predVar_iter(dt = dt_reg,
                                   model_formula = eqn,
                                   gee_id_var = "countryName",
                                   var_model_formula = "~ log((std_rt^2)/(mean_rt^2))",
                                   iterMax = 50, tol = 1e-5)$fit

#### ======================
#### REGRESSION PLOTS
#### ======================
# manually enter covariate names
variable_name <- c("Weekend or not", 
                   "Ref-Africa: Asia", 
                   "Ref-Africa: Europe", 
                   "Ref-Africa: North America",
                   "Ref-Africa: Oceania", 
                   "Ref-Africa: South America",
                   "Population", 
                   "Population density",
                   "Median age", 
                   "Diabetes prevalence", 
                   "CVD Death rate",
                   "# of Physicians per 1000", 
                   "# of hospital beds per 100k",
                   "Health expenditure as % of GDP",
                   "GDP per capita (PPP)",
                   "Life expectancy at birth",
                   "Days since lockdown",
                   "Lockdown or not",
                   "Government policy stringency",
                   "Mobility change in residential area")

# run for all three S.E. options
variance_types = c("vbeta", "vbeta.naiv")
                   # , "vbeta.ajs")
tags <- c("san", "naive")
          # , "jk")
labels <- c("Sandwich", "Model-based")
            # , "Jackknife")

# line up for graphing
fit_results <- list(fit_geelm_ar1, fit_geelm_ar1_ivw, fit_geelm_ar1_disp,
                    fit_geelm_ar1_predVar, fit_geelm_ar1_predVar_iter)
graph_names <- c("gee_no_weights_", "gee_ivw_", "gee_disp_", "gee_predVar_1_", "gee_predVar_iter_")
gee_types <- c("geese", "geese", "geese", "glm", "glm")

model_names <- c("No weights. ", "IVW (measurement error). ", "Estimate dispersion. ", 
                 "Single weight updates WLS. ", "Iterated WLS until convergence. ")

for (l in 1:5) {
  fit <- fit_results[[l]]
  fit_name <- graph_names[l]
  
  
  for (i in 1:2) {
    # 1. Plot results from geeglm
    plot_dt <- geeglm_coef_data(fit, var_type = variance_types[i], variable_name = variable_name, gee_obj = gee_types[l])
    plot <- ggplot(plot_dt) +
      geom_bar(aes(x=name, y=value), stat="identity", fill="royalblue", alpha=0.7) +
      geom_errorbar(aes(x=name, ymin=value-(1.96*se), ymax=value+(1.96*se), width=0.2),
                    colour="firebrick2", alpha=0.9, size=0.8) +
      # ylim(-5, 15) +
      coord_flip() +
      labs(x = 'Variables', y = expression("Coefficient" %+-% "1.96 SE"), 
           title = 'GEE-GLM Regression coefficients',
           subtitle = paste0("76 countries, 4086 observations. March 8 - May 9. \nCountry clusters, AR1 correlation structure. \n", 
                             model_names[l], labels[i], " SE.")) +
      theme(plot.title=element_text(hjust=0.5, size=30),
            plot.subtitle=element_text(size=15),
            axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0), hjust=0.5, size=15),
            axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0), vjust=0.5, size=15),
            axis.text.x=element_text(size=15),
            axis.text.y=element_text(size=15),
            legend.title=element_text(size=15),
            legend.text=element_text(size=15))
    ggsave(paste0(fit_name, tags[i], "_se_reg_coef_bar.png"), plot, width=15, height=10)
  }
}
