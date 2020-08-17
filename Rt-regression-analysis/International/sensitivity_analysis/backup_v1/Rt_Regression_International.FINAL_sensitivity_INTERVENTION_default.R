#----------------------------------------------------------------
# Project: Lin Lab - Covid19
# Analysis: R code for running the sensitivity analysis for INTERVENTION_ADDITIONAL parameters
# Author: Hui Li, Xihao Li, Derek Shyr, Zilin Li
# Requirements: cleaned data for regression
#----------------------------------------------------------------

rm(list=ls())

setwd("/n/holystore01/LABS/xlin/Lab/covid19/analysis/sensitivity/International")
# data management
library(data.table)
library(broom)
library(dplyr)
library(tidyr)
library(Matrix)
library(stringr)
library(usdm)

# graphing
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)


# model fitting
library(geepack)
library(splines)
library(lme4)

#### ===========================
#### CUSTOMIZABLE PARAMETERS
#### ===========================
OUT_TABLE_NAME = paste0("/n/holystore01/LABS/xlin/Lab/covid19/sensitivity/International/INTERVENTION_ADDITIONAL/International_Rt_Intervention_default_Output.csv")

# ** denotes parameters we may want to vary for sensitivity analysis. 

# "Outbreak" parameters.
OUTBREAK_PARAMS = list(
  MIN_CASES = 50 # ** Number of cases considered an "outbreak".
)

# Time trend (date) spline parameters 
SPLINE_PARAMS <- list(
  KNOTS_DATE = 8, # ** Knots for calendar date cubic spline.
  KNOTS_DS = 0 # ** Knots for days-since-outbreak cubic spline.
)

# Intervention lags
INTERVENTION_PARAMS <- list(
  DAYS_LAG = 14, # ** How many days before interventions affect Rt?
  VARS = c("workplace_closing", "stay_home_order", "mask_mandate"), #c("containmentIndex", "healthTestingIndex", "economicSupportIndex"), #"international_travel_control",
  VAR_TYPE = "scale" # factor: categorical variables; scale: continuous variables
)

# Serial interval (SI) parameters. 
SI_PARAMS = list(
  MEAN = 5.2, # ** SI mean.
  SD = 5.5, # ** SI standard deviation. 
  NTS = 30 # Maximum no. days infectious. 
)

# Confirmation-infection lag "deconvolution". 
LAG_PARAMS = list(
  USE_LAG = FALSE, # ** Use the random lag deconvolution?
  MAX_DAYS = 14, # Maximum lag period. 
  SHAPE = 7, # ** Lag distribution: Gamma rate. 
  RATE = 1 # ** Lag distribution: Gamma rate. 
)

# Standard error estimator
SE_ESTIMATOR = list(
  VAR_TYPE = "vbeta.naiv" # ** vbeta.naiv: Model-based SE. vbeta: Sandwich estimator.
)

# Date range
DATE_RANGE = list(
  MIN = "2020-03-05",
  MAX = "2020-07-20"
)

# Output parameters
OUTPUT_OPTION = list(
  CONVERT_COEF = TRUE,
  INTERVENTION_ONLY = TRUE
)

#### ===========================
#### LOAD CLEANED DATASET
#### ===========================
source("CalculateCumulativeIncidence.R")
source("deconv.R")

date_stamp <- "2020_08_09"
dt <- as.data.frame(fread(paste0("regression_intervention_", date_stamp, ".csv"), sep=",", header = TRUE, fill = TRUE))
setorder(dt, countryName, date)


#### ===========================
#### HELPER FUNCTIONS
#### ===========================
days_since <- function(start_date, dates, nz = TRUE){
  out <- as.numeric(as.Date(dates, format = "%Y-%m-%d") - as.Date(start_date, format = "%Y-%m-%d"))
  if(nz){out[out < 0] <- 0}
  return(out)
}

pval_mark <- function(pval){
  ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))
}

center <- function(x) scale(x, scale = FALSE, center = TRUE)

cap <- function(x, cp, f = `<=`) (x - cp)*f(x,cp) + cp

inv <- function(x, xn = names(x)){names(xn) <- x; xn}

# add new terms to a formula object
add_formula <- function(eqn, x, is.char = FALSE){
  if(!is.char){
	x <- gsub('"', "", deparse(substitute(x)))
  }
  as.formula(paste(as.character(eqn)[2] , '~', as.character(eqn)[3], '+', x))
}

tidyCoef <- function(x){
  require(data.table)
  out <- as.data.table(tidy(x))
  setnames(out, c("term", "estimate", "se", "stat", "pval"))
  out
}

# function to sum over sliding windows
slideSum <- function(x, n, fill = 0){
  x[is.na(x)] <- fill
  csx = cumsum(x)
  csx - data.table::shift(csx, n = n, type = 'lag', fill = fill)
}

# function to average over sliding windows
slideMean <- function(x, n, fill = 0, count_na = FALSE){
  if(count_na){
	denom <- n
  }else{
	# number of non-missing points in window
	denom <- slideSum(!is.na(x), n = n, fill = 0)
  }
  x[is.na(x)] <- fill
  csx = cumsum(x)
  out = (csx - data.table::shift(csx, n = n, type = 'lag', fill = fill))/denom
  out[denom == 0] = NA
  out
}


#### ===========================
#### PREPARE REGRESSION DATA
#### ===========================
# Subset to the specified date range
dt <- dt[dt$date <= DATE_RANGE$MAX & dt$date >= DATE_RANGE$MIN, ]

# Set any negative daily increase values to 0
dt$new_cases <- dt$positiveIncrease*(dt$positiveIncrease>=0)

# Random lag for case increase
if( LAG_PARAMS$USE_LAG ){
  dt$reported_new_cases <- as.numeric(dt$new_cases)
  dt$reported_positive <- as.numeric(dt$positive)
  dt$new_cases <- as.numeric(dt$new_cases + 1e-10)
  dt$positive <- as.numeric(dt$positive + 1e-10)
  dt <- as.data.table(dt)
  dt <- dt[order(countryName,date)]
  dt[,new_cases := deconvGammaLags(new_cases),by=countryName]
  dt[,positive := cumsum(new_cases),by=countryName]
}

dt <- subset(dt, !is.na(new_cases))
setorder(dt, countryName, date)

# Calculate the cumulative daily incidence
dt <- as.data.table(dt)
dt[, Lam := getLambda(
  x = new_cases,
  t = days_int, 
  mu = SI_PARAMS$MEAN, sd = SI_PARAMS$SD, NTS = SI_PARAMS$NTS
), by = countryName]


# Calculate days since "outbreak"
dt[, days_s_outbreak := days_since(min(date[positive >= OUTBREAK_PARAMS$MIN_CASES]), date)
, by = countryName]


# Update intervention and mobility metric covariates
dt <- dt[order(countryName,date)]

# Lagged intervention terms
intervention_names <- c("workplace_closing","gathering_restriction",
						 "stay_home_order","international_travel_control",
						 "mask_mandate", 
						 "containmentIndex","healthTestingIndex", "economicSupportIndex")

max_val_category <- c("workplace_closing3","gathering_restriction4",
                      "stay_home_order3","international_travel_control4",
                      "mask_mandate1")

dt[,
  (paste0(c(intervention_names), "_lag")) := lapply(.SD, data.table::shift, n = INTERVENTION_PARAMS$DAYS_LAG, fill = 0)
, by=countryName, .SDcols = c(intervention_names)]


# Additional QC: for each country
# drop the initial 7 days --> unreliable Lambda estimates
dt <- as.data.frame(dt)
dt_reg <- as.data.frame(dt %>%
						  group_by(countryName) %>%
						  arrange(countryName, date) %>%
						  mutate(days_within_country = sequence(n())))

dt_reg <- dt_reg[dt_reg$days_within_country > 7 & dt_reg$positive >= OUTBREAK_PARAMS$MIN_CASES, ]

# Remove NA values for RHS
dt_reg <- as.data.frame(dt_reg)
dt_reg <- dt_reg[complete.cases(dt_reg[c("new_cases",
										 "Lam",
										 "days_int", 
										 "days_s_outbreak",
										 "weekday",
										 "continent",
										 "population", 
										 "population_density", 
										 "median_age", 
										 "pct_diabetes", 
										 "cvd_death_rate", 
										 "Physicians_per_1000", 
										 "hospital_beds_per_100k", 
										 "Health_exp_pct_GDP_2016",
										 "gdp_per_capita_PPP", paste0(INTERVENTION_PARAMS$VARS, "_lag"))]), ]

dt_reg_mob <- dt_reg[complete.cases(dt_reg["mobility_s7"]), ]


#### ===========================
#### MAIN FUNCTIONS
#### ===========================

geeglm_fit_pois <- function(dt, model_formula, gee_id_var = "countryName", corstr = "ar1", estimate_type = "geese"){
  
  #' @param dt Data object for fitting models.
  #' @param model_formula Model formula for GEE GLM.
  #' @param gee_id_var Cluster variable name.
  #' @param corstr Correlation structure type.
  #' @param estimate_type What estimation method to use: geese is used by default. 
  
  # mold gee-id-var and sort by it
  dt <- as.data.frame(dt)
  dt$gee_id_var <- as.integer(factor(dt[[gee_id_var]], order = TRUE))
  dt <- as.data.table(dt)
  setorder(dt, gee_id_var)
  
  if (estimate_type == "glm") {
    fit <- geepack::geeglm(formula = model_formula,
                           data = dt, 
                           id = gee_id_var, 
                           family = poisson(),
                           corstr = corstr,
                           control = geese.control("epsilon" = 1e-4, 
                                                   "maxit" = 50000, 
                                                   trace = TRUE))
  } else if (estimate_type == "geese") {
    fit <- geepack::geese(formula = model_formula,
                          data = dt, 
                          id = gee_id_var, 
                          family = poisson(),
                          corstr = corstr,
                          control = geese.control("epsilon" = 1e-4, 
                                                  "maxit" = 50000, 
                                                  trace = TRUE))
  }
  
  return(fit)
}


#### ===========================
#### MODEL FITTING
#### ===========================

# Fixed Effects
fe_eqn = new_cases ~ 
  offset(log1p(Lam)) +
  factor(weekday) +
  factor(continent) + 
  center(population) +
  center(population_density) +
  center(median_age) +
  center(pct_diabetes) +
  center(cvd_death_rate)+ 
  center(Physicians_per_1000) +
  center(hospital_beds_per_100k) +
  center(Health_exp_pct_GDP_2016) +
  center(gdp_per_capita_PPP) +
  bs(days_s_outbreak, df = 3 + SPLINE_PARAMS$KNOTS_DS) +
  bs(days_int, df = 3 + SPLINE_PARAMS$KNOTS_DATE)

# Date range
as.Date("2020-01-01")+min(dt_reg$days_int) #03-15
as.Date("2020-01-01")+max(dt_reg$days_int) #07-19

# Anchor the order of categorical variable
dt_reg$continent <- as.factor(dt_reg$continent)
dt_reg$continent <- factor(dt_reg$continent,
                                 levels = c("Asia", "Africa","Europe",
                                            "Oceania", "North America",
                                            "Central America", "South America"))

dt_reg$weekday <- as.factor(dt_reg$weekday)
dt_reg$weekday <- factor(dt_reg$weekday, 
                               levels = c("Monday", "Tuesday", "Wednesday",
                                          "Thursday", "Friday", "Saturday", "Sunday"))

# FIT JOINT MODEL WITHOUT MOBILITY
eqn = add_formula(fe_eqn, paste(paste0(INTERVENTION_PARAMS$VAR_TYPE, "(", INTERVENTION_PARAMS$VARS, "_lag)"), collapse = ' + '), TRUE)

joint_fit <- geeglm_fit_pois(dt = dt_reg, 
                           model_formula = eqn,
                           gee_id_var = "countryName", 
                           corstr = "ar1", 
                           estimate_type = "glm")

# FIT MARGINAL MODEL WITHOUT MOBILITY
marginal_fits <- lapply(paste0(INTERVENTION_PARAMS$VARS, "_lag"), function(x){
	eqn = add_formula(fe_eqn, paste0(INTERVENTION_PARAMS$VAR_TYPE, "(", x, ")"), TRUE)
	geeglm_fit_pois(dt = dt_reg, 
                    model_formula = eqn,
                    gee_id_var = "countryName", 
                    corstr = "ar1", 
                    estimate_type = "glm")
})

names(marginal_fits) <- INTERVENTION_PARAMS$VARS

# FIT JOINT MODEL WITH MOBILITY
eqn = add_formula(fe_eqn, "mobility_s7", TRUE)
eqn = add_formula(eqn, paste(paste0(INTERVENTION_PARAMS$VAR_TYPE, "(", INTERVENTION_PARAMS$VARS, "_lag)"), collapse = ' + '), TRUE)

joint_wmob_fit <- geeglm_fit_pois(dt = dt_reg_mob, 
                           model_formula = eqn,
                           gee_id_var = "countryName", 
                           corstr = "ar1", 
                           estimate_type = "glm")

# FIT MARGINAL MODEL WITH MOBILITY
marginal_wmob_fits <- lapply(paste0(INTERVENTION_PARAMS$VARS, "_lag"), function(x){
	eqn = add_formula(fe_eqn, "mobility_s7", TRUE)
	eqn = add_formula(eqn, paste0(INTERVENTION_PARAMS$VAR_TYPE, "(", x, ")"), TRUE)
	geeglm_fit_pois(dt = dt_reg_mob, 
                    model_formula = eqn,
                    gee_id_var = "countryName", 
                    corstr = "ar1", 
                    estimate_type = "glm")
})

names(marginal_wmob_fits) <- INTERVENTION_PARAMS$VARS

fit_list <- list()
ii <- 1

fit_list[[ii]] <- tidyCoef(joint_fit)
fit_list[[ii]]$se <- sqrt(diag(joint_fit$geese[[SE_ESTIMATOR$VAR_TYPE]]))
fit_list[[ii]]$stat <- fit_list[[ii]]$estimate / fit_list[[ii]]$se
fit_list[[ii]]$pval <- (1-pnorm(q=abs(fit_list[[ii]]$stat), lower.tail = TRUE))*2
fit_list[[ii]]$pval_sym <- pval_mark(fit_list[[ii]]$pval)
fit_list[[ii]]$model <- "Interventions"
fit_list[[ii]]$type <- "Joint"
ii <- ii + 1

for(x in INTERVENTION_PARAMS$VARS){
  fit_list[[ii]] <- tidyCoef(marginal_fits[[x]])
  fit_list[[ii]]$se <- sqrt(diag(marginal_fits[[x]]$geese[[SE_ESTIMATOR$VAR_TYPE]]))
  fit_list[[ii]]$stat <- fit_list[[ii]]$estimate / fit_list[[ii]]$se
  fit_list[[ii]]$pval <- (1-pnorm(q=abs(fit_list[[ii]]$stat), lower.tail = TRUE))*2
  fit_list[[ii]]$pval_sym <- pval_mark(fit_list[[ii]]$pval)
  fit_list[[ii]]$model <- "Interventions"
  fit_list[[ii]]$type <- "Marginal"
  ii <- ii + 1
}

fit_list[[ii]] <- tidyCoef(joint_wmob_fit)
fit_list[[ii]]$se <- sqrt(diag(joint_wmob_fit$geese[[SE_ESTIMATOR$VAR_TYPE]]))
fit_list[[ii]]$stat <- fit_list[[ii]]$estimate / fit_list[[ii]]$se
fit_list[[ii]]$pval <- (1-pnorm(q=abs(fit_list[[ii]]$stat), lower.tail = TRUE))*2
fit_list[[ii]]$pval_sym <- pval_mark(fit_list[[ii]]$pval)
fit_list[[ii]]$model <- "Interventions + Mobility"
fit_list[[ii]]$type <- "Joint"
ii <- ii + 1

for(x in INTERVENTION_PARAMS$VARS){
  fit_list[[ii]] <- tidyCoef(marginal_wmob_fits[[x]])
  fit_list[[ii]]$se <- sqrt(diag(marginal_wmob_fits[[x]]$geese[[SE_ESTIMATOR$VAR_TYPE]]))
  fit_list[[ii]]$stat <- fit_list[[ii]]$estimate / fit_list[[ii]]$se
  fit_list[[ii]]$pval <- (1-pnorm(q=abs(fit_list[[ii]]$stat), lower.tail = TRUE))*2
  fit_list[[ii]]$pval_sym <- pval_mark(fit_list[[ii]]$pval)
  fit_list[[ii]]$model <- "Interventions + Mobility"
  fit_list[[ii]]$type <- "Marginal"
  ii <- ii + 1
}

fit_dt <- rbindlist(fit_list)

#### ===========================
#### CLEAN AND WRITE TO FILE
#### ===========================

if (OUTPUT_OPTION$INTERVENTION_ONLY) {
  remove_scale <- function(x){
    x <- gsub("scale\\(", "", x)
    x <- gsub("factor\\(", "", x)
    x <- gsub("_lag\\)", "", x)
    x <- gsub("[0-9]$", "", x)
  }
  
  fit_dt[, tag := remove_scale(term),]
  fit_out <- fit_dt[fit_dt$tag %in% c(intervention_names, "mobility_s") ]
  
  if (OUTPUT_OPTION$CONVERT_COEF & INTERVENTION_PARAMS$VAR_TYPE == "scale") {
    scale_factor <- sapply(INTERVENTION_PARAMS$VARS, 
                           function(x) (max(dt_reg[[paste0(x, "_lag")]]) - min(dt_reg[[paste0(x, "_lag")]]))/sd(dt_reg[[paste0(x, "_lag")]], na.rm=TRUE))
    fit_out$scaled_coeff <- exp(scale_factor[fit_out$tag]*fit_out$estimate) 
  }
  
} else {
  fit_out <- fit_dt
}

fwrite(fit_out, OUT_TABLE_NAME, sep = ',', col.names = TRUE, quote =TRUE)
