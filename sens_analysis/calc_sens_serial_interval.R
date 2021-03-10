#Serial interval sensitivity: Based on calc_pois_estimates.R

################################################################################
## Setup
################################################################################
library(splines)
library(magrittr)
library(dplyr)
library(data.table)
library(logger)

options(warn = 2)

source("../helpers/helper_functions.R")
source("pois_functions.R")

################################################################################
## Read files
################################################################################
log_info("Start reading files...")
start_date <- as.Date("2020-03-01")
jhu_counties <- read_jhu("County", start_date = start_date)[!(startsWith(county, "Out of") | startsWith(county, "Unassigned")), ]
jhu_counties[, Combined_Key := paste0(county, ", ", stateName)]
jhu_states <- read_jhu("State", start_date = start_date)
jhu_states[, Combined_Key := stateName]
jhu_global <- read_jhu("Global", start_date = start_date)
jhu_subnational_orig <- read_jhu("Subnational", start_date = start_date)[Province_State != "Unknown", ]

################################################################################
## Data munging for US census divisions
################################################################################
census_div <- fread("us_census_divisions.csv")
census_div_uid <- data.table(Division = unique(census_div$Division),
                             UID = 840 + seq_len(uniqueN(census_div$Division)))
census_div_merged <- merge(census_div, census_div_uid, by = "Division")
us_div_merged <- merge(jhu_states, census_div_merged, by.x = "stateName",
                       by.y = "State")
us_div <- us_div_merged[, .(positive = sum(positive), death = sum(death),
                            population = sum(population),
                            UID = as.integer(unique(UID.y))),
                        by = .(Division, date)]
us_div[, `:=` (positiveIncrease = positive -
                      shift(positive, n = 1, type = "lag", fill = 0),
               deathIncrease = death -
                      shift(death, n = 1, type = "lag", fill = 0)),
       by = Division]
us_div[, Combined_Key := Division]

knots_touse <- c(seq(as.Date("2020-03-30"), max(jhu_states$date), by = 30),
                 as.Date("2020-12-07")) %>%
  sort()

jhu_global_nonus <- jhu_global[UID != 840, ]
jhu_us <- jhu_global[UID == 840, ]

################################################################################
## Data munging for subnational: discard first date
################################################################################

jhu_subnational <- jhu_subnational_orig %>%
  group_by(UID) %>%
  filter(date > min(date)) %>%
  data.table()

################################################################################
## Run the fitting
################################################################################

days_lag <- 7L
run_sens_analysis <- function(si_mean_in, si_sd_in) {
  
fit_lst <- list()
system.time({
  fit_lst$state <-
    fit_poisson_by_unit(jhu_states, min_incr = 20,
                        days_per_knot = 30, min_positive = 50,
                        SI_mean = si_mean_in, SI_sd = si_sd_in,
                        family = "negbin", knots = knots_touse) %>%
    lag_dates(lag = days_lag)

  fit_lst$us_div <-
    fit_poisson_by_unit(us_div, min_incr = 20,
                        days_per_knot = 30, min_positive = 50,
                        SI_mean = si_mean_in, SI_sd = si_sd_in,
                        family = "negbin", knots = knots_touse) %>%
    lag_dates(lag = days_lag)

  global_fit_nonus <- fit_poisson_by_unit(jhu_global_nonus, min_incr = 20,
                                          days_per_knot = 30, min_positive = 50,
                                          SI_mean = si_mean_in, SI_sd = si_sd_in,
                                          family = "negbin") %>%
    lag_dates(lag = days_lag)
  
  us_fit <- fit_poisson_by_unit(jhu_us, min_incr = 20,
                                days_per_knot = 30, min_positive = 50,
                                SI_mean = si_mean_in, SI_sd = si_sd_in,
                                family = "negbin", knots = knots_touse) %>%
    lag_dates(lag = days_lag)
  fit_lst$global <- rbind(global_fit_nonus, us_fit)

  fit_lst$subnational <-
    fit_poisson_by_unit(jhu_subnational, min_incr = 20,
                        days_per_knot = 30, min_positive = 50,
                        SI_mean = si_mean_in, SI_sd = si_sd_in,
                        family = "negbin") %>%
    lag_dates(lag = days_lag)

  fit_lst$county <-
    fit_poisson_by_unit(jhu_counties, min_incr = 20,
                        days_per_knot = 30, min_positive = 50,
                        SI_mean = si_mean_in, SI_sd = si_sd_in,
                        family = "negbin") %>%
    lag_dates(lag = days_lag)
})

################################################################################
## Fit exceptions. Sorry for the spaghetti code
################################################################################

update_fit <- function(fit_lst, updated, loc_key,
                       resolution = c("county", "global", "subnational", "state"),
                       metric = c("rt", "case", "death")) {
  resolution <- match.arg(resolution)
  metric <- match.arg(metric)

  tmp <- fit_lst[[resolution]][Combined_Key == loc_key, ]
  stopifnot(nrow(tmp) > 0)
  if (metric == "rt") {
    fit_lst[[resolution]][Combined_Key == loc_key,
                          `:=` (rt = NA, rt_lower = NA, rt_upper = NA)]
    fit_lst[[resolution]][Combined_Key == loc_key & date %in% updated$date,
                          `:=` (rt = updated$rt,
                                rt_lower = updated$rt_lower,
                                rt_upper = updated$rt_upper)]

  } else if (metric == "case") {
    fit_lst[[resolution]][Combined_Key == loc_key,
                          `:=` (case_rate = NA, case_lower = NA, case_upper = NA)]
    fit_lst[[resolution]][Combined_Key == loc_key & date %in% updated$date,
                          `:=` (case_rate = updated$case_rate,
                                case_lower = updated$case_lower,
                                case_upper = updated$case_upper)]
  } else if (metric == "death") {
    fit_lst[[resolution]][Combined_Key == loc_key,
                          `:=` (death_rate = NA, death_lower = NA, death_upper = NA)]
    fit_lst[[resolution]][Combined_Key == loc_key & date %in% updated$date,
                          `:=` (death_rate = updated$death_rate,
                                death_lower = updated$death_lower,
                                death_upper = updated$death_upper)]
  } else {
    stop("Should not get here.")
  }
}

## Rt Exceptions

out <- with(jhu_counties[Combined_Key == "Okfuskee, Oklahoma"],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        min_positive = 50, min_incr = 1, do_remove = FALSE,
                        prop_zeros = 0,
                        knots = knots_touse, family = "negbin")
)[, .(date, rt = outcome_hat, rt_lower = ci_lower, rt_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Okfuskee, Oklahoma",
           resolution = "county", metric = "rt")

out <- with(jhu_counties[Combined_Key == "Rusk, Texas"],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        min_positive = 50, min_incr = 20,
                        knots = knots_touse, do_remove = TRUE,
                        prop_zeros = 0,
                        family = "negbin")
)[, .(date, rt = outcome_hat, rt_lower = ci_lower, rt_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Rusk, Texas", resolution = "county",
           metric = "rt")

out <- with(jhu_counties[Combined_Key == "Gillespie, Texas"],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        min_positive = 50, min_incr = 20,
                        knots = knots_touse,
                        do_remove = FALSE, prop_zeros = 0,
                        family = "negbin")
)[, .(date, rt = outcome_hat, rt_lower = ci_lower, rt_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Gillespie, Texas", resolution = "county",
           metric = "rt")

# Case Rate exceptions

out <- with(jhu_subnational[Combined_Key == "Baleares, Spain", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population, prop_zeros = 0,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 30,
                        family = "negbin")
)[, .(date, case_rate = outcome_hat, case_lower = ci_lower, case_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Baleares, Spain",
           resolution = "subnational", metric = "case")

out <- with(jhu_subnational[Combined_Key == "C. Valenciana, Spain", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population, prop_zeros = 0.2,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 30,
                        family = "negbin")
)[, .(date, case_rate = outcome_hat, case_lower = ci_lower, case_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "C. Valenciana, Spain",
           resolution = "subnational", metric = "case")


out <- with(jhu_global[Combined_Key == "Liaoning, China", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population, prop_zeros = 0,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 30, do_remove = FALSE,
                        family = "negbin")
)[, .(date, case_rate = outcome_hat, case_lower = ci_lower, case_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Liaoning, China",
           resolution = "global", metric = "case")

out <- with(jhu_global[Combined_Key == "Hubei, China", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population,
                        min_positive = 50, min_incr = 10, family = "negbin")
)[, .(date, case_rate = outcome_hat, case_lower = ci_lower, case_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Hubei, China",
           resolution = "global", metric = "case")

out <- with(jhu_counties[Combined_Key == "Kit Carson, Colorado", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population, prop_zeros = 0,
                        do_remove = FALSE,
                        min_positive = 50, min_incr = 2,
                        days_per_knot = 30,
                        family = "negbin")
)[, .(date, case_rate = outcome_hat, case_lower = ci_lower, case_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Kit Carson, Colorado",
           resolution = "county", metric = "case")

# Death Rate Exceptions

out <- with(jhu_counties[Combined_Key == "Anoka, Minnesota"],
            fit_poisson(date = date, new_counts = deathIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        knots = knots_touse, prop_zeros = 0,
                        family = "negbin")
)[, .(date, death_rate = outcome_hat, death_lower = ci_lower, death_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Anoka, Minnesota",
           resolution = "county", metric = "death")

out <- with(jhu_counties[Combined_Key == "Lee, Alabama"],
            fit_poisson(date = date, new_counts = deathIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        prop_zeros = 0,
                        family = "negbin", do_remove = FALSE)
)[, .(date, death_rate = outcome_hat, death_lower = ci_lower, death_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Lee, Alabama",
           resolution = "county", metric = "death")

out <- with(jhu_counties[Combined_Key == "Cumberland, Pennsylvania"],
            fit_poisson(date = date, new_counts = deathIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 60, prop_zeros = 0.3,
                        do_remove = TRUE,
                        family = "negbin")
)[, .(date, death_rate = outcome_hat, death_lower = ci_lower, death_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Cumberland, Pennsylvania",
           resolution = "county", metric = "death")

out <- with(jhu_global[Combined_Key == "Qatar"],
            fit_poisson(date = date, new_counts = deathIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        prop_zeros = 0,
                        family = "negbin", do_remove = FALSE)
)[, .(date, death_rate = outcome_hat, death_lower = ci_lower, death_upper = ci_upper)]
update_fit(fit_lst, out, loc_key = "Qatar",
           resolution = "global", metric = "death")


################################################################################
## Write out results.
################################################################################

base_fname <- paste0("../sens_analysis/results/jhu_%s_rt_case_death_rate_",si_mean_in,"_",si_sd_in,".csv")
for (fit_name in names(fit_lst)) {
  fwrite(fit_lst[[fit_name]], file = sprintf(base_fname, fit_name))
}
}

#For the Hubei parameters
run_sens_analysis(3.96, 4.75)
#For the Hong Kong parameters
run_sens_analysis(4.4, 3)
#For the multi-country parameters
run_sens_analysis(4.7, 2.9)
#For the Wuhan parameters
run_sens_analysis(7.5, 3.4)

sens_corr <- function(location){
  jhu_3.96_4.75 <- read.csv(paste0("../sens_analysis/results/jhu_",location,"_rt_case_death_rate_3.96_4.75.csv"), stringsAsFactors=FALSE)
  jhu_4.4_3 <- read.csv(paste0("../sens_analysis/results/jhu_",location,"_rt_case_death_rate_4.4_3.csv"), stringsAsFactors=FALSE)
  jhu_4.7_2.9 <- read.csv(paste0("../sens_analysis/results/jhu_",location,"_rt_case_death_rate_4.7_2.9.csv"), stringsAsFactors=FALSE)
  jhu_7.5_3.4 <- read.csv(paste0("../sens_analysis/results/jhu_",location,"_rt_case_death_rate_7.5_3.4.csv"), stringsAsFactors=FALSE)
  
  jhu_3.96_4.75 <- jhu_3.96_4.75[jhu_3.96_4.75$date %in% c('2020-04-01','2020-07-15', '2020-11-01', '2021-02-15'),c("UID","date","rt","case_rate","death_rate")]
  jhu_4.4_3 <- jhu_4.4_3[jhu_4.4_3$date %in% c('2020-04-01','2020-07-15', '2020-11-01', '2021-02-15'),c("UID","date","rt","case_rate","death_rate")]
  jhu_4.7_2.9 <- jhu_4.7_2.9[jhu_4.7_2.9$date %in% c('2020-04-01','2020-07-15', '2020-11-01', '2021-02-15'),c("UID","date","rt","case_rate","death_rate")]
  jhu_7.5_3.4 <- jhu_7.5_3.4[jhu_7.5_3.4$date %in% c('2020-04-01','2020-07-15', '2020-11-01', '2021-02-15'),c("UID","date","rt","case_rate","death_rate")]
  
  x_y <- merge(jhu_3.96_4.75, jhu_4.4_3, by=c('UID', 'date'))
  names(x_y)[3:5] <- paste0(names(x_y)[3:5], '_3.96_4.75') ; names(x_y)[6:8] <- paste0(names(x_y)[6:8], '_4.4_3')
  x_y_z <- merge(x_y, jhu_4.7_2.9, by=c('UID', 'date'))
  names(x_y_z)[9:11] <- paste0(names(x_y_z)[9:11], '_4.7_2.9')
  all_sens <- merge(x_y_z, jhu_7.5_3.4, by=c('UID', 'date'))
  names(all_sens)[12:14] <- paste0(names(all_sens)[12:14], '_7.5_3.4')
  
  date_sens_cols <- c()
  for (date_in in unique(all_sens$date)){
  all_sens_date  <- all_sens[all_sens$date == date_in,]
  rt_param_1 <- paste0(round(lm(all_sens_date$rt.x_3.96_4.75~all_sens_date$rt.y_4.4_3)$coefficients[2],3), ",", round(lm(all_sens_date$rt.x_3.96_4.75~all_sens_date$rt.y_4.4_3)$coefficients[1],3))
  rt_param_2 <- paste0(round(lm(all_sens_date$rt.x_3.96_4.75~all_sens_date$rt_4.7_2.9)$coefficients[2],3), ",", round(lm(all_sens_date$rt.x_3.96_4.75~all_sens_date$rt_4.7_2.9)$coefficients[1],3))
  rt_param_3 <- paste0(round(lm(all_sens_date$rt.x_3.96_4.75~all_sens_date$rt_7.5_3.4)$coefficients[2],3), ",", round(lm(all_sens_date$rt.x_3.96_4.75~all_sens_date$rt_7.5_3.4)$coefficients[1],3))

  date_sens_cols <- cbind( date_sens_cols, c(rt_param_1, rt_param_2, rt_param_3))
  }
  return(date_sens_cols)
}

cty <- sens_corr('county')
glb <- sens_corr('global')
stt <- sens_corr('state')
library(xtable)
xtable(rbind(cty,stt,glb))