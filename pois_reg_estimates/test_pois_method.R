library(splines)
library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(logger)
library(dplyr)

log_threshold(TRACE)

source("../preprocess_data/preprocess_jhu.R")
source("../helpers/helper_functions.R")
source("pois_functions.R")

click_plot <- function(plt_dat, place_name) {
  rt_plt_title <- sprintf("Rt for %s", place_name)
  rt_plt <- plt_dat %>%
    filter(!is.na(rt)) %>%
    ggplot(aes(x = date, y = rt, ymin = rt_lower, ymax = rt_upper)) +
    geom_ribbon(fill = "#9e9e9e") + geom_line() +
    #coord_cartesian(ylim = c(0, ymax_rt)) +
    geom_hline(yintercept = 1, lty = 2) +
    xlab("Date") + ylab("") + ggtitle(rt_plt_title) +
    theme_cowplot() +
    background_grid(major = "xy", minor = "xy") +
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 15))

  newcases_plt_title <- sprintf("New Cases per Million for %s", place_name)
  ymax_newcases <- NA
  newcases_plt <- plt_dat %>%
    ggplot(aes(x = date, y = case_rate, ymin = case_lower, ymax = case_upper)) +
    #ggplot(aes(x = date)) +
    geom_ribbon(fill = "#9e9e9e") + geom_line(aes(linetype = "Smoothed")) +
    geom_line(aes(y = positiveIncrease_percapita, linetype = "Unsmoothed")) +
    xlab("Date") + ylab("") + ggtitle(newcases_plt_title) +
    theme_cowplot() +
    background_grid(major = "xy", minor = "xy") +
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 15),
          legend.position = "bottom") +
    guides(linetype = guide_legend(nrow = 2))

  deaths_plt_title <- sprintf("New Deaths per Million for %s", place_name)
  deaths_plt <- plt_dat %>%
    ggplot(aes(x = date, y = death_rate, ymin = death_lower, ymax = death_upper)) +
    geom_ribbon(fill = "#9e9e9e") + geom_line(aes(linetype = "Smoothed")) +
    geom_line(aes(y = deathIncrease_percapita, linetype = "Unsmoothed")) +
    xlab("Date") + ylab("") + ggtitle(deaths_plt_title) +
    theme_cowplot() +
    coord_cartesian(ylim = c(0, NA)) +
    background_grid(major = "xy", minor = "xy") +
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 15),
          legend.position = "bottom") +
    guides(linetype = guide_legend(nrow = 2))
  final_plt <- plot_grid(rt_plt, newcases_plt, deaths_plt, ncol = 1,
                         align = "v", axis = "l")
  return(final_plt)
}

reformat_data <- function(dt) {
  scale_cols <- c("case_rate", "case_lower", "case_upper", "death_rate",
                  "death_lower", "death_upper")
  dt[, (scale_cols) := lapply(.SD, function(x) { x * 1e6 }),
     .SDcols = scale_cols]
  dt[, positive_percapita := 1e6 * positive / population]
  dt[, death_percapita := 1e6 * death / population]
  dt[, positiveIncrease_percapita := 1e6 * positiveIncrease / population]
  dt[, deathIncrease_percapita := 1e6 * deathIncrease / population]
}

plot_out <- function(out) {
  out %>%
    dplyr::filter(outcome_hat < Inf) %>%
    ggplot(aes(x = date, y = outcome_hat)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#9e9e9e") +
    geom_line() + geom_point() +
    #geom_hline(yintercept = 1, lty = 2) +
    background_grid(major = "xy", minor = "xy") +
    theme(text = element_text(size = 18),
          axis.text = element_text(size = 15))
}

# read files
days_lag = 7L
start_date <- as.Date("2020-03-01")
jhu_global <- read_jhu("Global", start_date = start_date)
setkey(jhu_global, UID)

source("pois_functions.R")
problem_uids <- c(558, 728, 140)
problem_df <- jhu_global[UID %in% problem_uids]
test_fit <- fit_poisson_by_unit(problem_df,
                                family = "negbin",
                                min_incr = 0,
                                days_per_knot = 30, min_positive = 50) %>%
  reformat_data() %>%
  lag_dates(lag = days_lag)

plt_lst <- lapply(problem_uids,
                  function(uid) {
                    dat <- test_fit[UID == uid, ]
                    country <- problem_df[UID == uid, Combined_Key][1]
                    click_plot(dat, country)
                  })
                    
plot_grid(plotlist = plt_lst, nrow = 1)

nic_dates <- jhu_global[Country_Region == "Nicaragua", date]
nic_counts <- jhu_global[Country_Region == "Nicaragua", positiveIncrease]
nic_pop <- jhu_global[Country_Region == "Nicaragua", population]
out <- fit_poisson(nic_dates, nic_counts, min_positive = 50, min_incr = 20)
plot_out(out)
out <- fit_poisson(nic_dates, nic_counts, min_positive = 50, min_incr = 1,
                   population = nic_pop)
plot_out(out)



start_date <- as.Date("2020-03-01")
jhu_counties <- read_jhu("County", start_date = start_date)[!(startsWith(county, "Out of") | startsWith(county, "Unassigned")), ]
jhu_counties[, Combined_Key := paste0(county, ", ", stateName)]
jhu_states <- read_jhu("State", start_date = start_date)
jhu_states[, Combined_Key := stateName]
jhu_global <- read_jhu("Global", start_date = start_date)
jhu_subnational_orig <- read_jhu("Subnational", start_date = start_date)[Province_State != "Unknown", ]
jhu_subnational <- jhu_subnational_orig %>%
  group_by(UID) %>%
  filter(date > min(date)) %>%
  data.table()

source("pois_functions.R")
#problem_uids <- c(316, 84000002, 84000050)
problem_uids <- c(84000009)
problem_states <- jhu_states[UID %in% problem_uids, ]
test_fit <- fit_poisson_by_unit(problem_states,
                                min_incr = 0,
                                days_per_knot = 30, min_positive = 100,
                                family = "negbin") %>%
  reformat_data() %>%
  lag_dates(days_lag)

test_fit %>%
  click_plot('Connecticut')

source("pois_functions.R")
knots_touse <- c(seq(as.Date("2020-03-30"), max(jhu_states$date), by = 30),
                 as.Date("2020-12-07")) %>%
  sort()
tmp <- jhu_states[UID == 84000022, ]
tmp
out <- with(tmp, fit_poisson(date = date, new_counts = positiveIncrease,
                             min_positive = 50, min_incr = 1,
                             knots = knots_touse, family = "negbin")
)
out


new_counts <- jhu_states[UID == 84000050, deathIncrease]
date <- jhu_states[UID == 84000050, date]
population <- jhu_states[UID == 84000050, population]
out <- fit_poisson(date, new_counts, population,
                   lagged_counts = FALSE, min_positive = 50, min_incr = 1,
                   family = "negbin")

new_counts <- jhu_states[UID == 84000009, positiveIncrease]
date <- jhu_states[UID == 84000009, date]
population <- jhu_states[UID == 84000009, population]
out <- fit_poisson(date, new_counts, population,
                   lagged_counts = FALSE, min_positive = 50, min_incr = 0,
                   family = "negbin")

providence_ri <- jhu_counties[Combined_Key == "Providence, Rhode Island"]
new_counts <- providence_ri[, positiveIncrease]
date <- providence_ri[, date]
out <- fit_poisson(date, new_counts,
                   lagged_counts = FALSE, min_positive = 50, min_incr = 10,
                   family = "negbin")

problem_county <- jhu_counties[Combined_Key == "Carbon, Wyoming"]
out <- with(problem_county,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        knots = knots_touse, family = "negbin")
)


options(warn = 2)
china_dat <- jhu_global[UID == 156, ]
new_counts <- china_dat[, positiveIncrease]
date <- china_dat[, date]
out <- fit_poisson(date, new_counts,
                   lagged_counts = FALSE, min_positive = 50, min_incr = 20,
                   family = "negbin")

source("pois_functions.R")
test_fit <- fit_poisson_by_unit(china_dat, lagged_counts = FALSE,
                                min_incr = 20,
                                days_per_knot = 30, min_positive = 50,
                                family = "negbin") %>%
  reformat_data() %>%
  lag_dates()

hubei_dat <- jhu_global[Combined_Key == "Hubei, China", ]
hubei_dat
out <- with(hubei_dat,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population,
                        min_positive = 50, min_incr = 10, family = "negbin")
)

tmp <- jhu_subnational[Combined_Key == "England, United Kingdom", ]
tmp
setorderv(tmp, "date")
out <- with(tmp,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        days_per_knot = 30,
                        min_positive = 50, min_incr = 20, family = "negbin",
                        knots = NULL)
)

tmp <- jhu_subnational[Combined_Key == "Northern Ireland, United Kingdom", ]
tmp
setorderv(tmp, "date")
out <- with(tmp,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        days_per_knot = 30,
                        min_positive = 50, min_incr = 20, family = "negbin",
                        knots = NULL)
)

out <- with(tmp,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population,
                        days_per_knot = 30,
                        min_positive = 50, min_incr = 20, family = "negbin",
                        knots = NULL)
)

problem_county <- jhu_counties[Combined_Key == "Okfuskee, Oklahoma"]
out <- with(problem_county,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        min_positive = 50, min_incr = 1, do_remove = FALSE,
                        prop_zeros = 0,
                        knots = knots_touse, family = "negbin")
)
plot_out(out)


problem_county <- jhu_counties[Combined_Key == "Anoka, Minnesota"]
out <- with(problem_county,
            fit_poisson(date = date, new_counts = deathIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        knots = knots_touse, prop_zeros = 0,
                        family = "negbin")
)
plot_out(out)


problem_county <- jhu_counties[Combined_Key == "Lee, Alabama"]
out <- with(problem_county,
            fit_poisson(date = date, new_counts = deathIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        prop_zeros = 0,
                        family = "negbin", do_remove = FALSE)
)
plot_out(out)

problem_county <- jhu_counties[Combined_Key == "Rusk, Texas"]
out <- with(problem_county,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        min_positive = 50, min_incr = 20,
                        knots = knots_touse, do_remove = TRUE,
                        prop_zeros = 0,
                        family = "negbin")
)
plot_out(out)

problem_county <- jhu_counties[Combined_Key == "Gillespie, Texas"]
out <- with(problem_county,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        min_positive = 50, min_incr = 20,
                        knots = knots_touse,
                        do_remove = FALSE, prop_zeros = 0,
                        family = "negbin")
)
plot_out(out)

problem_county <- jhu_counties[Combined_Key == "Lincoln, Arkansas"]
out <- with(problem_county,
            fit_poisson(date = date, new_counts = positiveIncrease,
                        min_positive = 50, min_incr = 20,
                        days_per_knot = 75, prop_zeros = 0.3,
                        do_remove = TRUE,
                        family = "negbin")
)
plot_out(out)

out <- with(jhu_counties[Combined_Key == "Cumberland, Pennsylvania"],
            fit_poisson(date = date, new_counts = deathIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 60, prop_zeros = 0.3,
                        do_remove = TRUE,
                        family = "negbin")
)
plot_out(out)

out <- with(jhu_counties[Combined_Key == "Jefferson, Washington"],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 60, prop_zeros = 0.3,
                        do_remove = TRUE,
                        family = "negbin")
)
plot_out(out)

out <- with(jhu_counties[Combined_Key == "Madison, Texas"],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        min_positive = 50, min_incr = 20,
                        days_per_knot = 30, prop_zeros = 0.3,
                        do_remove = TRUE,
                        family = "negbin")
)
plot_out(out)

out <- with(jhu_subnational[Combined_Key == "Baleares, Spain", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population, prop_zeros = 0,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 30,
                        family = "negbin")
)
plot_out(out)

out <- with(jhu_subnational[Combined_Key == "C. Valenciana, Spain", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population, prop_zeros = 0.2,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 30,
                        family = "negbin")
)
plot_out(out)

out <- with(jhu_global[Combined_Key == "Liaoning, China", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        population = population, prop_zeros = 0.2,
                        min_positive = 50, min_incr = 1,
                        days_per_knot = 30,
                        family = "negbin")
)
plot_out(out)

out <- with(jhu_global[Combined_Key == "Madagascar", ],
            fit_poisson(date = date, new_counts = positiveIncrease,
                        prop_zeros = 0.2,
                        min_positive = 50, min_incr = 20,
                        days_per_knot = 30,
                        family = "negbin")
)
plot_out(out)
