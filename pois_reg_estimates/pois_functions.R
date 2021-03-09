library(splines)
library(data.table)
library(logger)
library(lubridate)
library(MASS)
library(glm2)

read_jhu <- function(level = c("County", "State", "Global", "Subnational"),
                     start_date) {
  level <- match.arg(level)
  base_url <- "https://hsph-covid-study.s3.us-east-2.amazonaws.com/JHU_Cleaned/"
  filename <- sprintf("JHU_COVID-19_%s.csv", level)
  zip_filename <- paste0(filename, ".zip")
  data_url <- paste0(base_url, zip_filename)

  # download AWS zipped csv to temp file, and read it in
  temp_file <- tempfile()
  download.file(data_url, temp_file)
  cmd <- sprintf("unzip -p %s", temp_file)
  dat <- data.table::fread(cmd = cmd)[date >= start_date, ]
  dat[, date := as.Date(date, format = "%Y-%m-%d")]
  return(dat)
}

# Weighted shift to "correct" for random lags
random_lag <- function(y, lags, weights) {
  stopifnot(length(weights) == length(lags))
  stopifnot(min(weights) >= 0)
  out <- rep(0, length(y))

  weights <- weights / sum(weights)

  for (i in seq_along(weights)) {
    if (weights[i] > 0) {
      out <- out + weights[i] * data.table::shift(y, n = -lags[i], fill = NA)
    }
  }
  return(out)
}

# Weighted shift to "correct" for random Gamma-distributed lags
gamma_lag <- function(y, max_days = 14, shape = 5, rate = 1) {
  l <- 0:max_days
  w <- pgamma(l + 1, shape, rate) - pgamma(l, shape, rate)
  w <- w / sum(w)
  return(random_lag(y, l, w))
}

# generate model matrix for weekday adjustment
weekday_matrix <- function(dates) {
  model.matrix(~weekdays(dates))[, -1]
}

warning_handler <- function(w) {
  if (any(grepl("fitted rates numerically 0 occurred", w))) {
    log_warn(sprintf("Ignoring warning: %s", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
}

call_glm <- function(formula, data,
                     family = c("quasipoisson", "negbin", "poisson"),
                     chisq_cutoff = 0.01, refit = TRUE) {
  family <- match.arg(family)
  fit_pois <- (family != "negbin")
  if (family == "negbin") {
    ret <- NULL
    ret <- tryCatch({
      # use more iterations and better fitter
      withCallingHandlers(glm.nb(formula, data = data, method = "glm.fit2",
                                 control = glm.control(maxit = 1000)),
                          warning = warning_handler)
      },
      warning = function(w) {
        log_warn(sprintf("Warning in glm.nb: %s", conditionMessage(w)))
        NULL
      },
      error = function(e) {
        log_warn(sprintf("Error in glm.nb: %s", conditionMessage(e)))
        NULL
      }
    )
    if (is.null(ret) && refit) {
      fit_pois <- TRUE
      log_warn("NB fit returned error. Refitting with poisson.")
    } else if (!is.null(ret)) {
      log_trace(sprintf("NB fit success: theta = %0.3f", ret$theta))
    } else {
      log_warn("NB fit returned error.")
    }
  }

  if (fit_pois) {
    # try Poisson first
    ret <- withCallingHandlers({
      stats::glm(formula, data = data, family = "poisson")
    }, warning = warning_handler)

    # check for overdispersion---if none, stick with Poisson
    overdisp_stat <- pchisq(deviance(ret), df.residual(ret), lower.tail = FALSE)

    if (overdisp_stat < chisq_cutoff) {
      log_info(sprintf("Overdispersion detected. %f deviance on %d df, pchisq = %e. Using quasipoisson...",
                       deviance(ret), df.residual(ret), overdisp_stat))
      ret <- withCallingHandlers({
        stats::glm(formula, data = data, family = "quasipoisson")
      }, warning = warning_handler)
    } else {
      log_trace(sprintf("No overdispersion detected. %f deviance on %d df, pchisq = %e.",
                        deviance(ret), df.residual(ret), overdisp_stat))
    }
  }
  return(ret)
}

fit_spline_glm <- function(formula_str, data, days_per_knot = 30,
                           knots = NULL, refit = FALSE,
                           family = c("quasipoisson", "negbin", "poisson")) {
  family <- match.arg(family)
  if (is.null(knots)) {
    log_info(sprintf("Using knots every %d days.", days_per_knot))
    formula_splines <-
      paste(formula_str,
            sprintf("+ bs(date, knots = seq_days(date, %d))", days_per_knot))
    model_formula <- as.formula(formula_splines)
  } else {
    log_info(sprintf("Using custom knots: %s", paste(knots, collapse = "; ")))
    model_formula <- update.formula(as.formula(formula_str), ~ . + bs(date, knots = knots))
  }

  fit <- tryCatch({
    call_glm(model_formula, data, family = family, refit = refit)
    },
    warning = function(w) {
      # usually will get a warning if GLM didn't converge
      log_warn("Handled warning in GLM fitting: ", conditionMessage(w))
      NULL
    },
    error = function(e) {
      log_warn("Handled error in GLM fitting: ", conditionMessage(e))
      NULL
    }
  )
  knots_str <- ifelse(is.null(knots), as.character(days_per_knot), "custom")
  if (is.null(fit)) {
    return(NULL)
  }
  fit_family <- fit$family$family
  if (startsWith(fit_family, "Negative Binomial")) {
    log_info(sprintf("Fit successfully with NB: %s knots; theta = %f",
                     knots_str, fit$theta))

  } else if (fit_family == "quasipoisson") {
    disp <- summary(fit)$dispersion
    log_info(sprintf("Fit successfully with quasipoisson: %s knots; phi = %f",
                     knots_str, disp))
  } else {
    log_info(sprintf("Fit successfully with poisson: %s knots", knots_str))
  }
  return(fit)
}

remove_outliers <- function(full_data) {
  data_copy <- copy(full_data)
  data_copy[, rollmean_counts := lagMean(new_counts, 7)]
  idx <- (with(data_copy, new_counts > 5 * rollmean_counts))
  if (sum(idx) < 10) {
    log_info(sprintf("Removing counts on %d days: %s", sum(idx),
                    paste(data_copy[idx, date], collapse = "; ")))

    data_copy <- data_copy[!idx, ]
  }
  setorderv(data_copy, cols = "new_counts", order = -1L)
  while(nrow(data_copy) > 2 &&
        data_copy$new_counts[1] > 2 * data_copy$new_counts[2]) {
    log_info(sprintf("Removing count %d on %s",
             data_copy$new_counts[1], format(data_copy$date[1], "%Y-%m-%d")))
    data_copy$new_counts[1] <- NA
    data_copy <- na.omit(data_copy)
  }
  return(data_copy)
}

#' Estimate Rt or case/death rates using Poisson model with splines
#'
#' Estimate Rt or case/death rates by fitting a Poisson model to daily new cases
#' at a particular location.
#'
#' @param date Vector of dates
#' @param new_counts Vector of new counts, one for each date.
#' @param population Population for the location of interest where the new
#' counts occurred. Should either be length 1 or have length equal to new_counts
#' and date. If supplied, the case/death rate will be calculated. If not
#' supplied (default), Rt will be calculated.
#' @param days_per_knot For the spline, how many days apart should each knot be
#' spaced?
#' @param min_positive Don't calculate Rt for dates when the cumulative number
#' of cases is below this number.
#' @param min_incr Don't calculate Rt on dates when the average number of cases
#' in the past 7 days is below this number.
#' @param prop_zeros If the proportion of zeros in the data (after removing
#' outliers) is above this threshold, then calculate a weekly Rt.
#' @param family The GLM family to be fitted.
#' @param SI_mean Mean of serial interval distribution
#' @param SI_sd Standard deviation of serial interval distribution
#' @param SI_NTS Max number of days back that a primary case can infect a
#' secondary case.
#' @param min_days Don't fit the Poisson model if the number of dates observed
#' is less than this number.
#' @param conf_level Confidence level for confidence intervals.
#' @param knots If non-null, use these dates as knots instead of the default
#' knots.
fit_poisson <- function(date, new_counts = NULL,
                        population = NULL, days_per_knot = 30,
                        min_positive = 50, min_incr = 10,
                        prop_zeros = 0.2,
                        family = c("negbin", "quasipoisson", "poisson"),
                        SI_mean = 5.2, SI_sd = 5.5, SI_NTS = 30,
                        min_days = 14, do_remove = TRUE,
                        conf_level = 0.95, knots = NULL) {

  stopifnot(0 < conf_level && conf_level < 1)
  stopifnot(SI_mean > 0 && SI_sd > 0)
  stopifnot(SI_NTS > 0)
  stopifnot(days_per_knot > 0)
  stopifnot(length(date) == length(new_counts))
  family <- match.arg(family)
  days_per_knot_use <- days_per_knot

  if (min(new_counts) < 0) {
    log_warn(sprintf("Set %d/%d negative case counts to 0.\n",
                     sum(new_counts < 0), length(new_counts)))
    new_counts[new_counts < 0] <- 0
  }

  if (is.null(population)) {
    full_data <- data.table(new_counts = new_counts, date = date)
  } else {
    full_data <- data.table(new_counts = new_counts, date = date,
                            population = population)
  }
  if (do_remove) {
    full_data <- remove_outliers(full_data)
  }

  data.table::setorder(full_data, date)
  out_error <- data.table::data.table(date = date,
                                      outcome_hat = NA, ci_lower = NA,
                                      ci_upper = NA, offset_var = NA)
  if (nrow(full_data) == 0) {
    log_info("All data were outliers, returning error...")
    return(out_error)
  }

  if (nrow(full_data) <= min_days) {
    log_warn("Too few days, Poisson GLM not fitted.")
    return(out_error)
  }

  full_data[, `:=`(
               total_counts = cumsum(new_counts),
               rollmean_counts = lagMean(new_counts, 7),
               week_num = paste0(lubridate::year(date), "_",
                                 lubridate::week(date)))]
  # set offset variable: either Lambda from EpiEstim or log(population)
  if (is.null(population)) {
    full_data[, offset_var := log1p(getLambda(
      x = new_counts,
      t = as.integer(date),
      mu = SI_mean, sd = SI_sd, NTS = SI_NTS))]
    full_data[, offset_orig := exp(offset_var) - 1]
  } else {
    full_data[, `:=`(offset_var = log(population),
                     offset_orig = population)]
  }

  week_data <- full_data[,
                         .(sum_counts = sum(new_counts), week_start = min(date),
                           week_end = max(date)), by = week_num]
  if (is.null(population)) {
    # set offset to be previous week's count
    week_data[, offset_orig := as.numeric(data.table::shift(sum_counts, n = 1, type = "lag"))]
    week_data[, offset2 := (sum_counts + shift(offset_orig, n = 1, type = "lag")) / 2]

    # if offset is 0, make it the avg of current week and 2 weeks ago
    week_data[offset_orig == 0, offset_orig := offset2]
    week_data[, offset_var := log(offset_orig)]
  } else {
    # set offset to be population
    pop_weekly <- full_data[, .(pop = mean(population)), by = week_num]
    week_data[, offset_orig := pop_weekly$pop]
    week_data[, offset_var := log(pop_weekly$pop)]
  }

  data_subset <- full_data[total_counts >= min_positive &
                           rollmean_counts >= min_incr &
                           offset_var > 0, ]

  if (nrow(data_subset) <= min_days) {
    log_warn("Too few days, Poisson GLM not fitted.")
    return(out_error)
  }
  # check if data is irregularly reported
  irregular_dat <- isTRUE(nrow(data_subset == 0)) ||
    isTRUE(mean(data_subset$new_counts <= 0, na.rm = TRUE) >= prop_zeros)
  log_trace(sprintf("Proportion of zeros: %0.3f",
                    mean(data_subset$new_counts <= 0, na.rm = TRUE)))

  # data_fit is either the weekly data or the subsetted full data.
  if (irregular_dat) {
    log_info("Irregular data, using weekly counts instead.")
    if (nrow(data_subset) == 0) {
      start_date <- min(full_data$date)
      end_date <- max(full_data$date)
    } else {
      start_date <- min(data_subset$date)
      end_date <- max(data_subset$date)
    }
    # start time 0 as date when there are over min_positive cases
    row_start <- 1 + with(week_data, which(week_start <= start_date &
                                           start_date <= week_end))
    row_end <- with(week_data, which(week_start <= end_date &
                                     end_date <= week_end))
    stopifnot(length(row_start) == 1)
    stopifnot(length(row_end) == 1)
    if (row_start >= row_end) {
      log_warn(sprintf("Start date %s and end date %s are in same week, no fitting performed.",
                       start_date, end_date))
      return(out_error)
    }
    data_fit <- week_data[row_start:row_end, ]
    data_fit[, date := week_start]
    formula_str <- "sum_counts ~ offset(offset_var)"
    if (days_per_knot < 60) {
      log_info("Setting days_per_knot to 60...")
      days_per_knot_use <- 60
    }
    knots_touse <- NULL
  } else {
    if (nrow(data_subset) <= min_days) {
      log_warn("Too few days, Poisson GLM not fitted.")
      return(out_error)
    }
    formula_str <- "new_counts ~ offset(offset_var)"
    data_fit <- data_subset
    knots_touse <- knots
  }


  if (nrow(data_fit) <= 5) {
    log_warn(sprintf("Data only has %d observations, no fitting performed.",
                     nrow(data_fit)))
    return(out_error)
  }

  fit <- fit_spline_glm(formula_str, data_fit,
                        days_per_knot = days_per_knot_use, family = family,
                        knots = knots_touse, refit = (days_per_knot_use >= 120))

  days_per_knot_use <- days_per_knot_use + 15
  while(days_per_knot_use <= 120 &&
        (is.null(fit) || length(fit$coefficients) > fit$rank)) {
    # sometimes we'll get a rank-deficient fit. In that case, try fitting with
    # more days per knot.
    do_refit <- days_per_knot_use >= 120
    log_warn(sprintf("Rank deficient fit. Fitting with %d days per knot instead",
                     days_per_knot_use))
    fit <- fit_spline_glm(formula_str, data_fit,
                          days_per_knot = days_per_knot_use, family = family,
                          refit = do_refit)
    days_per_knot_use <- days_per_knot_use + 15
  }

  if (is.null(fit) || !fit$converged) {
    log_warn("Returning error...")
    return(out_error)
  }

  # run prediction to get estimates and confidence intervals
  pred <- tryCatch({
      stats::predict.glm(fit, newdata = data_fit, type = "link", se.fit = TRUE)
    },
    warning = function(w) {
      log_warn("Handled warning in predict.glm: ", conditionMessage(w))
      NULL
    },
    error = function(e) {
      log_warn("Handled error in predict.glm: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(pred)) {
    log_warn("Returning error...")
    return(out_error)
  }

  # data table of predictions
  zstar <- qnorm(1 - (1 - conf_level) / 2)
  res_dt <- data.table(date = data_fit$date,
                       outcome_hat = exp(pred$fit) / data_fit$offset_orig,
                       ci_lower = exp(pred$fit - zstar * pred$se.fit) /
                         data_fit$offset_orig,
                       ci_upper = exp(pred$fit + zstar * pred$se.fit) /
                         data_fit$offset_orig,
                       offset_var = data_fit$offset_var)

  if (irregular_dat) {
    if (!is.null(population)) {
      # make weekly case/death rates into daily
      res_dt[, outcome_hat := outcome_hat / 7]
      res_dt[, ci_lower := ci_lower / 7]
      res_dt[, ci_upper := ci_upper / 7]
    }
    # make weekly predictions into daily predictions.
    res_dt[, week_num := paste0(lubridate::year(date), "_",
                                lubridate::week(date))]
    df1 <- full_data[, .(date, week_num)]
    df2 <- res_dt[, !"date"]
    # merge on week_num and drop week_num column
    out <- df1[df2, on = "week_num"][, !"week_num"]
  } else {
    out <- res_dt
  }
  return(out)
}

fit_rt_case_death <- function(cur_uid, data, rt_opts, case_opts, death_opts) {
  data_subset <- data[UID == cur_uid]
  log_info(sprintf("Working on UID %d (%s)...", cur_uid,
                   data_subset$Combined_Key[1]))

  log_trace("Working on Rt regression...")
  rt_fit_params <- c(list(date = data_subset$date,
                          new_counts = data_subset$positiveIncrease,
                          population = NULL), rt_opts)
  rt_fit <- do.call(fit_poisson, rt_fit_params)[,
                    .(date, rt = outcome_hat, rt_lower = ci_lower,
                      rt_upper = ci_upper)]

  log_trace("Working on new cases regression...")
  case_fit_params <- c(list(date = data_subset$date,
                            new_counts = data_subset$positiveIncrease,
                            population = data_subset$population), case_opts)
  case_fit <- do.call(fit_poisson, case_fit_params)[,
                      .(date, case_rate = outcome_hat,
                        case_lower = ci_lower, case_upper = ci_upper)]

  log_trace("Working on new deaths regression...")
  death_fit_params <- c(list(date = data_subset$date,
                             new_counts = data_subset$deathIncrease,
                             population = data_subset$population), death_opts)
  death_fit <- do.call(fit_poisson, death_fit_params)[,
                        .(date, death_rate = outcome_hat,
                          death_lower = ci_lower, death_upper = ci_upper)]

  tmp <- merge(rt_fit, case_fit, by = "date", all = TRUE)
  fit_df <- merge(tmp, death_fit, by = "date", all = TRUE)
  fit_df$UID <- cur_uid
  return(fit_df)
}

fit_poisson_by_unit <- function(data, days_per_knot = 30, min_positive = 50,
                                min_incr = 10,
                                family = c("negbin", "quasipoisson", "poisson"),
                                SI_mean = 5.2, SI_sd = 5.5, SI_NTS = 30,
                                min_days = 14, knots = NULL) {

  uniq_uids <- unique(data$UID)
  n_uids <- length(uniq_uids)
  data.table::setkey(data, UID, date)

  rt_opts <- list(days_per_knot = days_per_knot,
                  min_positive = min_positive,
                  min_incr = min_incr,
                  family = family,
                  SI_mean = SI_mean,
                  SI_sd = SI_sd,
                  SI_NTS = SI_NTS,
                  min_days = min_days,
                  knots = knots)

  case_opts <- list(days_per_knot = days_per_knot,
                    min_positive = 50,
                    min_incr = 1,
                    family = family,
                    SI_mean = SI_mean,
                    SI_sd = SI_sd,
                    SI_NTS = SI_NTS,
                    min_days = min_days,
                    knots = knots)

  death_opts <- case_opts

  fit_df_lst <- lapply(X = uniq_uids, FUN = fit_rt_case_death, data = data,
                       rt_opts = rt_opts, case_opts = case_opts,
                       death_opts = death_opts)

  all_fitted <- do.call(rbind, fit_df_lst)
  all_fitted[, date_str := format(date, "%Y-%m-%d")]
  all_fitted[, date := NULL]
  data[, date_str := format(date, "%Y-%m-%d")]
  data.table::setkey(data, UID, date_str)
  data.table::setkey(all_fitted, UID, date_str)

  ret <- merge(data, all_fitted, all = TRUE)
  ret[, date_str := NULL]
  return(ret)
}

# Adjusts Rt by a fixed lag. Performs this in-place on a data.table.
lag_dates <- function(dt, lag) {
  stopifnot("date" %in% colnames(dt))
  stopifnot(lag > 0)
  stopifnot(is.integer(lag))
  dt[, date_lag := date - as.difftime(lag, units = "days")]
}
