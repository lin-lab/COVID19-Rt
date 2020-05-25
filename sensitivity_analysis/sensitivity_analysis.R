#Function for sensitivity analysis

# mean_serial_vec <- c(5.2, 7.5, 4.7, 3.96, 4.4); std_serial_vec <- c(5.1, 3.4, 2.9, 4.75, 3.0)
sensitivity_si_EpiEstim <- function(dates, positive_increase, mean_serial_vec, std_serial_vec, handle_negatives=TRUE){
  if(handle_negatives){
    positive_increase <- ifelse(positive_increase < 0, 0, positive_increase)
    warning("Warning: There are negative values provided for the daily positive increase")
  } else {
    stop("There are negative values provided for the daily positive increase")
  }
  out_r_est_list <- list()
  covid_positive <- data.frame(dates=dates,I=positive_increase)
  covid_positive$dates <- as.Date(covid_positive$dates)
  for (i in 1:length(mean_serial_vec)){
    r_est <- EpiEstim::estimate_R(covid_positive, method = "parametric_si",
                                  config = make_config(list(mean_si = mean_serial_vec[i], std_si = std_serial_vec[i]))) #, t_start=t_start, t_end=t_end)))$R
    out_r_est_list[[i]] <- r_est
  }
  return(out_r_est_list)
}

estimate_rt_EpiEstim_state <- function(state_data_combined, mean_serial, std_serial, handle_negatives=TRUE){
  out_r_est_state <- list()
  for (stateName in unique(state_data_combined$stateName)){
    state_data <- state_data_combined[state_data_combined$stateName == stateName,]
    state_data$positiveIncrease[is.na(state_data$positiveIncrease)] <- 0
    try(out_r_est_state[[stateName]] <- estimate_rt_EpiEstim(dates = state_data$date,
                                                             positive_increase = state_data$positiveIncrease,
                                                             mean_serial = mean_serial,
                                                             std_serial = std_serial,
                                                             handle_negatives = handle_negatives))
    try(out_r_est_state[[stateName]]["state"] <- stateName)
  }
  return(out_r_est_state)
}


# dates_list <- list(jhu_states$date, covidtracking_states$date, yu_states$date)
# positive_increase_list <- list(jhu_states$positiveIncrease, covidtracking_states$positiveIncrease, yu_states$positiveIncrease)
sensitivity_data_EpiEstim <- function(dates_list, positive_increase_list, mean_si = 5.2, std_si = 5.1, handle_negatives=TRUE){
  out_r_est_list <- list()
  for (i in 1:length(dates_list)){
    if(handle_negatives){
      positive_increase <- ifelse(positive_increase_list[[i]] < 0, 0, positive_increase_list[[i]])
      warning("Warning: There are negative values provided for the daily positive increase")
      positive_increase[is.na(positive_increase)] <- 0
    } else {
      stop("There are negative values provided for the daily positive increase")
    }
    
    covid_positive <- data.frame(dates=dates_list[[i]],I=positive_increase)
    covid_positive$dates <- as.Date(covid_positive$dates)
    r_est <- EpiEstim::estimate_R(covid_positive, method = "parametric_si",
                                  config = make_config(list(mean_si = mean_si, std_si = std_si))) #, t_start=t_start, t_end=t_end)))$R
    out_r_est_list[[i]] <- r_est
  }
  return(out_r_est_list)
}


estimate_rt_R0_state <- function(state_data_combined, mean_serial, std_serial, handle_negatives=TRUE){
  out_r_est_state <- list()
  for (stateName in unique(state_data_combined$stateName)){
    state_data <- state_data_combined[state_data_combined$stateName == stateName,]
    state_data$positiveIncrease[is.na(state_data$positiveIncrease)] <- 0
    try(out_r_est_state[[stateName]] <- estimate_rt_R0(dates = state_data$date,
                                                       positive_increase = state_data$positiveIncrease,
                                                       gen_mean = mean_serial,
                                                       gen_sd = std_serial,
                                                       handle_negatives = handle_negatives))
    try(out_r_est_state[[stateName]]["state"] <- stateName)
  }
  return(out_r_est_state)
}

