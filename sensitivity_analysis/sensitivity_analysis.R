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
    out_r_est <- r_est
    out_r_est_list[[i]] <- out_r_est
  }
  return(out_r_est_list)
}

estimate_rt_EpiEstim_state <- function(state_data_combined, mean_serial, std_serial, handle_negatives=TRUE){
  out_r_est_state <- list()
  for (stateName in unique(state_data_combined$stateName)){
    state_data <- state_data_combined[state_data_combined$stateName == stateName,]
    try(out_r_est_state[[stateName]] <- estimate_rt_EpiEstim(dates = state_data$date,
                                                             positive_increase = state_data$positiveIncrease,
                                                             mean_serial = mean_serial,
                                                             std_serial = std_serial,
                                                             handle_negatives = handle_negatives))
    try(out_r_est_state[[stateName]]["state"] <- stateName)
  }
  return(out_r_est_state)
}

