#Load in relevant libraries
library(EpiEstim); library(R0)

#Inputs: date vector, count of positive increases, mean serial interval, std serial interval, handling of negative positive increase
estimate_rt_EpiEstim <- function(dates, positive_increase, mean_serial, std_serial, handle_negatives=TRUE){
  if(handle_negatives){
    positive_increase <- ifelse(positive_increase < 0, 0, positive_increase)
    warning("Warning: There are negative values provided for the daily positive increase")
  } else {
    stop("There are negative values provided for the daily positive increase")
  }
  r_est <- EpiEstim::estimate_R(positive_increase, method = "parametric_si",
                                config = make_config(list(mean_si = mean_serial, std_si = std_serial)))$R #, t_start=t_start, t_end=t_end)))$R
  out_r_est <- data.frame(interval_start = dates[r_est$t_start],
                          interval_end = dates[r_est$t_end],
                          mean_rt = r_est$`Mean(R)`,
                          std_rt = r_est$`Std(R)`,
                          ci_lower = r_est$`Quantile.0.025(R)`,
                          ci_upper = r_est$`Quantile.0.975(R)`)
  return(out_r_est)
}


#Inputs: date vector, count of positive increases, mean serial interval, std serial interval, handling of negative positive increase
estimate_rt_R0 <- function( dates, positive_increase, gen_mean, gen_sd, handle_negatives=TRUE ){
  if(handle_negatives){
    positive_increase <- ifelse(positive_increase < 0, 0, positive_increase)
    warning("Warning: There are negative values provided for the daily positive increase")
  } else {
    stop("There are negative values provided for the daily positive increase")
  }
  names(positive_increase) <- dates
  gen_time <- generation.time("gamma", c(gen_mean, gen_sd)) # from DOI 10.7326/M20-0504 - update to consider weibull, etc.
  r_est <- estimate.R(epid=positive_increase, GT=gen_time,
                       methods=c("TD"), begin=1L, end=as.integer(length(positive_increase)))
  out_r_est <- data.frame(unname(cbind(names(r_est$estimates$TD$R), r_est$estimates$TD$R, r_est$estimates$TD$conf.int)))
  names(out_r_est) <- c('date', 'rt_est', 'rt_ci_lower', 'rt_ci_upper')
  return(out_r_est)
}
