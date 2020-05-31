#Load in relevant libraries
library(readr)
library(EpiEstim)
library(R0)
library(corrplot)
library(RColorBrewer)
source("./COVID19-Rt/sensitivity_analysis/ggcorplot.R")
source('./COVID19-Rt/sensitivity_analysis/sensitivity_analysis.R')
source('./COVID19-Rt/preprocess_data/preprocess_jhu.R')
source('./COVID19-Rt/preprocess_data/preprocess_covidtracking.R')
source('./COVID19-Rt/preprocess_data/preprocess_nyt.R')
source('./COVID19-Rt/preprocess_data/preprocess_yu.R')
source('./COVID19-Rt/estimate_rt/estimate_rt_master.R')


######## JHU state
jhu_states_pre <- load_jhu(level='State',pull_date = '2020-05-12', start_date = '2020-03-18', end_date = '2020-05-08')
non_us_stateName <- unique(jhu_states_pre$stateName)[c(3,10,14,15,40,45,54)] # non US states
jhu_states_pre <- jhu_states_pre[jhu_states_pre$stateName %in% setdiff(unique(jhu_states_pre$stateName), non_us_stateName),]

jhu_states_post <- load_jhu(level='State',pull_date = '2020-05-12', start_date = '2020-03-25', end_date = '2020-05-08')
non_us_stateName <- unique(jhu_states_post$stateName)[c(3,10,14,15,40,45,54)] # non US states
jhu_states_post <- jhu_states_post[jhu_states_post$stateName %in% setdiff(unique(jhu_states_post$stateName), non_us_stateName),]

##### load Lancet Inf Dis calcuated R
state_result_Lancet_combined <- read.csv('./COVID19-Rt/sensitivity_analysis/Lancet_Inf_Dis/Rt_Lancet_moving_average.csv')
state_result_Lancet_combined$date <- as.Date(state_result_Lancet_combined$date)
state_result_Lancet_combined <- state_result_Lancet_combined[state_result_Lancet_combined$date >= "2020-03-25",]

## Rt plot
for (stateName in unique(jhu_states_pre$stateName)){
  print(stateName)
  state_data_pre <- jhu_states_pre[jhu_states_pre$stateName == stateName,]
  state_data_pre$positiveIncrease[is.na(state_data_pre$positiveIncrease)] <- 0
  state_result_EpiEstim <- estimate_rt_EpiEstim(dates = state_data_pre$date,
                                                positive_increase = state_data_pre$positiveIncrease,
                                                mean_serial = 5.2,
                                                std_serial = 5.1,
                                                handle_negatives = TRUE)
  state_result_EpiEstim <- state_result_EpiEstim[,c("interval_end","mean_rt","ci_lower","ci_upper")]
  colnames(state_result_EpiEstim) <- c("date", "R", "Rt_lower", "Rt_upper")
  state_result_EpiEstim$date <- as.Date(state_result_EpiEstim$date)
  state_result_EpiEstim$R <- as.numeric(as.character(state_result_EpiEstim$R))
  state_result_EpiEstim$Rt_lower <- as.numeric(as.character(state_result_EpiEstim$Rt_lower))
  state_result_EpiEstim$Rt_upper <- as.numeric(as.character(state_result_EpiEstim$Rt_upper))
  state_result_EpiEstim$Rt_method <- "EpiEstim"
  
  state_data_post <- jhu_states_post[jhu_states_post$stateName == stateName,]
  state_data_post$positiveIncrease[is.na(state_data_post$positiveIncrease)] <- 0
  state_result_R0 <- estimate_rt_R0(dates = state_data_post$date,
                                    positive_increase = state_data_post$positiveIncrease,
                                    gen_mean = 5.2,
                                    gen_sd = 5.1,
                                    handle_negatives = TRUE)
  rownames(state_result_R0) <- NULL
  colnames(state_result_R0) <- c("date", "R", "Rt_lower", "Rt_upper")
  state_result_R0$date <- as.Date(state_result_R0$date)
  state_result_R0$R <- as.numeric(as.character(state_result_R0$R))
  state_result_R0$Rt_lower <- as.numeric(as.character(state_result_R0$Rt_lower))
  state_result_R0$Rt_upper <- as.numeric(as.character(state_result_R0$Rt_upper))
  state_result_R0$Rt_method <- "TD"
  
  state_result_Lancet <- state_result_Lancet_combined[state_result_Lancet_combined$stateName == stateName,]
  state_result_Lancet <- state_result_Lancet[,c("date","Rt","Rt_0.025","Rt_0.975")]
  colnames(state_result_Lancet) <- c("date", "R", "Rt_lower", "Rt_upper")
  state_result_Lancet$date <- as.Date(state_result_Lancet$date)
  state_result_Lancet$R <- as.numeric(as.character(state_result_Lancet$R))
  state_result_Lancet$Rt_lower <- as.numeric(as.character(state_result_Lancet$Rt_lower))
  state_result_Lancet$Rt_upper <- as.numeric(as.character(state_result_Lancet$Rt_upper))
  state_result_Lancet$Rt_method <- "Lancet"
  
  state_result_combined <- cbind(EpiEstim = state_result_EpiEstim$R,
                                 R0 = state_result_R0$R,
                                 Lancet = state_result_Lancet$R)
    
  png(paste0("./COVID19-Rt/sensitivity_analysis/results_3methods_moving_average_state_longitudinal/slope_heatmap_Apr24_3methods_",stateName,".png"), width = 8, height = 8, units = 'in', res = 300)
  print(ggcorplot(data = data.frame(state_result_combined),
                  var_text_size = 5, 
                  cor_text_limits = c(5,10)))
  
  dev.off()
}

