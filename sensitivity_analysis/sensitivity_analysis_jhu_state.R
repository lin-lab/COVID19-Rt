#Load in relevant libraries
library(readr)
library(EpiEstim)
library(corrplot)
library(RColorBrewer)
source("./COVID19-Rt/sensitivity_analysis/ggcorplot.R")
source('./COVID19-Rt/sensitivity_analysis/sensitivity_analysis.R')
source('./COVID19-Rt/preprocess_data/preprocess_jhu.R')
source('./COVID19-Rt/estimate_rt/estimate_rt_master.R')


######## JHU state
jhu_states <- load_jhu(level='State',pull_date = '2020-05-12', start_date = '2020-03-08', end_date = '2020-05-08')
non_us_stateName <- unique(jhu_states$stateName)[c(3,10,14,15,40,45,54)] # non US states
jhu_states <- jhu_states[jhu_states$stateName %in% setdiff(unique(jhu_states$stateName), non_us_stateName),]

## Rt plot
for (stateName in unique(jhu_states$stateName)){
  state_data <- jhu_states[jhu_states$stateName == stateName,]
  mean_serial_vec <- c(5.2, 7.5, 4.7, 3.96, 4.4); std_serial_vec <- c(5.1, 3.4, 2.9, 4.75, 3.0)
  state_result <- sensitivity_si_EpiEstim(dates = state_data$date,
                                          positive_increase = state_data$positiveIncrease,
                                          mean_serial_vec = mean_serial_vec,
                                          std_serial_vec = std_serial_vec)
  png(paste0("./COVID19-Rt/sensitivity_analysis/results_jhu_state/state_rt/estimate_R_",stateName,".png"), width=8, height=6, units = 'in', res = 300)
  print(estimate_R_plots(state_result, what = "R",
                         options_R = list(col = c("red", "orange", "blue", "green","purple")), legend = TRUE))
  dev.off()
}


######## JHU state sensitivity analysis
state_Rt_EpiEstim_5.2_5.1 <- estimate_rt_EpiEstim_state(state_data_combined = jhu_states,
                                                        mean_serial = 5.2,
                                                        std_serial = 5.1)

state_Rt_EpiEstim_7.5_3.4 <- estimate_rt_EpiEstim_state(state_data_combined = jhu_states,
                                                        mean_serial = 7.5,
                                                        std_serial = 3.4)

state_Rt_EpiEstim_4.7_2.9 <- estimate_rt_EpiEstim_state(state_data_combined = jhu_states,
                                                        mean_serial = 4.7,
                                                        std_serial = 2.9)

state_Rt_EpiEstim_3.96_4.75 <- estimate_rt_EpiEstim_state(state_data_combined = jhu_states,
                                                          mean_serial = 3.96,
                                                          std_serial = 4.75)

state_Rt_EpiEstim_4.4_3.0 <- estimate_rt_EpiEstim_state(state_data_combined = jhu_states,
                                                        mean_serial = 4.4,
                                                        std_serial = 3.0)

### Apr 1
Rt_Apr1_EpiEstim_5.2_5.1 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_5.2_5.1,function(x){x[x$interval_end=="2020-04-01","mean_rt"]}))))
Rt_Apr1_EpiEstim_7.5_3.4 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_7.5_3.4,function(x){x[x$interval_end=="2020-04-01","mean_rt"]}))))
Rt_Apr1_EpiEstim_4.7_2.9 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_4.7_2.9,function(x){x[x$interval_end=="2020-04-01","mean_rt"]}))))
Rt_Apr1_EpiEstim_3.96_4.75 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_3.96_4.75,function(x){x[x$interval_end=="2020-04-01","mean_rt"]}))))
Rt_Apr1_EpiEstim_4.4_3.0 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_4.4_3.0,function(x){x[x$interval_end=="2020-04-01","mean_rt"]}))))


### Heatmaps
data_for_heatmap_Apr1 <- cbind(si_5.2_5.1 = Rt_Apr1_EpiEstim_5.2_5.1,
                               si_7.5_3.4 = Rt_Apr1_EpiEstim_7.5_3.4,
                               si_4.7_2.9 = Rt_Apr1_EpiEstim_4.7_2.9,
                               si_3.96_4.75 = Rt_Apr1_EpiEstim_3.96_4.75,
                               si_4.4_3.0 = Rt_Apr1_EpiEstim_4.4_3.0)

M1 <- cor(data_for_heatmap_Apr1)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_jhu_state/corr_heatmap_Apr1_jhu_si.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_jhu_state/slope_heatmap_Apr1_jhu_si.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_Apr1),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### Apr 24
Rt_Apr24_EpiEstim_5.2_5.1 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_5.2_5.1,function(x){x[x$interval_end=="2020-04-24","mean_rt"]}))))
Rt_Apr24_EpiEstim_7.5_3.4 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_7.5_3.4,function(x){x[x$interval_end=="2020-04-24","mean_rt"]}))))
Rt_Apr24_EpiEstim_4.7_2.9 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_4.7_2.9,function(x){x[x$interval_end=="2020-04-24","mean_rt"]}))))
Rt_Apr24_EpiEstim_3.96_4.75 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_3.96_4.75,function(x){x[x$interval_end=="2020-04-24","mean_rt"]}))))
Rt_Apr24_EpiEstim_4.4_3.0 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_4.4_3.0,function(x){x[x$interval_end=="2020-04-24","mean_rt"]}))))


data_for_heatmap_Apr24 <- cbind(si_5.2_5.1 = Rt_Apr24_EpiEstim_5.2_5.1,
                                si_7.5_3.4 = Rt_Apr24_EpiEstim_7.5_3.4,
                                si_4.7_2.9 = Rt_Apr24_EpiEstim_4.7_2.9,
                                si_3.96_4.75 = Rt_Apr24_EpiEstim_3.96_4.75,
                                si_4.4_3.0 = Rt_Apr24_EpiEstim_4.4_3.0)

M1 <- cor(data_for_heatmap_Apr24)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_jhu_state/corr_heatmap_Apr24_jhu_si.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_jhu_state/slope_heatmap_Apr24_jhu_si.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_Apr24),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### May 3
Rt_May3_EpiEstim_5.2_5.1 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_5.2_5.1,function(x){x[x$interval_end=="2020-05-03","mean_rt"]}))))
Rt_May3_EpiEstim_7.5_3.4 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_7.5_3.4,function(x){x[x$interval_end=="2020-05-03","mean_rt"]}))))
Rt_May3_EpiEstim_4.7_2.9 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_4.7_2.9,function(x){x[x$interval_end=="2020-05-03","mean_rt"]}))))
Rt_May3_EpiEstim_3.96_4.75 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_3.96_4.75,function(x){x[x$interval_end=="2020-05-03","mean_rt"]}))))
Rt_May3_EpiEstim_4.4_3.0 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_4.4_3.0,function(x){x[x$interval_end=="2020-05-03","mean_rt"]}))))


data_for_heatmap_May3 <- cbind(si_5.2_5.1 = Rt_May3_EpiEstim_5.2_5.1,
                               si_7.5_3.4 = Rt_May3_EpiEstim_7.5_3.4,
                               si_4.7_2.9 = Rt_May3_EpiEstim_4.7_2.9,
                               si_3.96_4.75 = Rt_May3_EpiEstim_3.96_4.75,
                               si_4.4_3.0 = Rt_May3_EpiEstim_4.4_3.0)

M1 <- cor(data_for_heatmap_May3)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_jhu_state/corr_heatmap_May3_jhu_si.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_jhu_state/slope_heatmap_May3_jhu_si.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_May3),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### Stay at home order
lockdown <- read.csv("~/Dropbox/COVID19_Project/sensitivity_analysis/lockdown_us.csv")
head(lockdown)
unique(lockdown$State)
lockdown <- lockdown[!duplicated(lockdown$State),]

Rt_stayathome_EpiEstim_5.2_5.1 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_5.2_5.1,function(x){x[as.character(x$interval_end)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"mean_rt"]}))))
Rt_stayathome_EpiEstim_7.5_3.4 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_7.5_3.4,function(x){x[as.character(x$interval_end)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"mean_rt"]}))))
Rt_stayathome_EpiEstim_4.7_2.9 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_4.7_2.9,function(x){x[as.character(x$interval_end)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"mean_rt"]}))))
Rt_stayathome_EpiEstim_3.96_4.75 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_3.96_4.75,function(x){x[as.character(x$interval_end)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"mean_rt"]}))))
Rt_stayathome_EpiEstim_4.4_3.0 <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_4.4_3.0,function(x){x[as.character(x$interval_end)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"mean_rt"]}))))


data_for_heatmap_stayathome <- cbind(si_5.2_5.1 = Rt_stayathome_EpiEstim_5.2_5.1,
                                     si_7.5_3.4 = Rt_stayathome_EpiEstim_7.5_3.4,
                                     si_4.7_2.9 = Rt_stayathome_EpiEstim_4.7_2.9,
                                     si_3.96_4.75 = Rt_stayathome_EpiEstim_3.96_4.75,
                                     si_4.4_3.0 = Rt_stayathome_EpiEstim_4.4_3.0)

M1 <- cor(data_for_heatmap_stayathome)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_jhu_state/corr_heatmap_stayathome_jhu_si.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_jhu_state/slope_heatmap_stayathome_jhu_si.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_stayathome),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()

