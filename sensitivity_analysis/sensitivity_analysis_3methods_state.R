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
state_result_Lancet_combined <- read.csv('./COVID19-Rt/sensitivity_analysis/Lancet_Inf_Dis/Rt_Lancet.csv')
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
  
  state_result_combined <- rbind(state_result_EpiEstim, state_result_R0, state_result_Lancet)

  png(paste0("./COVID19-Rt/sensitivity_analysis/results_3methods_state/state_rt/estimate_R_",stateName,".png"), width=8, height=6, units = 'in', res = 300)
  print(ggplot(data = state_result_combined, aes(x=date, y=R, colour=Rt_method)) +
          geom_line() +
          geom_ribbon(aes(ymin=Rt_lower, ymax=Rt_upper, fill = Rt_method), colour = NA, alpha=0.25) + 
          scale_colour_manual(values=c(EpiEstim="red",TD="orange",Lancet="blue")) + 
          scale_fill_manual(values=c(EpiEstim="red",TD="orange",Lancet="blue")) + 
          geom_abline(intercept = 1,slope=0,linetype = 3) + 
          ggtitle("Estimated R"))
  dev.off()
}


######## 3 methods state sensitivity analysis
state_Rt_EpiEstim <- estimate_rt_EpiEstim_state(state_data_combined = jhu_states_pre,
                                                        mean_serial = 5.2,
                                                        std_serial = 5.1)

state_Rt_R0 <- estimate_rt_R0_state(state_data_combined = jhu_states_post,
                                                        mean_serial = 5.2,
                                                        std_serial = 5.1)

state_Rt_Lancet <- state_result_Lancet_combined
non_us_stateName <- unique(state_Rt_Lancet$stateName)[c(3,10,14,15,40,45,54)] # non US states
state_Rt_Lancet <- state_Rt_Lancet[state_Rt_Lancet$stateName %in% setdiff(unique(state_Rt_Lancet$stateName), non_us_stateName),]


### Apr 1
Rt_Apr1_EpiEstim <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim,function(x){x[x$interval_end=="2020-04-01","mean_rt"]}))))
Rt_Apr1_R0 <- as.numeric(as.character(unlist(lapply(state_Rt_R0,function(x){x[x$date=="2020-04-01","rt_est"]}))))
Rt_Apr1_Lancet <- state_Rt_Lancet[state_Rt_Lancet$date == "2020-04-01",c("Rt")]


### Heatmaps
data_for_heatmap_Apr1 <- cbind(EpiEstim = Rt_Apr1_EpiEstim,
                               R0 = Rt_Apr1_R0,
                               Lancet = Rt_Apr1_Lancet)

M1 <- cor(data_for_heatmap_Apr1)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_3methods_state/corr_heatmap_Apr1_3methods.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_3methods_state/slope_heatmap_Apr1_3methods.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_Apr1),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### Apr 24
Rt_Apr24_EpiEstim <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim,function(x){x[x$interval_end=="2020-04-24","mean_rt"]}))))
Rt_Apr24_R0 <- as.numeric(as.character(unlist(lapply(state_Rt_R0,function(x){x[x$date=="2020-04-24","rt_est"]}))))
Rt_Apr24_Lancet <- state_Rt_Lancet[state_Rt_Lancet$date == "2020-04-24",c("Rt")]


data_for_heatmap_Apr24 <- cbind(EpiEstim = Rt_Apr24_EpiEstim,
                                R0 = Rt_Apr24_R0,
                                Lancet = Rt_Apr24_Lancet)

M1 <- cor(data_for_heatmap_Apr24)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_3methods_state/corr_heatmap_Apr24_3methods.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_3methods_state/slope_heatmap_Apr24_3methods.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_Apr24),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### May 3
Rt_May3_EpiEstim <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim,function(x){x[x$interval_end=="2020-05-03","mean_rt"]}))))
Rt_May3_R0 <- as.numeric(as.character(unlist(lapply(state_Rt_R0,function(x){x[x$date=="2020-05-03","rt_est"]}))))
Rt_May3_Lancet <- state_Rt_Lancet[state_Rt_Lancet$date == "2020-05-03",c("Rt")]


data_for_heatmap_May3 <- cbind(EpiEstim = Rt_May3_EpiEstim,
                               R0 = Rt_May3_R0,
                               Lancet = Rt_May3_Lancet)

M1 <- cor(data_for_heatmap_May3)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_3methods_state/corr_heatmap_May3_3methods.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_3methods_state/slope_heatmap_May3_3methods.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_May3),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### Stay at home order
lockdown <- read.csv("~/Dropbox/COVID19_Project/sensitivity_analysis/lockdown_us.csv")
head(lockdown)
unique(lockdown$State)
lockdown <- lockdown[!duplicated(lockdown$State),]

Rt_stayathome_EpiEstim <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim,function(x){x[as.character(x$interval_end)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"mean_rt"]}))))
Rt_stayathome_R0 <- as.numeric(as.character(unlist(lapply(state_Rt_R0,function(x){x[as.character(x$date)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"rt_est"]}))))
lockdown$Date <- as.Date(lockdown$Date)
Rt_stayathome_Lancet <- left_join(lockdown,state_Rt_Lancet,by=c("State"="stateName",
                                                                "Date"="date"))$Rt
Rt_stayathome_Lancet <- Rt_stayathome_Lancet[!is.na(Rt_stayathome_Lancet)]


data_for_heatmap_stayathome <- cbind(EpiEstim = Rt_stayathome_EpiEstim,
                                     R0 = Rt_stayathome_R0,
                                     Lancet = Rt_stayathome_Lancet)

M1 <- cor(data_for_heatmap_stayathome)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_3methods_state/corr_heatmap_stayathome_3methods.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(-1,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_3methods_state/slope_heatmap_stayathome_3methods.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_stayathome),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()

