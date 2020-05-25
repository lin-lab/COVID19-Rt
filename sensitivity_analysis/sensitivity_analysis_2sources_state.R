#Load in relevant libraries
library(readr)
library(EpiEstim)
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
jhu_states <- load_jhu(level='State',pull_date = '2020-05-12', start_date = '2020-03-08', end_date = '2020-05-04')
non_us_stateName <- unique(jhu_states$stateName)[c(3,10,14,15,40,45,54)] # non US states
jhu_states <- jhu_states[jhu_states$stateName %in% setdiff(unique(jhu_states$stateName), non_us_stateName),]

######## nyt state
## placeholder

######## yu state
yu_states <- load_yu(level='State',pull_date = '05_05_2020', start_date = '2020-03-08', end_date = '2020-05-04')
yu_states[yu_states$stateName == "District Of Columbia","stateName"] <- "District of Columbia"

## Rt plot
for (stateName in unique(jhu_states$stateName)){
  state_data_jhu <- jhu_states[jhu_states$stateName == stateName,]
  ## placeholder for nyt
  state_data_yu <- yu_states[yu_states$stateName == stateName,]
  
  dates_list <- list(state_data_jhu$date, state_data_yu$date)
  positive_increase_list <- list(state_data_jhu$positiveIncrease, state_data_yu$positiveIncrease)
  state_result <- sensitivity_data_EpiEstim(dates_list = dates_list,
                                            positive_increase_list = positive_increase_list
                                            )
  png(paste0("./COVID19-Rt/sensitivity_analysis/results_2sources_state/state_rt/estimate_R_",stateName,".png"), width=8, height=6, units = 'in', res = 300)
  print(estimate_R_plots(state_result, what = "R",
                         options_R = list(col = c("red", "blue")), legend = TRUE))
  dev.off()
}


######## 2 sources state sensitivity analysis
state_Rt_EpiEstim_jhu <- estimate_rt_EpiEstim_state(state_data_combined = jhu_states,
                                                        mean_serial = 5.2,
                                                        std_serial = 5.1)

state_Rt_EpiEstim_yu <- estimate_rt_EpiEstim_state(state_data_combined = yu_states,
                                                        mean_serial = 5.2,
                                                        std_serial = 5.1)

### Apr 1
Rt_Apr1_EpiEstim_jhu <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_jhu,function(x){x[x$interval_end=="2020-04-01","mean_rt"]}))))
Rt_Apr1_EpiEstim_yu <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_yu,function(x){x[x$interval_end=="2020-04-01","mean_rt"]}))))


### Heatmaps
data_for_heatmap_Apr1 <- cbind(jhu = Rt_Apr1_EpiEstim_jhu,
                               yu = Rt_Apr1_EpiEstim_yu)

M1 <- cor(data_for_heatmap_Apr1)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_2sources_state/corr_heatmap_Apr1_2sources.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_2sources_state/slope_heatmap_Apr1_2sources.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_Apr1),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### Apr 24
Rt_Apr24_EpiEstim_jhu <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_jhu,function(x){x[x$interval_end=="2020-04-24","mean_rt"]}))))
Rt_Apr24_EpiEstim_yu <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_yu,function(x){x[x$interval_end=="2020-04-24","mean_rt"]}))))


data_for_heatmap_Apr24 <- cbind(jhu = Rt_Apr24_EpiEstim_jhu,
                               yu = Rt_Apr24_EpiEstim_yu)

M1 <- cor(data_for_heatmap_Apr24)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_2sources_state/corr_heatmap_Apr24_2sources.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_2sources_state/slope_heatmap_Apr24_2sources.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_Apr24),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### May 3
Rt_May3_EpiEstim_jhu <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_jhu,function(x){x[x$interval_end=="2020-05-03","mean_rt"]}))))
Rt_May3_EpiEstim_yu <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_yu,function(x){x[x$interval_end=="2020-05-03","mean_rt"]}))))


data_for_heatmap_May3 <- cbind(jhu = Rt_May3_EpiEstim_jhu,
                                yu = Rt_May3_EpiEstim_yu)

M1 <- cor(data_for_heatmap_May3)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_2sources_state/corr_heatmap_May3_2sources.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_2sources_state/slope_heatmap_May3_2sources.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_May3),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()


### Stay at home order
lockdown <- read.csv("~/Dropbox/COVID19_Project/sensitivity_analysis/lockdown_us.csv")
head(lockdown)
unique(lockdown$State)
lockdown <- lockdown[!duplicated(lockdown$State),]

Rt_stayathome_EpiEstim_jhu <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_jhu,function(x){x[as.character(x$interval_end)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"mean_rt"]}))))
Rt_stayathome_EpiEstim_yu <- as.numeric(as.character(unlist(lapply(state_Rt_EpiEstim_yu,function(x){x[as.character(x$interval_end)==as.character(lockdown[lockdown$State == x$state[1],"Date"]),"mean_rt"]}))))


data_for_heatmap_stayathome <- cbind(jhu = Rt_stayathome_EpiEstim_jhu,
                               yu = Rt_stayathome_EpiEstim_yu)

M1 <- cor(data_for_heatmap_stayathome)
M1[which(M1 > 1)] <- 1
png("./COVID19-Rt/sensitivity_analysis/results_2sources_state/corr_heatmap_stayathome_2sources.png", width = 8, height = 8, units = 'in', res = 300)
corrplot(M1,
         method = "circle",
         tl.col = "black",
         type = "upper",
         cl.lim = c(0,1),
         addCoef.col = "white", number.digits = 2, number.cex = 1.25,
         col = colorRampPalette(c("blue", "white", "red"))(100), cl.cex=1)
dev.off()

png("./COVID19-Rt/sensitivity_analysis/results_2sources_state/slope_heatmap_stayathome_2sources.png", width = 8, height = 8, units = 'in', res = 300)
ggcorplot(data = data.frame(data_for_heatmap_stayathome),
          var_text_size = 5, 
          cor_text_limits = c(5,10))

dev.off()

