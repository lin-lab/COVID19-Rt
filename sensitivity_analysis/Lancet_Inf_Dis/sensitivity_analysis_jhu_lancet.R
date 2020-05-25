rm(list=ls())
gc()

#Function for sensitivity analysis
#Load in relevant libraries
library (readr)
#Inputs: level (County, State, or USA), pull date from jhu ("%Y-%m-%d"), start date for the data period ("%Y-%m-%d"), end date for the data period ("%Y-%m-%d")
load_jhu <- function( level, pull_date=format(Sys.Date(), "%Y-%m-%d"), start_date, end_date ){
  jhu_url <- paste0("https://raw.githubusercontent.com/lin-lab/COVID-data-cleaning/master/jhu_data/cleaned_data/JHU_COVID-19_",level,"_",pull_date,".csv")
  input_data <- read.csv(url(jhu_url), stringsAsFactors = F)
  jhu_out <- input_data[ as.Date(input_data$date) >= start_date & as.Date(input_data$date) <= end_date, ]
  return(jhu_out)
}

#Load data state-level
jhu_state <- load_jhu(level = "State", pull_date = "2020-05-12", start_date = "2020-03-08", end_date = "2020-05-08")
jhu_state$date <- as.Date(jhu_state$date) 
# parameter (mean=5.2,sd=5.1)
days <- length(unique(jhu_state$date))
mean_gamma <- 5.2
SD_gamma <- 5.1
shape_gamma = (mean_gamma/SD_gamma)^2
rate_gamma = mean_gamma/SD_gamma^2

para_gamma <- cbind(rep(shape_gamma,days),rep(rate_gamma,days))
# save parameter file
write.table(para_gamma,col.names=FALSE,row.names=FALSE,"/n/holystore01/LABS/xlin/Lab/covid19_sensitivity_analysis/Tg_5.2_5.1.txt")

## save Rt value
for (stateName in unique(jhu_state$stateName)){
  state_data <- jhu_state[jhu_state$stateName == stateName,]
  # handle negative increase
  state_data$positiveIncrease <- ifelse(state_data$positiveIncrease < 0, 0, state_data$positiveIncrease)
  # save data
  state_data <- data.frame(date=state_data$date,increase=state_data$positiveIncrease,import=rep(0,dim(state_data)[1]))
  
  write.table(state_data,col.names=FALSE,row.names=FALSE,paste0("/n/holystore01/LABS/xlin/Lab/covid19_sensitivity_analysis/state/state.txt"))

  system(paste0("./a.out 0 100000 0.05 Tg_5.2_5.1.txt /n/holystore01/LABS/xlin/Lab/covid19_sensitivity_analysis/state/state.txt >/n/holystore01/LABS/xlin/Lab/covid19_sensitivity_analysis/state/Rt_",gsub("[[:space:]]", "", stateName),".txt"))
}

## save Rt summary
Rt_US <- c()
for (stateName in unique(jhu_state$stateName)){

	print(stateName)
	Rt <- read.table(paste0("/n/holystore01/LABS/xlin/Lab/covid19_sensitivity_analysis/state/Rt_",gsub("[[:space:]]", "", stateName),".txt"),header=FALSE)
	Rt_mean <- apply(Rt,2,mean)
	Rt_lower <- apply(Rt,2, function(z) quantile(z,0.025))
	Rt_upper <- apply(Rt,2, function(z) quantile(z,0.975))
	
	# dates
	state_data <- jhu_state[jhu_state$stateName == stateName,]
	# Rt summary
	Rt_state <- data.frame(stateName = state_data$stateName, date=state_data$date, 
	Rt=Rt_mean, Rt_0.025=Rt_lower, Rt_0.975=Rt_upper)
	
	Rt_US <- rbind(Rt_US,Rt_state)
}

write.csv(Rt_US,"/n/holystore01/LABS/xlin/Lab/covid19_sensitivity_analysis/Rt_Lancet.csv")


