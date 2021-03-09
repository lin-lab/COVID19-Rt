#goal: fit international cases regression
#load in the data
library(data.table)
library(splines)
library(geepack)
library(dplyr)
source("../../helpers/helper_functions.R")
source("gee_fit_helpers.R")

setwd("/Users/zhangh24/GoogleDrive/Covid/")
data <- as.data.frame(fread("./data/regression_2020_06_28.csv"))
#constrained the analysis between Marth 15th to April 30th
data$date <- as.Date(data$date, format = "%m/%d/%y")
data$weekday = as.factor(data$weekday)
levels(data$weekday) = c("Friday","Monday","Tuesday","Wednesday",
                         "Thursday","Saturday","Sunday")
data$continent = factor(data$continent,levels = c("Asia","Africa","Europe","North America",
                                                  "Oceania","South America","Central America")
)


idx <- which(data$positive>=50)

#idx <- which(is.na(data$lockdown))
data.clean = data[idx,]

#data.clean$date[1] is 2020-03-16
# calculate days since an event with specified start_date
#every analysis started from March-16th
#data.clean$days_int = days_since("01/01/2020",data.clean$date)
#data.clean$lock_down_days = days_since(data.clean$lockdown_date,data.clean$date,nz=T)

#data.clean = data.clean %>%
#  mutate(lock_down_binary= ifelse(lock_down_days>0,1,0))

#gee model





data.clean <- data.clean %>% filter(positiveIncrease>=0)
#daily death
data.complete = data.clean[complete.cases(data.clean[
  c("positiveIncrease",
    "countryName",
    "population",
    "positive",
    "weekday",
    "continent",
    "population_density",
    "median_age",
    "pct_diabetes",
    "cvd_death_rate",
    "Physicians_per_1000",
    "hospital_beds_per_100k",
    "Health_exp_pct_GDP_2016",
    "gdp_per_capita_PPP",
    "containmentIndex",
    "healthTestingIndex",
    "economicSupportIndex",
    "mobility_avg7",
    "days_int")]),]
countryremove <- names(which(table(data.complete$countryName)<=15))
data.complete = data.complete %>%
  filter(countryName%in%countryremove==F)

# Prepare knots values for B-spline
day_knots <- seq_days(data.complete$days_s_50p, 30)



#daily death regression~offset(log(population))
cases_eqn =  positiveIncrease~
  offset(log(population))+
  log(positive)+
  weekday+
  factor(continent)+
  scale(population_density)+
  scale(median_age)+
  scale(pct_diabetes)+
  scale(cvd_death_rate)+
  scale(Physicians_per_1000)+
  scale(hospital_beds_per_100k)+
  scale(Health_exp_pct_GDP_2016)+
  scale(gdp_per_capita_PPP)+
  scale(mobility_avg7)+
  scale(containmentIndex) +
  scale(healthTestingIndex) +
  scale(economicSupportIndex) +
  bs(days_s_50p, knots = day_knots)


fit = gee_fit(dt = data.complete, eqn = cases_eqn,
              gee_id_var = "countryName", corstr = "ar1")
coefficients(summary(fit))



plot_dt <- reg_coef_data(fit)
name_plot_dt = plot_dt$name
plot_dt$name = factor(name_plot_dt,levels=name_plot_dt)

library(ggplot2)
# p_daily_positiveIncrease <- ggplot(plot_dt) +
#   theme_Publication()+
#   geom_bar(aes(x=name, y=value), stat="identity", fill="royalblue", alpha=0.7) +
#   geom_errorbar(aes(x=name, ymin=value-(1.96*se), ymax=value+(1.96*se), width=0.2),
#                 colour="firebrick2", alpha=0.9, size=0.8) +
#   coord_flip() +
#   # ylim(-5, 15) +
#   labs(x = 'Variables', y = expression("Coefficient" %+-% "1.96 SE"),
#        title = 'Longitudinal model',
#        subtitle = "Daily cases {~offset(log_population)}") +
#   theme(plot.title=element_text(hjust=0.5, size=30),
#         plot.subtitle=element_text(size=20),
#         axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0), hjust=0.5, size=15),
#         axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0), vjust=0.5, size=15),
#         axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15),
#         legend.title=element_text(size=15),
#         legend.text=element_text(size=15))
# png(file = "./result/daily_cases_regression_offsetpopulation_ar1.png",width=8,height=6,units = "in",res=300)
# p_daily_positiveIncrease
# dev.off()

daily_cases_result <- data.frame(Variable = plot_dt$name,
                                 Estimate = round(plot_dt$value,2),
                                 CI = plot_dt$CI,
                                 P = plot_dt$P)
colnames(daily_cases_result)[3] = "95%CI"
write.csv(daily_cases_result,file = "./result/daily_cases_result.csv",row.names = F)

