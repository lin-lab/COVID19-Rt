#Load in relevant libraries
library (readr)

#Inputs: level (County or State), pull date from Covid tracking ("md" - note not daily), start date for the data period ("%Y-%m-%d"), end date for the data period ("%Y-%m-%d")
load_covid_tracking <- function( level, pull_date="0505", start_date, end_date ){
  covid_tracking_url <- paste0("https://raw.githubusercontent.com/lin-lab/covid-data-cleaning/master/covid-tracking_data/covidtracking_cleaned/covidtracking_",level,"_",pull_date,".csv")
  input_data <- read.csv(url(covid_tracking_url), stringsAsFactors = F)
  covid_tracking_out <- input_data[ as.Date(input_data$date) >= start_date & as.Date(input_data$date) <= end_date, ]
  return(covid_tracking_out)
}
