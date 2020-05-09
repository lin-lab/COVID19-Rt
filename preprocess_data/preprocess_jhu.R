#Load in relevant libraries
library (readr)

#Inputs: level (County, State, or USA), pull date from jhu ("%Y-%m-%d"), start date for the data period ("%Y-%m-%d"), end date for the data period ("%Y-%m-%d")
load_jhu <- function( level, pull_date=format(Sys.Date(), "%Y-%m-%d"), start_date, end_date ){
  jhu_url <- paste0("https://raw.githubusercontent.com/lin-lab/COVID-data-cleaning/master/jhu_data/cleaned_data/JHU_COVID-19_",level,"_",pull_date,".csv")
  input_data <- read.csv(url(jhu_url), stringsAsFactors = F)
  jhu_out <- input_data[ as.Date(input_data$date) >= start_date & as.Date(input_data$date) <= end_date, ]
  return(jhu_out)
}
