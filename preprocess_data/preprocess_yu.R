#Load in relevant libraries
library (readr)

#Inputs: level (County or State), pull date from Yu ("%m_%d_%y" - note not daily), start date for the data period ("%Y-%m-%d"), end date for the data period ("%Y-%m-%d")
load_yu <- function( level, pull_date="05_05_2020", start_date, end_date ){
  yu_url <- paste0("https://raw.githubusercontent.com/lin-lab/COVID-data-cleaning/master/yu_group/cleaned_data/yu_Group_",level,".",pull_date,".csv")
  input_data <- read.csv(url(yu_url), stringsAsFactors = F)
  yu_out <- input_data[ as.Date(input_data$date) >= start_date & as.Date(input_data$date) <= end_date, ]
  return(yu_out)
}
