#Inputs: level (County, State, or USA), pull date from jhu ("%Y-%m-%d"), start date for the data period ("%Y-%m-%d"), end date for the data period ("%Y-%m-%d")
load_jhu <- function(level, start_date, end_date){
  #jhu_url <- paste0("https://raw.githubusercontent.com/lin-lab/COVID-data-cleaning/master/jhu_data/cleaned_data/JHU_COVID-19_", level, ".csv")
  #jhu_path <- url(jhu_url)
  jhu_path <- paste0("raw_data/JHU_COVID-19_", level, ".csv")
  input_data <- read.csv(jhu_path, stringsAsFactors = FALSE)
  jhu_out <- input_data[ as.Date(input_data$date) >= as.Date(start_date) & as.Date(input_data$date) <= as.Date(end_date), ]
  return(jhu_out)
}
