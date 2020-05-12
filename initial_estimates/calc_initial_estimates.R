source('./COVID19-Rt/preprocess_data/preprocess_jhu.R')
source('./COVID19-Rt/estimate_rt/estimate_rt_master.R')

jhu_counties <- load_jhu(level='County',pull_date = '2020-05-11', start_date = '2020-03-01', end_date = '2020-05-11')
jhu_states <- load_jhu(level='State',pull_date = '2020-05-11', start_date = '2020-03-01', end_date = '2020-05-11')

state_rt <- c()
for (ii in 1:length(unique(jhu_states$stateName))){
  jhu_state_sub <- jhu_states[jhu_states$stateName==unique(jhu_states$stateName)[ii],]
  state_rt <- rbind(state_rt,
  cbind(estimate_rt_EpiEstim(jhu_state_sub$date, jhu_state_sub$positiveIncrease, mean_serial = 5.2, std_serial = 5.1), unique(jhu_states$stateName)[ii]))
}
names(state_rt)[c(2,5)] <- c('date','stateName')
state_rt$date <- as.character(state_rt$date)
jhu_state_out <- merge(jhu_states, state_rt, by=c('stateName','date'))


county_rt <- c()
for (ii in 1:length(unique(jhu_counties$FIPS))){
  jhu_county_sub <- jhu_counties[jhu_counties$FIPS==unique(jhu_counties$FIPS)[ii],]
  county_rt <- rbind(county_rt,
                    cbind(estimate_rt_EpiEstim(jhu_county_sub$date, jhu_county_sub$positiveIncrease, mean_serial = 5.2, std_serial = 5.1), unique(jhu_counties$FIPS)[ii]))
}
names(county_rt)[c(2,5)] <- c('date','FIPS')
county_rt$date <- as.character(county_rt$date)
jhu_county_out <- merge(jhu_counties, county_rt, by=c('FIPS','date'))

write.table(jhu_county_out, "./COVID19-Rt/initial_estimates/jhu_county_rt.csv", quote = F, row.names = F, col.names = T, sep=',')
write.table(jhu_state_out, "./COVID19-Rt/initial_estimates/jhu_state_rt.csv", quote = F, row.names = F, col.names = T, sep=',')
