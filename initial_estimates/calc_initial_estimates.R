source('./COVID19-Rt/preprocess_data/preprocess_jhu.R')
source('./COVID19-Rt/estimate_rt/estimate_rt_master.R')

#Load in data from data-cleaning github
jhu_counties <- load_jhu(level='County',pull_date = '2020-05-12', start_date = '2020-03-01', end_date = '2020-05-11')
jhu_states <- load_jhu(level='State',pull_date = '2020-05-12', start_date = '2020-03-01', end_date = '2020-05-11')
jhu_global <- load_jhu(level='Global',pull_date = '2020-05-12', start_date = '2020-03-01', end_date = '2020-05-11')

#Perform Rt estimation across all states, indexing by stateName
state_rt_lst <- list()
uniq_states <- unique(jhu_states$stateName)
for (ii in 1:length(uniq_states)) {
  cur_state <- uniq_states[ii]
  jhu_state_sub <- jhu_states[jhu_states$stateName == cur_state, ]
  state_rt_lst[[ii]] <-
    cbind(estimate_rt_EpiEstim(jhu_state_sub$date,
                               jhu_state_sub$positiveIncrease,
                               mean_serial = 5.2, std_serial = 5.1),
          cur_state)
}
state_rt <- do.call(rbind, state_rt_lst)
names(state_rt)[c(2,5)] <- c('date','stateName')
state_rt$date <- as.character(state_rt$date)
jhu_state_out <- merge(jhu_states, state_rt, by=c('stateName','date'))

#Perform Rt estimation across counties, indexing by UID
county_rt_lst <- list()
uniq_uid <- unique(jhu_counties$UID)
for (ii in 1:length(uniq_uid)) {
  cur_uid <- uniq_uid[ii]
  jhu_county_sub <- jhu_counties[jhu_counties$UID == cur_uid,]
  county_rt_lst[[ii]] <-
    cbind(estimate_rt_EpiEstim(jhu_county_sub$date,
                               jhu_county_sub$positiveIncrease,
                               mean_serial = 5.2, std_serial = 5.1),
          cur_uid)
}
county_rt <- do.call(rbind, county_rt_lst)
names(county_rt)[c(2,5)] <- c('date','UID')
county_rt$date <- as.character(county_rt$date)
jhu_county_out <- merge(jhu_counties, county_rt, by=c('UID','date'))

#Perform Rt estimation across countries, indexing by UID
global_rt_lst <- list()
uniq_uid <- unique(jhu_global$UID)
for (ii in 1:length(uniq_uid)) {
  cur_uid <- uniq_uid[ii]
  jhu_global_sub <- jhu_global[jhu_global$UID == cur_uid,]
  global_rt_lst[[ii]] <-
    cbind(estimate_rt_EpiEstim(jhu_global_sub$date,
                               jhu_global_sub$positiveIncrease,
                               mean_serial = 5.2, std_serial = 5.1),
          cur_uid)
}
global_rt <- do.call(rbind, global_rt_lst)
names(global_rt)[c(2,5)] <- c('date','UID')
global_rt$date <- as.character(global_rt$date)
jhu_global_out <- merge(jhu_global, global_rt, by=c('UID','date'))

write.table(jhu_county_out, "./COVID19-Rt/initial_estimates/jhu_county_rt.csv", quote = F, row.names = F, col.names = T, sep=',')
write.table(jhu_state_out, "./COVID19-Rt/initial_estimates/jhu_state_rt.csv", quote = F, row.names = F, col.names = T, sep=',')
write.table(jhu_global_out, "./COVID19-Rt/initial_estimates/jhu_global_rt.csv", quote = F, row.names = F, col.names = T, sep=',')
