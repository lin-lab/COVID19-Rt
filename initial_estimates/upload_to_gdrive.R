library(googledrive)

# upload current data
for (ext in c("county", "global", "state")) {
  local_path <- sprintf("jhu_%s_rt.tsv", ext)
  cat(sprintf("Uploading file %s...\n", local_path))
  upload_res <- drive_upload(local_path,
                             path = "COVID19 Rt/",
                             overwrite = TRUE)
  print(upload_res)
}

# upload current time
cur_time <- Sys.time()
last_updated <- sprintf("Last updated: %s", format(cur_time, usetz = TRUE))
write(last_updated, "last_updated.txt")
upload_time <- drive_upload("last_updated.txt",
                            path = "COVID19 Rt/",
                            overwrite = TRUE)
upload_time
