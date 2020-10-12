#!/bin/bash

Rscript --vanilla calc_initial_estimates.R

# zip files so we can upload to github
zip jhu_county_rt.tsv.zip jhu_county_rt.tsv
zip jhu_global_rt.tsv.zip jhu_global_rt.tsv
zip jhu_state_rt.tsv.zip jhu_state_rt.tsv

Rscript --vanilla upload_to_gdrive.R
