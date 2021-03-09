#!/bin/bash

set -euxo pipefail

# set path to Rscript
RSCRIPT=/opt/R/4.0.4/bin/Rscript
if [[ ! -f ${RSCRIPT} ]]
then
    RSCRIPT=/usr/bin/Rscript
fi

$RSCRIPT -e "rmarkdown::render('examine_results.Rmd')"

date_str=$(date +"%Y-%m-%d_%H%M%S")
html_fname=examine_results.html
out_fname=examine_results_${date_str}.html
out_dir=$HOME/Dropbox/phd_research/covid19/covid-rt-diagnostics

cp $html_fname ${out_dir}/${out_fname}
