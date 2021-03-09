#!/bin/bash

set -euxo pipefail

# set path to Rscript
RSCRIPT=/opt/R/4.0.4/bin/Rscript
if [[ ! -f ${RSCRIPT} ]]
then
    RSCRIPT=/usr/bin/Rscript
fi

# call script to calculate estimates
${RSCRIPT} --vanilla calc_pois_estimates.R
bash examine_results.sh

# zip results and upload to amazon s3
cd results
for ext in county state global subnational
do
  csv=jhu_${ext}_rt_case_death_rate.csv
  out_zip=${csv}.zip
  zip - ${csv} > ${out_zip}
  aws s3 cp ${out_zip} s3://hsph-covid-study/pois_metrics_values/
done
