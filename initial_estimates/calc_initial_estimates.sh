#!/bin/bash

set -euxo pipefail

# set path to Rscript
RSCRIPT=/opt/R/4.0.4/bin/Rscript
if [[ ! -f ${RSCRIPT} ]]
then
    RSCRIPT=/usr/bin/Rscript
fi

base_url=https://hsph-covid-study.s3.us-east-2.amazonaws.com/JHU_Cleaned
out_dir=raw_data

if [ ! -d ${out_dir} ]
then
  mkdir -v ${out_dir}
fi

# download zipped files and unzip
for ext in County State Global
do
    filename=JHU_COVID-19_${ext}.csv
    out_zip=${out_dir}/${filename}.zip
    url=${base_url}/${filename}.zip
    wget --output-document=$out_zip $url && unzip -o $out_zip -d ${out_dir} && rm -f $out_zip
done

# call script to calculate estimates
$RSCRIPT --vanilla calc_initial_estimates.R

# zip results and upload to amazon s3
for ext in county state global
do
  tsv=jhu_${ext}_rt.tsv
  out_zip=${tsv}.zip
  zip - ${tsv} > ${out_zip}
  aws s3 cp ${out_zip} s3://hsph-covid-study/Rt-values/
done
