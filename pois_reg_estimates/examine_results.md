---
title: "Diagnostics for Rt, Case Rate, and Death Rate Calculations"
author: "Andy Shi"
date: "2021-11-16"
output:
    html_document:
        keep_md: true
---


```r
library(dplyr)
library(purrr)
library(readr)
library(knitr)
library(stringr)
library(ggplot2)
library(cowplot)
```

First, we load libraries and read in the data.


```r
colspec <- cols(
  UID = col_integer(),
  county = col_character(),
  stateName = col_character(),
  date = col_date(format = ""),
  FIPS = col_integer(),
  positiveIncrease = col_integer(),
  deathIncrease = col_integer(),
  positive = col_integer(),
  death = col_integer(),
  population = col_integer(),
  Combined_Key = col_character(),
  rt = col_double(),
  rt_lower = col_double(),
  rt_upper = col_double(),
  case_rate = col_double(),
  case_lower = col_double(),
  case_upper = col_double(),
  death_rate = col_double(),
  death_lower = col_double(),
  death_upper = col_double(),
  date_lag = col_date(format = "")
)

state_fit <- read_csv("results/jhu_state_rt_case_death_rate.csv",
                      col_types = colspec) %>%
  mutate(geo_level = "state")
county_fit <- read_csv("results/jhu_county_rt_case_death_rate.csv",
                       col_types = colspec) %>%
  mutate(geo_level = "county")
global_fit <- read_csv("results/jhu_global_rt_case_death_rate.csv",
                       guess_max = 1e5) %>%
  mutate(geo_level = "country")
subnational_fit <- read_csv("results/jhu_subnational_rt_case_death_rate.csv",
                            guess_max = 1e4) %>%
  mutate(geo_level = "subnational")

fit_all <- bind_rows(state_fit, county_fit, global_fit, subnational_fit)
```

# How Big is the Upper CI

## Quantiles

Here we calculate quantiles of upper CI by resolution (country, state, county,
subnational).


```r
check_probs <- c(0.5, 0.8, 0.9, 0.99, 0.999, 0.9999)
check_quantile <- purrr::partial(quantile, probs = check_probs, na.rm = TRUE)

quantile_df <- fit_all %>%
  group_by(geo_level) %>%
  summarize(across(ends_with("_upper"), check_quantile)) %>%
  mutate(quantile = check_probs)
```

```
## `summarise()` has grouped output by 'geo_level'. You can override using the `.groups` argument.
```

```r
kable(quantile_df, n = Inf)
```



|geo_level   |   rt_upper| case_upper|  death_upper| quantile|
|:-----------|----------:|----------:|------------:|--------:|
|country     |   1.298344|  0.0000275| 1.000000e-06|   0.5000|
|country     |   1.730980|  0.0001915| 4.300000e-06|   0.8000|
|country     |   2.218589|  0.0003666| 7.500000e-06|   0.9000|
|country     |   5.762953|  0.0016983| 2.260000e-05|   0.9900|
|country     |  14.193032|        Inf| 1.235203e+07|   0.9990|
|country     |  19.929478|        Inf| 2.389281e+15|   0.9999|
|county      |   1.548273|  0.0002768| 5.300000e-06|   0.5000|
|county      |   2.014690|  0.0006879| 1.310000e-05|   0.8000|
|county      |   2.441380|  0.0009832| 2.000000e-05|   0.9000|
|county      |   5.709148|  0.0019430| 5.210000e-05|   0.9900|
|county      |  37.403941|  0.0040558| 1.486000e-04|   0.9990|
|county      | 652.437692|  0.0573105|          Inf|   0.9999|
|state       |   1.337868|  0.0001847| 3.200000e-06|   0.5000|
|state       |   1.665815|  0.0004798| 7.900000e-06|   0.8000|
|state       |   1.938941|  0.0006787| 1.200000e-05|   0.9000|
|state       |   3.635303|  0.0012220| 2.440000e-05|   0.9900|
|state       |  11.491222|  0.0019428| 4.620000e-05|   0.9990|
|state       |  42.351092|  0.0019893| 5.910000e-05|   0.9999|
|subnational |   1.245656|  0.0000890| 3.700000e-06|   0.5000|
|subnational |   1.547706|  0.0002587| 8.100000e-06|   0.8000|
|subnational |   1.769602|  0.0003988| 1.180000e-05|   0.9000|
|subnational |   3.041493|  0.0009111| 2.920000e-05|   0.9900|
|subnational |   7.397128|  0.0015763| 5.030000e-05|   0.9990|
|subnational |  16.386878|  0.0041992|          Inf|   0.9999|

## Worst Offenders

Here we show which locations have the most unstable estimates.


```r
top_rates <- function(df, var = c("rt_upper", "case_upper", "death_upper")) {
  var <- match.arg(var)
  ret <- df %>%
    filter(get(var) < Inf) %>%
    filter(!is.na(get(var))) %>%
    group_by(UID) %>%
    slice_max(n = 1, order_by = get(var), with_ties = FALSE) %>%
    ungroup() %>%
    select(Combined_Key, !!var) %>%
    arrange(desc(get(var))) %>%
    head(n = 20)
  return(ret)
}

top_rates_lst <- lapply(c("rt_upper", "case_upper", "death_upper"),
                        function(v) { top_rates(fit_all, v) })
kable(top_rates_lst)
```



<table class="kable_wrapper">
<tbody>
  <tr>
   <td> 

|Combined_Key           |   rt_upper|
|:----------------------|----------:|
|Livingston, Missouri   | 8660.06441|
|Karnes, Texas          | 7243.46641|
|South Dakota           | 3166.81817|
|Rhea, Tennessee        | 1564.82295|
|Jasper, Iowa           | 1325.56839|
|Dodge, Nebraska        |  765.69284|
|Madison, Texas         |  744.29717|
|Letcher, Kentucky      |  652.43769|
|DeWitt, Texas          |  504.37087|
|Laurel, Kentucky       |  274.13571|
|Marion, Iowa           |  179.68160|
|Scotts Bluff, Nebraska |  132.92298|
|Titus, Texas           |  124.47142|
|Grimes, Texas          |  119.67576|
|Walker, Texas          |  112.61004|
|Lincoln, Arkansas      |  108.67827|
|Lancaster, Nebraska    |  107.17901|
|Bastrop, Texas         |  103.13029|
|Moore, Texas           |   93.83180|
|Woodward, Oklahoma     |   86.37044|

 </td>
   <td> 

|Combined_Key                           |    case_upper|
|:--------------------------------------|-------------:|
|Grenada                                | 2.671400e+270|
|Jiangsu, China                         | 2.740080e+240|
|British Virgin Islands, United Kingdom |  3.093824e+05|
|Crowley, Colorado                      |  1.699501e+05|
|Hubei, China                           |  2.090313e+00|
|Marion, Ohio                           |  6.677810e-02|
|Aleutians East, Alaska                 |  5.731050e-02|
|Rock, Nebraska                         |  2.340740e-02|
|Bent, Colorado                         |  2.110620e-02|
|Daviess, Missouri                      |  1.961520e-02|
|Goliad, Texas                          |  1.802780e-02|
|Lincoln, Arkansas                      |  1.800700e-02|
|La Salle, Texas                        |  1.563990e-02|
|Lee, Arkansas                          |  1.530870e-02|
|Livingston, Missouri                   |  1.509000e-02|
|San Juan, Colorado                     |  1.354970e-02|
|Dakota, Nebraska                       |  1.101220e-02|
|Norton, Kansas                         |  9.891400e-03|
|Terrell, Texas                         |  9.702200e-03|
|Emmons, North Dakota                   |  9.129400e-03|

 </td>
   <td> 

|Combined_Key             |  death_upper|
|:------------------------|------------:|
|Elkhart, Indiana         | 3.646890e+60|
|Hawaii, Hawaii           | 5.071304e+25|
|Nova Scotia, Canada      | 4.688131e+15|
|McDowell, North Carolina | 1.508332e+01|
|Okmulgee, Oklahoma       | 7.138200e-03|
|Custer, Oklahoma         | 3.306000e-04|
|Caddo, Oklahoma          | 3.233000e-04|
|Luna, New Mexico         | 2.350000e-04|
|Bee, Texas               | 2.291000e-04|
|Lesotho                  | 1.928000e-04|
|Perry, Kentucky          | 1.807000e-04|
|Lawrence, Alabama        | 1.781000e-04|
|Pitt, North Carolina     | 1.771000e-04|
|Orange, New York         | 1.667000e-04|
|Burke, North Carolina    | 1.621000e-04|
|Preble, Ohio             | 1.564000e-04|
|Hancock, West Virginia   | 1.548000e-04|
|Rockland, New York       | 1.511000e-04|
|Herkimer, New York       | 1.496000e-04|
|Belmont, Ohio            | 1.495000e-04|

 </td>
  </tr>
</tbody>
</table>

# Negative Binomial Diagnostics

Here we check out what the theta values were for the negative binomial
regression. The command line finds the most recent log file and extracts the
theta values.

First, we show the histogram of theta values where the negative binomial
successfully fitted.


```r
# get the most recent log file, grep for thetas, print them out
cmd <- "find -name '*.log.gz' | sort -r | head -n 1 | xargs -I {} zcat {} | grep -Po 'NB fit success: theta = ([0-9.]+)' | cut --delim=' ' -f6"

theta_vals <- as.numeric(system(cmd, intern = TRUE))
summary(theta_vals)
```

```
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## 0.000e+00 2.000e+00 3.000e+00 3.161e+15 7.000e+00 6.588e+18
```

```r
data.frame(theta_vals = theta_vals) %>%
  ggplot(aes(x = theta_vals)) +
  geom_histogram(bins = 50, color = "black") +
  scale_x_log10() +
  xlab(expression(theta)) +
  ggtitle("Histogram of nb theta for successful fits") +
  theme(text = element_text(size = 16))
```

![](examine_results_files/figure-html/negbin_theta_success-1.png)<!-- -->

# Poisson Diagnostics

When the negative binomial fitting fails, we revert to either Poisson or
quasipoisson. We first fit a Poisson, then check for overdispersion. If there's
overdispersion, we fit the quasipoisson.


```r
cmd <- "find -name '*.log.gz' | sort -r | head -n 1 | xargs -I {} zcat {} | grep 'pchisq'"

chisq_cmd <- system(cmd, intern = TRUE)
pchisq_matches <- str_match(chisq_cmd, "pchisq = ([0-9]\\.[0-9\\-e\\+]+)")
pchisq_vals <- as.double(pchisq_matches[, 2])

summary(pchisq_vals)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.0000  0.1668  0.3573  0.7594  1.0000
```

```r
ggplot(data.frame(pchisq_vals = pchisq_vals),
       aes(x = pchisq_vals)) +
  geom_histogram(binwidth = 0.1, color = "black", fill = "grey") +
  theme_cowplot() +
  xlab("Chisq p-values for testing overdispersion") +
  ylab("Count") +
  ggtitle("Histogram of chisq p-values")
```

![](examine_results_files/figure-html/pois_pchisq-1.png)<!-- -->

```r
cat(sprintf("Overdispersion detected in %f cases.\n", mean(pchisq_vals < 0.01)))
```

```
## Overdispersion detected in 0.350000 cases.
```


```r
resid_df_matches <- str_match(chisq_cmd, "([0-9\\.]+) deviance on ([0-9]+) df")
resid_vals <- as.double(resid_df_matches[, 2])
df_vals <- as.double(resid_df_matches[, 3])
resids_df <- data.frame(resids = resid_vals, df = df_vals, pval = pchisq_vals,
                        overdispersed = (pchisq_vals < 0.01))

ggplot(resids_df, aes(x = df, y = resids, color = overdispersed)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, lty = 2, lwd = 1, color = "black") +
  scale_x_log10() + scale_y_log10() +
  ggtitle("Residual Deviance and Residual DF") +
  xlab("Residual degrees of freedom") + ylab("Residual deviance") +
  scale_color_discrete("Identified as\nOverdispersed?") +
  theme_cowplot()
```

![](examine_results_files/figure-html/resids_df-1.png)<!-- -->

# Session Info


```r
sessionInfo()
```

```
## R version 4.0.4 (2021-02-15)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.3 LTS
## 
## Matrix products: default
## BLAS/LAPACK: /opt/OpenBLAS/OpenBLAS_0.3.13/lib/libopenblas-r0.3.13.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] cowplot_1.1.1 ggplot2_3.3.3 stringr_1.4.0 knitr_1.33    readr_1.4.0  
## [6] purrr_0.3.4   dplyr_1.0.6  
## 
## loaded via a namespace (and not attached):
##  [1] highr_0.9         pillar_1.6.1      bslib_0.2.5.1     compiler_4.0.4   
##  [5] jquerylib_0.1.4   tools_4.0.4       digest_0.6.27     jsonlite_1.7.2   
##  [9] evaluate_0.14     lifecycle_1.0.0   tibble_3.1.2      gtable_0.3.0     
## [13] pkgconfig_2.0.3   rlang_0.4.11      rstudioapi_0.13   cli_2.5.0        
## [17] DBI_1.1.1         yaml_2.2.1        xfun_0.23         withr_2.4.2      
## [21] generics_0.1.0    vctrs_0.3.8       sass_0.4.0        hms_1.1.0        
## [25] grid_4.0.4        tidyselect_1.1.1  glue_1.4.2        R6_2.5.0         
## [29] fansi_0.5.0       rmarkdown_2.8     farver_2.1.0      magrittr_2.0.1   
## [33] ps_1.6.0          scales_1.1.1      ellipsis_0.3.2    htmltools_0.5.1.1
## [37] assertthat_0.2.1  colorspace_2.0-1  labeling_0.4.2    utf8_1.2.1       
## [41] stringi_1.6.2     munsell_0.5.0     crayon_1.4.1
```
