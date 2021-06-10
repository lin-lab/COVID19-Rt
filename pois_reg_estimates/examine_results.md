---
title: "Diagnostics for Rt, Case Rate, and Death Rate Calculations"
author: "Andy Shi"
date: "2021-06-10"
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
|country     |   1.285111|  0.0000208| 1.100000e-06|   0.5000|
|country     |   1.682777|  0.0001616| 4.400000e-06|   0.8000|
|country     |   2.190031|  0.0003167| 7.700000e-06|   0.9000|
|country     |   5.799516|  0.0010317| 2.130000e-05|   0.9900|
|country     |  12.921135|  0.0232764| 1.442181e+05|   0.9990|
|country     | 824.772653| 17.3087670|          Inf|   0.9999|
|county      |   1.364336|  0.0002592| 6.700000e-06|   0.5000|
|county      |   1.716119|  0.0006338| 1.530000e-05|   0.8000|
|county      |   2.083550|  0.0009224| 2.280000e-05|   0.9000|
|county      |   5.507032|  0.0018769| 5.800000e-05|   0.9900|
|county      |  56.347493|  0.0039405| 1.351000e-04|   0.9990|
|county      | 562.536081|  0.0098686| 7.007607e+29|   0.9999|
|state       |   1.257915|  0.0001581| 3.400000e-06|   0.5000|
|state       |   1.479482|  0.0004002| 8.300000e-06|   0.8000|
|state       |   1.648093|  0.0006511| 1.230000e-05|   0.9000|
|state       |   3.222350|  0.0011914| 4.180000e-05|   0.9900|
|state       |   7.068664|  0.0017126|          Inf|   0.9990|
|state       |  15.000427|  0.0019750|          Inf|   0.9999|
|subnational |   1.244412|  0.0000909| 3.900000e-06|   0.5000|
|subnational |   1.509709|  0.0002776| 8.400000e-06|   0.8000|
|subnational |   1.703627|  0.0004170| 1.240000e-05|   0.9000|
|subnational |   3.024352|  0.0008469| 2.970000e-05|   0.9900|
|subnational |   7.202996|  0.0014500| 4.850000e-05|   0.9990|
|subnational |  32.702677|  0.0049077|          Inf|   0.9999|

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
|Lesotho                | 1136.89200|
|Madison, Texas         |  744.29717|
|Medina, Texas          |  314.30089|
|Rapides, Louisiana     |  147.95422|
|Orange, Texas          |  140.73089|
|Lincoln, Arkansas      |  108.67827|
|South Sudan            |  107.13500|
|Bastrop, Texas         |  106.28627|
|Titus, Texas           |  104.05893|
|Carson City, Nevada    |   99.19458|
|Moore, Texas           |   93.83180|
|Woodward, Oklahoma     |   86.37044|
|Columbiana, Ohio       |   83.09680|
|Madagascar             |   70.29292|
|Morgan, Kentucky       |   69.16967|
|Scotts Bluff, Nebraska |   65.47305|
|Chicot, Arkansas       |   65.44827|
|Bexar, Texas           |   63.92303|
|Iberville, Louisiana   |   63.44578|

 </td>
   <td> 

|Combined_Key            | case_upper|
|:-----------------------|----------:|
|Heilongjiang, China     | 18.4129651|
|Hubei, China            |  2.0402384|
|Marion, Ohio            |  0.0802726|
|Johnson, Wyoming        |  0.0476055|
|Lee, Arkansas           |  0.0413296|
|Bent, Colorado          |  0.0412814|
|Daviess, Missouri       |  0.0234203|
|La Salle, Texas         |  0.0183337|
|Livingston, Missouri    |  0.0164446|
|Putnam, Missouri        |  0.0154255|
|Lincoln, Colorado       |  0.0130275|
|Maverick, Texas         |  0.0116431|
|Crowley, Colorado       |  0.0105167|
|Nunavut, Canada         |  0.0105015|
|Dakota, Nebraska        |  0.0101876|
|Knox, Texas             |  0.0098752|
|Finney, Kansas          |  0.0092955|
|Norton, Kansas          |  0.0089420|
|Bon Homme, South Dakota |  0.0086269|
|Scurry, Texas           |  0.0083727|

 </td>
   <td> 

|Combined_Key           |  death_upper|
|:----------------------|------------:|
|Nevada, California     | 2.931638e+34|
|Nova Scotia, Canada    | 6.723059e+06|
|Stanly, North Carolina | 2.872210e-02|
|Luna, New Mexico       | 2.350000e-04|
|Phelps, Missouri       | 2.176000e-04|
|Lea, New Mexico        | 1.850000e-04|
|Duplin, North Carolina | 1.781000e-04|
|Muskogee, Oklahoma     | 1.721000e-04|
|Montgomery, New York   | 1.562000e-04|
|Penobscot, Maine       | 1.531000e-04|
|Herkimer, New York     | 1.496000e-04|
|Cayuga, New York       | 1.360000e-04|
|Warren, New Jersey     | 1.352000e-04|
|Barren, Kentucky       | 1.351000e-04|
|Mifflin, Pennsylvania  | 1.261000e-04|
|Fayette, Pennsylvania  | 1.139000e-04|
|Lawrence, Alabama      | 1.129000e-04|
|Columbia, New York     | 1.113000e-04|
|Jefferson, Ohio        | 1.109000e-04|
|Houston, Alabama       | 1.108000e-04|

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
## 0.000e+00 2.000e+00 4.000e+00 2.473e+15 7.000e+00 6.889e+18
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
## 0.00000 0.01919 0.22323 0.36969 0.68786 1.00000
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
## Overdispersion detected in 0.225519 cases.
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
## Running under: Ubuntu 20.04.2 LTS
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
