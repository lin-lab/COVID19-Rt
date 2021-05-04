---
title: "Diagnostics for Rt, Case Rate, and Death Rate Calculations"
author: "Andy Shi"
date: "2021-05-04"
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



|geo_level   |    rt_upper| case_upper|  death_upper| quantile|
|:-----------|-----------:|----------:|------------:|--------:|
|country     |    1.288587|  0.0000197| 1.000000e-06|   0.5000|
|country     |    1.693487|  0.0001578| 4.300000e-06|   0.8000|
|country     |    2.262666|  0.0003152| 7.400000e-06|   0.9000|
|country     |    6.018337|  0.0010027| 1.930000e-05|   0.9900|
|country     |   14.291991|  0.0304641| 2.920000e-05|   0.9990|
|country     | 1136.892001| 17.3087670| 4.860000e-05|   0.9999|
|county      |    1.368138|  0.0002865| 7.100000e-06|   0.5000|
|county      |    1.720311|  0.0006737| 1.590000e-05|   0.8000|
|county      |    2.090650|  0.0009592| 2.340000e-05|   0.9000|
|county      |    5.011461|  0.0019265| 5.820000e-05|   0.9900|
|county      |   17.667874|  0.0041178| 1.285000e-04|   0.9990|
|county      |  124.489828|  0.0096900| 2.931638e+34|   0.9999|
|state       |    1.266551|  0.0001731| 3.700000e-06|   0.5000|
|state       |    1.476043|  0.0004267| 8.800000e-06|   0.8000|
|state       |    1.651092|  0.0006742| 1.270000e-05|   0.9000|
|state       |    3.114665|  0.0012146| 4.790000e-05|   0.9900|
|state       |    7.016983|  0.0017493|          Inf|   0.9990|
|state       |   16.814506|  0.0019737|          Inf|   0.9999|
|subnational |    1.258501|  0.0000900| 3.800000e-06|   0.5000|
|subnational |    1.534421|  0.0002757| 7.800000e-06|   0.8000|
|subnational |    1.745034|  0.0004166| 1.060000e-05|   0.9000|
|subnational |    2.980330|  0.0008446| 2.310000e-05|   0.9900|
|subnational |    6.539594|  0.0014033| 4.500000e-05|   0.9990|
|subnational |   14.185460|  0.0028883|          Inf|   0.9999|

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

|Combined_Key                       |   rt_upper|
|:----------------------------------|----------:|
|Lesotho                            | 1136.89200|
|Madison, Texas                     |  744.29717|
|Medina, Texas                      |  314.30089|
|Rapides, Louisiana                 |  167.92356|
|Orange, Texas                      |  140.73089|
|Lincoln, Arkansas                  |  108.67827|
|South Sudan                        |  107.13500|
|Bastrop, Texas                     |  106.28627|
|Titus, Texas                       |  104.05893|
|Moore, Texas                       |  103.19862|
|Carson City, Nevada                |   99.19458|
|Woodward, Oklahoma                 |   86.37044|
|Columbiana, Ohio                   |   83.09680|
|San Andres y Providencia, Colombia |   69.22465|
|Morgan, Kentucky                   |   69.16967|
|Scotts Bluff, Nebraska             |   65.47305|
|Iberville, Louisiana               |   63.44578|
|Burnet, Texas                      |   60.71425|
|Livingston, Louisiana              |   59.38281|
|Churchill, Nevada                  |   58.25152|

 </td>
   <td> 

|Combined_Key            | case_upper|
|:-----------------------|----------:|
|Livingston, Missouri    | 42.7531026|
|Heilongjiang, China     | 18.4129651|
|Hubei, China            |  2.0402384|
|Marion, Ohio            |  0.0821275|
|Lee, Arkansas           |  0.0418891|
|Bent, Colorado          |  0.0412709|
|Johnson, Wyoming        |  0.0311574|
|Spink, South Dakota     |  0.0232611|
|La Salle, Texas         |  0.0196529|
|Lincoln, Colorado       |  0.0136576|
|French Guiana, France   |  0.0119647|
|Crowley, Colorado       |  0.0105218|
|Dakota, Nebraska        |  0.0102191|
|Knox, Texas             |  0.0098752|
|Finney, Kansas          |  0.0098556|
|Crockett, Texas         |  0.0096837|
|Norton, Kansas          |  0.0089420|
|Bon Homme, South Dakota |  0.0086269|
|Scurry, Texas           |  0.0080993|
|Lincoln, Arkansas       |  0.0080449|

 </td>
   <td> 

|Combined_Key            |  death_upper|
|:-----------------------|------------:|
|Nevada, California      | 2.931638e+34|
|Stanly, North Carolina  | 2.872210e-02|
|Phelps, Missouri        | 2.176000e-04|
|Duplin, North Carolina  | 1.781000e-04|
|Jefferson, Pennsylvania | 1.635000e-04|
|Montgomery, New York    | 1.562000e-04|
|Penobscot, Maine        | 1.531000e-04|
|Herkimer, New York      | 1.496000e-04|
|Cayuga, New York        | 1.360000e-04|
|Warren, New Jersey      | 1.352000e-04|
|Mifflin, Pennsylvania   | 1.261000e-04|
|Genesee, New York       | 1.258000e-04|
|Lawrence, Alabama       | 1.129000e-04|
|Columbia, New York      | 1.113000e-04|
|Jefferson, Ohio         | 1.109000e-04|
|Houston, Alabama        | 1.108000e-04|
|Covington, Alabama      | 1.102000e-04|
|Lycoming, Pennsylvania  | 1.087000e-04|
|Hill, Texas             | 1.085000e-04|
|Somerset, Pennsylvania  | 1.076000e-04|

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
## 0.000e+00 2.000e+00 4.000e+00 4.118e+15 7.000e+00 1.296e+19
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
## 0.00000 0.02416 0.23252 0.37534 0.70419 1.00000
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
## Overdispersion detected in 0.209091 cases.
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
## [1] cowplot_1.1.1 ggplot2_3.3.3 stringr_1.4.0 knitr_1.31    readr_1.4.0  
## [6] purrr_0.3.4   dplyr_1.0.5  
## 
## loaded via a namespace (and not attached):
##  [1] highr_0.8         pillar_1.5.1      bslib_0.2.4       compiler_4.0.4   
##  [5] jquerylib_0.1.3   tools_4.0.4       digest_0.6.27     jsonlite_1.7.2   
##  [9] evaluate_0.14     lifecycle_1.0.0   tibble_3.1.0      gtable_0.3.0     
## [13] pkgconfig_2.0.3   rlang_0.4.10      rstudioapi_0.13   cli_2.3.1        
## [17] DBI_1.1.1         yaml_2.2.1        xfun_0.22         withr_2.4.1      
## [21] generics_0.1.0    vctrs_0.3.6       sass_0.3.1        hms_1.0.0        
## [25] grid_4.0.4        tidyselect_1.1.0  glue_1.4.2        R6_2.5.0         
## [29] fansi_0.4.2       rmarkdown_2.7     farver_2.1.0      magrittr_2.0.1   
## [33] ps_1.6.0          scales_1.1.1      ellipsis_0.3.1    htmltools_0.5.1.1
## [37] assertthat_0.2.1  colorspace_2.0-0  labeling_0.4.2    utf8_1.2.1       
## [41] stringi_1.5.3     munsell_0.5.0     crayon_1.4.1
```
