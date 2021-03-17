---
title: "Diagnostics for Rt, Case Rate, and Death Rate Calculations"
author: "Andy Shi"
date: "2021-03-17"
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



|geo_level   |  rt_upper| case_upper| death_upper| quantile|
|:-----------|---------:|----------:|-----------:|--------:|
|country     |  1.293816|  0.0000167|   0.0000010|   0.5000|
|country     |  1.698601|  0.0001358|   0.0000042|   0.8000|
|country     |  2.259367|  0.0002716|   0.0000071|   0.9000|
|country     |  5.997080|  0.0009877|   0.0000179|   0.9900|
|country     | 15.633333|  0.1520483|   0.0000263|   0.9990|
|country     | 31.927538| 17.3087670|   0.0000336|   0.9999|
|county      |  1.373630|  0.0003198|   0.0000077|   0.5000|
|county      |  1.730567|  0.0007310|   0.0000168|   0.8000|
|county      |  2.127936|  0.0010098|   0.0000246|   0.9000|
|county      |  5.071113|  0.0019928|   0.0000586|   0.9900|
|county      | 15.270051|  0.0040094|   0.0001079|   0.9990|
|county      | 68.009671|  0.0096837|   0.0002853|   0.9999|
|state       |  1.262864|  0.0001718|   0.0000040|   0.5000|
|state       |  1.458186|  0.0004665|   0.0000095|   0.8000|
|state       |  1.622523|  0.0007225|   0.0000129|   0.9000|
|state       |  3.150629|  0.0013326|         Inf|   0.9900|
|state       |  7.190881|  0.0018633|         Inf|   0.9990|
|state       | 15.771895|  0.0020739|         Inf|   0.9999|
|subnational |  1.256279|  0.0000876|   0.0000038|   0.5000|
|subnational |  1.538107|  0.0002502|   0.0000077|   0.8000|
|subnational |  1.771998|  0.0003808|   0.0000105|   0.9000|
|subnational |  3.117860|  0.0007988|   0.0000225|   0.9900|
|subnational |  6.592524|  0.0013496|   0.0000566|   0.9990|
|subnational | 15.752657|  0.0019517|         Inf|   0.9999|

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

|Combined_Key             |  rt_upper|
|:------------------------|---------:|
|Madison, Texas           | 744.29717|
|Medina, Texas            | 314.30089|
|Bahamas                  | 145.72467|
|Bastrop, Texas           | 110.28951|
|Lincoln, Arkansas        | 108.67827|
|Titus, Texas             | 104.05893|
|Woodward, Oklahoma       |  84.16939|
|Columbiana, Ohio         |  83.09680|
|Morgan, Kentucky         |  69.16967|
|Iberville, Louisiana     |  63.44578|
|Burnet, Texas            |  60.71425|
|Sao Tome and Principe    |  57.62610|
|Karnes, Texas            |  56.34749|
|Reeves, Texas            |  53.58266|
|Howard, Texas            |  52.12505|
|Central African Republic |  49.89349|
|Washington, Louisiana    |  49.56265|
|Atascosa, Texas          |  48.40705|
|Marion, Ohio             |  46.44294|
|Ascension, Louisiana     |  43.81949|

 </td>
   <td> 

|Combined_Key            | case_upper|
|:-----------------------|----------:|
|Heilongjiang, China     | 18.4129651|
|Hubei, China            |  2.0402384|
|Marion, Ohio            |  0.0847478|
|Livingston, Missouri    |  0.0489374|
|Lee, Arkansas           |  0.0418891|
|Johnson, Wyoming        |  0.0311574|
|Lyon, Kentucky          |  0.0303053|
|La Salle, Texas         |  0.0196973|
|Bent, Colorado          |  0.0126872|
|Crowley, Colorado       |  0.0108799|
|Dakota, Nebraska        |  0.0101297|
|Finney, Kansas          |  0.0098849|
|Crockett, Texas         |  0.0096837|
|Norton, Kansas          |  0.0089420|
|Madison, Texas          |  0.0080956|
|Lincoln, Arkansas       |  0.0080732|
|Moore, Texas            |  0.0080524|
|Comanche, Kansas        |  0.0080249|
|Pawnee, Kansas          |  0.0077960|
|Bon Homme, South Dakota |  0.0076169|

 </td>
   <td> 

|Combined_Key            | death_upper|
|:-----------------------|-----------:|
|Stanly, North Carolina  |   0.0010251|
|Decatur, Indiana        |   0.0003501|
|Flathead, Montana       |   0.0002853|
|Phelps, Missouri        |   0.0002176|
|Duplin, North Carolina  |   0.0001781|
|Jefferson, Pennsylvania |   0.0001635|
|Penobscot, Maine        |   0.0001531|
|Holmes, Ohio            |   0.0001524|
|Herkimer, New York      |   0.0001496|
|Warren, New Jersey      |   0.0001352|
|Mifflin, Pennsylvania   |   0.0001261|
|Genesee, New York       |   0.0001258|
|Daviess, Kentucky       |   0.0001191|
|Kamchatka Krai, Russia  |   0.0001133|
|Lawrence, Alabama       |   0.0001129|
|Montgomery, New York    |   0.0001119|
|Columbia, New York      |   0.0001113|
|Houston, Alabama        |   0.0001108|
|Covington, Alabama      |   0.0001102|
|Marion, Illinois        |   0.0001093|

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
## 0.000e+00 2.000e+00 4.000e+00 2.378e+15 7.000e+00 5.990e+18
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
##  0.0000  0.0272  0.2234  0.3569  0.6126  1.0000
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
## Overdispersion detected in 0.202492 cases.
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
## [6] purrr_0.3.4   dplyr_1.0.4  
## 
## loaded via a namespace (and not attached):
##  [1] pillar_1.4.7      compiler_4.0.4    highr_0.8         tools_4.0.4      
##  [5] digest_0.6.27     evaluate_0.14     lifecycle_1.0.0   tibble_3.0.6     
##  [9] gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10      DBI_1.1.1        
## [13] cli_2.3.0         rstudioapi_0.13   yaml_2.2.1        xfun_0.21        
## [17] withr_2.4.1       generics_0.1.0    vctrs_0.3.6       hms_1.0.0        
## [21] grid_4.0.4        tidyselect_1.1.0  glue_1.4.2        R6_2.5.0         
## [25] rmarkdown_2.6     farver_2.0.3      magrittr_2.0.1    scales_1.1.1     
## [29] ps_1.5.0          ellipsis_0.3.1    htmltools_0.5.1.1 assertthat_0.2.1 
## [33] colorspace_2.0-0  labeling_0.4.2    stringi_1.5.3     munsell_0.5.0    
## [37] crayon_1.4.1
```
