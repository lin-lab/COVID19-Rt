---
title: "Diagnostics for Rt, Case Rate, and Death Rate Calculations"
author: "Andy Shi"
date: "2021-03-30"
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



|geo_level   |    rt_upper| case_upper| death_upper| quantile|
|:-----------|-----------:|----------:|-----------:|--------:|
|country     |    1.292805|  0.0000177|   0.0000010|   0.5000|
|country     |    1.698443|  0.0001428|   0.0000042|   0.8000|
|country     |    2.282170|  0.0002861|   0.0000071|   0.9000|
|country     |    6.389453|  0.0010128|   0.0000184|   0.9900|
|country     |   16.204128|  0.1335427|   0.0000272|   0.9990|
|country     | 1271.292788| 17.3087670|   0.0000381|   0.9999|
|county      |    1.374532|  0.0003080|   0.0000076|   0.5000|
|county      |    1.729290|  0.0007144|   0.0000166|   0.8000|
|county      |    2.122974|  0.0009937|   0.0000242|   0.9000|
|county      |    5.062541|  0.0019803|   0.0000585|   0.9900|
|county      |   15.723140|  0.0040280|   0.0001066|   0.9990|
|county      |  113.161688|  0.0096127|   0.0010251|   0.9999|
|state       |    1.267797|  0.0001721|   0.0000040|   0.5000|
|state       |    1.463159|  0.0004531|   0.0000096|   0.8000|
|state       |    1.623592|  0.0006939|   0.0000132|   0.9000|
|state       |    3.056127|  0.0012846|   0.0000555|   0.9900|
|state       |    6.814701|  0.0017283|         Inf|   0.9990|
|state       |    9.762749|  0.0019752|         Inf|   0.9999|
|subnational |    1.257974|  0.0000877|   0.0000038|   0.5000|
|subnational |    1.536128|  0.0002559|   0.0000077|   0.8000|
|subnational |    1.767536|  0.0003890|   0.0000105|   0.9000|
|subnational |    3.034128|  0.0008125|   0.0000229|   0.9900|
|subnational |    6.408077|  0.0012976|   0.0000566|   0.9990|
|subnational |   13.925664|  0.0015355|         Inf|   0.9999|

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

|Combined_Key          |   rt_upper|
|:---------------------|----------:|
|Lesotho               | 1271.29279|
|Madison, Texas        |  744.29717|
|Medina, Texas         |  314.30089|
|Rapides, Louisiana    |  167.92356|
|Bahamas               |  144.21988|
|Orange, Texas         |  140.73089|
|Carson City, Nevada   |  116.84549|
|South Sudan           |  115.57681|
|Bastrop, Texas        |  110.28951|
|Lincoln, Arkansas     |  108.67827|
|Titus, Texas          |  104.05893|
|Moore, Texas          |  103.19862|
|Woodward, Oklahoma    |   84.16939|
|Columbiana, Ohio      |   83.09680|
|Morgan, Kentucky      |   69.16967|
|Iberville, Louisiana  |   63.44578|
|Burnet, Texas         |   60.71425|
|Sao Tome and Principe |   57.62610|
|Karnes, Texas         |   56.34749|
|Reeves, Texas         |   53.58266|

 </td>
   <td> 

|Combined_Key              | case_upper|
|:-------------------------|----------:|
|Heilongjiang, China       | 18.4129651|
|Hubei, China              |  2.0402384|
|Marion, Ohio              |  0.0838188|
|Livingston, Missouri      |  0.0526679|
|Lee, Arkansas             |  0.0418891|
|Johnson, Wyoming          |  0.0311574|
|La Salle, Texas           |  0.0196973|
|Bent, Colorado            |  0.0126872|
|Crowley, Colorado         |  0.0103265|
|Wallis and Futuna, France |  0.0103020|
|Dakota, Nebraska          |  0.0101501|
|Finney, Kansas            |  0.0098849|
|Crockett, Texas           |  0.0096837|
|Norton, Kansas            |  0.0089420|
|Ralls, Missouri           |  0.0086793|
|Bon Homme, South Dakota   |  0.0086349|
|Lincoln, Arkansas         |  0.0080790|
|Comanche, Kansas          |  0.0080249|
|Pawnee, Kansas            |  0.0077960|
|Dickey, North Dakota      |  0.0075346|

 </td>
   <td> 

|Combined_Key            | death_upper|
|:-----------------------|-----------:|
|Stanly, North Carolina  |   0.0010251|
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
|Kamchatka Krai, Russia  |   0.0001140|
|Lawrence, Alabama       |   0.0001129|
|Montgomery, New York    |   0.0001119|
|Columbia, New York      |   0.0001113|
|Houston, Alabama        |   0.0001108|
|Covington, Alabama      |   0.0001102|
|Marion, Illinois        |   0.0001093|
|Lycoming, Pennsylvania  |   0.0001087|
|Hill, Texas             |   0.0001085|

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
## 0.000e+00 2.000e+00 4.000e+00 2.489e+15 7.000e+00 5.167e+18
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
## 0.00000 0.02185 0.21605 0.35378 0.61861 1.00000
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
## Overdispersion detected in 0.216867 cases.
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
