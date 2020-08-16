# R code for running the sensitivity analysis for HERD parameters
# Author: Corbin Quick, Xihao Li, Derek Shyr, Zilin Li

setwd("/n/holystore01/LABS/xlin/Lab/covid19/analysis/sensitivity/US")
require(data.table)
require(stringr)
require(ggplot2)
require(Matrix)

require(lme4)
require(splines)
require(geepack)
require(broom)

HERD_PARAMS_choices <- data.frame(DAYS_IMMUNE = c(60,60,60,60,60,90,90,90,90,90,365,365,365,365,365),
                                  ASC_RATIO = c(1,2,4,8,16,1,2,4,8,16,1,2,4,8,16))

for (choice in 1:dim(HERD_PARAMS_choices)[1]){
  rm(list=setdiff(ls(), c("HERD_PARAMS_choices", "choice")))
  gc()
  
  # ----------------------------------------------------------------
  #  Fixed parameters for Rt model fitting. 
  # ----------------------------------------------------------------
  
  OUT_TABLE_NAME = paste0("/n/holystore01/LABS/xlin/Lab/covid19/sensitivity/US/HERD/US_Rt_Intervention_Output_","DAYS_IMMUNE=",HERD_PARAMS_choices[choice,"DAYS_IMMUNE"],"_ASC_RATIO=",HERD_PARAMS_choices[choice,"ASC_RATIO"],".csv")
  
  # ** denotes parameters we may want to vary for sensitivity analysis. 
  
  # "Outbreak" parameters.
  OUTBREAK_PARAMS = list(
    MIN_CASES = 50 # ** Number of cases considered an "outbreak".
  )
  
  # "Herd immunity" parameters
  HERD_PARAMS <- list(
    DAYS_IMMUNE = HERD_PARAMS_choices[choice,"DAYS_IMMUNE"], # ** How long does immunity last after infection?
    ASC_RATIO = HERD_PARAMS_choices[choice,"ASC_RATIO"], # ** Ascertainment ratio (total:confirmed)
    MAX_P_INFECTED = 0.98 # Max prop. residents infected (must be < 1)
  )
  
  # Time trend (date) spline parameters 
  SPLINE_PARAMS <- list(
    KNOTS_DATE = 8, # ** Knots for calendar date cubic spline.
    KNOTS_DS = 0 # ** Knots for days-since-outbreak cubic spline.
  )
  
  # Intervention lags
  INTERVENTION_PARAMS <- list(
    DAYS_LAG = 14 # ** How many days before interventions affect Rt?
  )
  
  # Serial interval (SI) parameters. 
  SI_PARAMS = list(
    MEAN = 5.2, # ** SI mean.
    SD = 5.5, # ** SI standard deviation. 
    NTS = 30 # Maximum no. days infectious. 
  )
  
  # Confirmation-infection lag "deconvolution". 
  LAG_PARAMS = list(
    USE_LAG = FALSE, # ** Use the random lag deconvolution?
    MAX_DAYS = 14, # Maximum lag period. 
    SHAPE = 7, # ** Lag distribution: Gamma rate. 
    RATE = 1 # ** Lag distribution: Gamma rate. 
  )
  
  
  # ----------------------------------------------------------------
  #  Define misc. functions.
  # ----------------------------------------------------------------
  
  do_or_dont_die <- function(expr) tryCatch(expr, error = function(x) NULL, finally = function(x) NULL)
  
  center <- function(x) scale(x, scale = FALSE, center = TRUE)
  
  cap <- function(x, cp, f = `<=`) (x - cp)*f(x,cp) + cp
  
  inv <- function(x, xn = names(x)){names(xn) <- x; xn}
  
  # add new terms to a formula object
  add_formula <- function(eqn, x, is.char = FALSE){
    if(!is.char){
      x <- gsub('"', "", deparse(substitute(x)))
    }
    as.formula(paste(as.character(eqn)[2] , '~', as.character(eqn)[3], '+', x))
  }
  
  tidyCoef <- function(x){
    require(data.table)
    out <- as.data.table(tidy(x))
    setnames(out, c("term", "estimate", "se", "stat", "pval"))
    out
  }
  
  # function to sum over sliding windows
  slideSum <- function(x, n, fill = 0){
    x[is.na(x)] <- fill
    csx = cumsum(x)
    csx - data.table::shift(csx, n = n, type = 'lag', fill = fill)
  }
  
  # function to average over sliding windows
  slideMean <- function(x, n, fill = 0, count_na = FALSE){
    if(count_na){
      denom <- n
    }else{
      # number of non-missing points in window
      denom <- slideSum(!is.na(x), n = n, fill = 0)
    }
    x[is.na(x)] <- fill
    csx = cumsum(x)
    out = (csx - data.table::shift(csx, n = n, type = 'lag', fill = fill))/denom
    out[denom == 0] = NA
    out
  }
  
  
  # ----------------------------------------------------------------
  #  Process input data and source code.
  # ----------------------------------------------------------------
  
  dm2 <- readRDS("input_Rt_31Jul20.rds")
  
  source("CalculateCumulativeIncidence.R")
  source("write_fit_docx.R")
  
  setorder(dm2, FIPS, date)
  
  
  # ----------------------------------------------------------------
  #  Confirmation-infection lag time: infection count "deconvolution"
  # ----------------------------------------------------------------
  
  deconvLags <- function(y, lags, weights){
    require(data.table)
    
    if( length(weights) != length(lags)) stop("Mismatch.")
    if( min(weights) < 0 ) stop("Negative weight")
    
    out <- rep(0, length(y))
    
    weights <- weights/sum(weights)
    
    for(i in 1:length(weights)){
      if(weights[i] > 0){
        out <- out + weights[i] * shift(y, n = -lags[i], fill = NA)
      }
    }
    out
  }
  
  getGammaLags <- function(max_days = LAG_PARAMS$MAX_DAYS, shape = LAG_PARAMS$SHAPE, rate = LAG_PARAMS$RATE){
    l <- 0:max_days
    w <- pgamma(l + 1, shape, rate) - pgamma(l, shape, rate)
    data.table(l,w)
  }
  
  deconvGammaLags <- function(y, max_days = LAG_PARAMS$MAX_DAYS, shape = LAG_PARAMS$SHAPE, rate = LAG_PARAMS$RATE){
    l <- 0:max_days
    w <- pgamma(l + 1, shape, rate) - pgamma(l, shape, rate)
    w <- w/sum(w)
    deconvLags(y, l, w)
  }
  
  if( LAG_PARAMS$USE_LAG ){
    dm2$reported_new_cases <- as.numeric(dm2$new_cases)
    dm2$reported_positive <- as.numeric(dm2$positive)
    dm2$new_cases <- as.numeric(dm2$new_cases + 1e-10)
    dm2$positive <- as.numeric(dm2$positive + 1e-10)
    dm2 <- dm2[order(FIPS,date)]
    dm2[,new_cases := deconvGammaLags(new_cases),by=FIPS]
    dm2[,positive := cumsum(new_cases),by=FIPS]
  }
  
  dm2 <- subset(dm2, !is.na(new_cases))[order(FIPS,date)]
  
  
  # ----------------------------------------------------------------
  #  Calculate Infectious Potential, "Lambda"
  # ----------------------------------------------------------------
  
  # Calculate the cumulative daily incidence
  dm2[,Lam := getLambda(
    x = new_cases,
    t = days_int, 
    mu = SI_PARAMS$MEAN, sd = SI_PARAMS$SD, NTS = SI_PARAMS$NTS
  ),by=FIPS]
  
  
  # ----------------------------------------------------------------
  #  Calculate days since "outbreak" within each county
  # ----------------------------------------------------------------
  
  dm2[,
      days_s_outbreak := if( max(positive) < OUTBREAK_PARAMS$MIN_CASES ){
        as.numeric(0)
      }else{
        cap(as.numeric(date - min(date[positive >= OUTBREAK_PARAMS$MIN_CASES])), 0, `>=`)
      }
      ,by = FIPS]
  
  
  # ----------------------------------------------------------------
  #  Update intervention and mobility metric covariates
  # ----------------------------------------------------------------
  
  dm2 <- dm2[order(FIPS,date)]
  
  # Baseline mobility levels in each county
  dm2[,sg_left_home_baseline := mean(sg_left_home[date <= "2020-03-10" & date >= "2020-02-01"]),by=FIPS]
  dm2[,sg_time_home_baseline := mean(sg_time_home[date <= "2020-03-10" & date >= "2020-02-01"]),by=FIPS]
  
  # Lagged intervention terms
  CSP_names <- c('mask_score_CSP','business_score_CSP','bar_score_CSP','rest_score_CSP','gatherings_CSP','sho_CSP','bar_rr_CSP')
  CSPb_names <- paste0(CSP_names, 'b')
  
  dm2[,
      (paste0(c(CSP_names, CSPb_names), "_lag")) := lapply(.SD, shift, n = INTERVENTION_PARAMS$DAYS_LAG, fill = 0)
      ,by=FIPS, .SDcols = c(CSP_names, CSPb_names)]
  
  
  # ----------------------------------------------------------------
  #  Calculate % population "immune" to infection in each county
  # ----------------------------------------------------------------
  
  dm2[,
      n_immune := shift(slideSum(new_cases, HERD_PARAMS$DAYS_IMMUNE), n =1, fill = 0)
      ,by= FIPS]
  
  dm2[,
      p_immune := cap(HERD_PARAMS$ASC_RATIO * n_immune, HERD_PARAMS$MAX_P_INFECTED*pop_size)/pop_size
      ,by= FIPS]
  
  
  # ----------------------------------------------------------------
  #  Function to estimate Rt regression model
  # ----------------------------------------------------------------
  
  # Fit overdispersed Poisson model using geepack::geeglm. 
  # Betas from this model can be interpreted as 
  # differences in log(Rt). 
  
  fitPoissonRt <- function(eqn, data, idvar = "FIPS", beta0 = NULL, lag_variables = c(), lag_size = 0, lag_fill = 0, min_positives = 0, min_Lam = 0){
    
    if( length(lag_variables) > 0 & lag_size > 0 ){
      
      if( var(data[,length(date),by=FIPS]$V1) > 0 ){
        stop("Incomplete dates for some counties. Cannot lag variables after filtering by date.")
      }
      
      data <- as.data.table(as.data.frame(data))[order(FIPS,date)]
      data[,(paste0(lag_variables, "_lag")) := lapply(.SD, shift, n = lag_size, fill = lag_fill),.SDcols = lag_variables,by = FIPS]
    }
    
    data = subset(data, positive >= min_positives & Lam >= min_Lam)
    
    lm_fit <- lm(eqn, data = data, na.action = na.exclude)
    
    # Need to exclude observations with any 
    # missing values before fitting GEE.
    data <- subset(data, !is.na(predict(lm_fit)))
    data[["ID"]] <- as.integer(factor(data[[idvar]], order = TRUE))
    data <- data[order(ID,date)]
    
    geepack::geeglm( 
      eqn, 
      data = data, 
      id = ID, 
      family = poisson(),
      #		start = beta0,
      corstr = "ar1",
      control = geese.control("maxit" = 50000, trace = TRUE)
    )
  }
  
  
  # ----------------------------------------------------------------
  #  Define covariate names
  # ----------------------------------------------------------------
  
  # Nice covariate names
  nice_names <- c(
    "sg_left_home_baseline" = "SafeGraph % Left Home (Baseline)",
    "sg_left_home" = "SafeGraph % Left Home",
    "sho_CSPb" = "Stay-Home Mandate",
    "mask_score_CSPb" = "Mask Mandate",
    "business_score_CSPb" = "Non-Essential Business Closure",
    "gatherings_CSPb" = "Social Gathering Restriction",
    "bar_rr_CSPb" = "Bar and Restaraunt Closure",
    "pop_size" = "County Population Size",
    "HouseholdSize" = "Average Household Size",
    "Pct_Afr" = "% African American",
    "Pct_Lat" = "% Latinx American",
    "Pct_Asn" = "% Asian American",
    "Pct_HighSchool" = "% Highschool Graduates",
    "p_urban" = "Pr. Urban Population",
    "p_poverty" = "Pr. Poverty"
  )
  
  rm_prefixes <- c("^stateName", "^weekday", "^division", "^region")
  rm_suffixes <- c("_s[0-9]*$", "_l[0-9]*$", "_lag$")
  
  niceNames <- function(x){
    for(i in c(rm_prefixes, rm_suffixes)){
      x <- gsub(i, "", x)
    }
    original <- x
    x <- gsub(".*\\(", "", x)
    x <- gsub("\\).*", "", x)
    wr <- x %in% names(nice_names)
    x[wr] <- nice_names[x[wr]]
    x[!wr] <- original[!wr]
    x <- gsub("_", " ", x)
    x
  }
  
  
  # ----------------------------------------------------------------
  #  Fit models
  # ----------------------------------------------------------------
  
  csp_ivs <- paste0(c("sho_CSPb", "mask_score_CSPb", "business_score_CSPb", "gatherings_CSPb", "bar_rr_CSPb"), "_lag")
  
  eqn_base <- new_cases ~ offset(log1p(Lam)) + offset(log(1 - p_immune)) + 
    center(sg_left_home_baseline) + 
    center(log(pop_size)) + center(HouseholdSize) +
    center(Pct_Afr) + center(Pct_Lat) + center(Pct_Asn) + center(Pct_HighSchool) + 
    center(p_urban) + center(p_poverty) + 
    stateName + bs(days_s_outbreak, df = 3 + SPLINE_PARAMS$KNOTS_DS) +  
    weekday + bs(date, df = 3 + SPLINE_PARAMS$KNOTS_DATE) 
  
  eqn_wmob <- add_formula(eqn_base, "center(sg_left_home_s7)", TRUE)
  
  joint_fit <- fitPoissonRt(
    add_formula(eqn_base, paste(csp_ivs, collapse = ' + '), TRUE)
    , data = dm2, min_positive = OUTBREAK_PARAMS$MIN_CASES)
  
  marginal_fits <- lapply(csp_ivs, function(x){
    fitPoissonRt(add_formula(eqn_base, x, TRUE), data = dm2, min_positive = OUTBREAK_PARAMS$MIN_CASES)
  })
  names(marginal_fits) <- csp_ivs
  
  joint_wmob_fit <- fitPoissonRt(
    add_formula(eqn_wmob, paste(csp_ivs, collapse = ' + '), TRUE)
    , data = dm2, min_positive = OUTBREAK_PARAMS$MIN_CASES)
  
  marginal_wmob_fits <- lapply(csp_ivs, function(x){
    fitPoissonRt(add_formula(eqn_wmob, x, TRUE), data = dm2, min_positive = OUTBREAK_PARAMS$MIN_CASES)
  })
  names(marginal_wmob_fits) <- csp_ivs
  
  
  fit_list <- list()
  ii <- 1
  
  fit_list[[ii]] <- tidyCoef(joint_fit)
  fit_list[[ii]]$model <- "Interventions"
  fit_list[[ii]]$type <- "Joint"
  ii <- ii + 1
  
  for(x in csp_ivs){
    fit_list[[ii]] <- tidyCoef(marginal_fits[[x]])
    fit_list[[ii]]$model <- "Interventions"
    fit_list[[ii]]$type <- "Marginal"
    ii <- ii + 1
  }
  
  fit_list[[ii]] <- tidyCoef(joint_wmob_fit)
  fit_list[[ii]]$model <- "Interventions + Mobility"
  fit_list[[ii]]$type <- "Joint"
  ii <- ii + 1
  
  for(x in csp_ivs){
    fit_list[[ii]] <- tidyCoef(marginal_wmob_fits[[x]])
    fit_list[[ii]]$model <- "Interventions + Mobility"
    fit_list[[ii]]$type <- "Marginal"
    ii <- ii + 1
  }
  
  fit_dt <- rbindlist(fit_list)
  
  fit_dt[,Term := niceNames(term),]
  
  fwrite(fit_dt, OUT_TABLE_NAME, sep = ',', col.names = TRUE, quote =TRUE)
  
  # ggplot(subset(fit_dt, term %in% csp_ivs),
  # aes(y = estimate, ymin = estimate - 2* se, ymax = estimate + 2* se, x = Term, colour = paste(type, model))
  # ) + geom_hline(yintercept = 0, colour = 'black') + geom_point(position = position_dodge(width=0.5)) + geom_errorbar(width = 0, position = position_dodge(width=0.5)) + coord_flip() + xlab(NULL) + ylab(NULL) + theme_minimal() + guides(colour = guide_legend(title=NULL))
  
  
}

