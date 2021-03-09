gee_fit <- function(eqn, dt, gee_id_var, corstr = "exchangeable"){
  dt <- as.data.frame(dt)
  dt$gee_id_var <- as.integer(factor(dt[[gee_id_var]], order = TRUE))

  dt <- as.data.table(dt)
  setorder(dt, gee_id_var)
  fit = geepack::geeglm(formula = eqn, data = dt,
                        id = gee_id_var, corstr = corstr,
                        family =  poisson(link = "log"))
  # fit = geepack::geese(formula = eqn, data = dt,
  #                       id = gee_id_var, corstr = corstr,
  #                       family =  poisson(link = "log"))
  return(fit)
}


reg_coef_data <- function(fit, se_name = "Std.err", pattern="scale"){
  coef_table <- as.data.table(coef(summary(fit)), keep = "Variable")

  coef_table <- coef_table %>% filter(Variable %like% pattern |
                                        Variable %like% "factor" |
                                        Variable %like% "lock_down"|
                                        Variable %like% "log")
  y <- coef_table[, "Estimate"]
  y_se <- coef_table[, se_name]
  library(magrittr)
  coef_name  = coef_table$Variable %<>% gsub(pattern = "scale",
                                             replacement="") %>%
    gsub(pattern = "weekday",
         replacement="") %>%
    gsub(pattern = "continent",
         replacement="") %>%
    gsub(pattern = "factor",
         replacement="") %>%
    gsub(pattern = "\\(",
         replacement="") %>%
    gsub(pattern = "\\)",
         replacement="") %>%
    gsub(pattern = "\\)",
         replacement="")
  colnames(coef_table)[5] <- "P"
  P = coef_table %>% select(P)

  CI <- paste0("[",round(y-1.96*y_se,2),",",round(y+1.96*y_se,2),"]")



  plot_dt <- data.frame(name=coef_name,
                        value=y,
                        se=y_se,
                        CI= CI,
                        P = P)
  return(plot_dt)
  return(plot_dt)
}
