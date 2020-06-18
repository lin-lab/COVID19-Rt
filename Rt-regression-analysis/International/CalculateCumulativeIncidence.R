
# the following is taken directly from EpiEstim:
vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}

# serial interval function, taken directly from EpiEstim:
discr_si <- function(k, mu, sigma) 
{
  if (sigma < 0) {
    stop("sigma must be >=0.")
  }
  if (mu <= 1) {
    stop("mu must be >1")
  }
  if (any(k < 0)) {
    stop("all values in k must be >=0.")
  }

  a <- ((mu - 1) / sigma)^2
  b <- sigma^2 / (mu - 1)

  cdf_gamma <- function(k, a, b) stats::pgamma(k, shape = a, scale = b)

  res <- k * cdf_gamma(k, a, b) + 
    (k - 2) * cdf_gamma(k - 2, a, b) - 2 * (k - 1) * cdf_gamma(k - 1, a, b)
  res <- res + a * b * (2 * cdf_gamma(k - 1, a + 1, b) - 
                          cdf_gamma(k - 2, a + 1, b) - cdf_gamma(k, a + 1, b))
  res <- vnapply(res, function(e) max(0, e))

  return(res)
}

getLambda <- function(x, t, mu, sd, NTS = 30){

	# x : daily increase in number of positives
	# t : date or integer time

	ti = as.integer(t)

	to = t - min(t) + 1

	if( !all( to == sort(to) ) ) stop("Input must be sorted")

	# full range of dates, in case some are missing
	tt = min(to):max(to)

	# fill in missing values with zeros
	xx = 0*tt
	xx[to] <- x

	# calculate cumulative incidence
	Lambda = colSums(do.call(rbind, lapply(1:NTS, function(i) discr_si(i,mu,sd)*data.table::shift(xx, n = i, fill = 0))))
	
	# map these back to dates present in input
	Lambda = Lambda[to]
	
	# return
	Lambda
}



