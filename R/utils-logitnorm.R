

#' The logit function
#'
#' @param x a number between 0 and 1
#'
#' @return a number between -Inf and Inf
#' @export
logit = function(x) {
  log(x/(1-x))
}

#' The inverse logit function
#'
#' @param y a number between -Inf and Inf
#'
#' @return a number between 0 and 1
#' @export
inv_logit = function(y) {
  1/(1+exp(-y))
}

# # adapted from logitnorm package
# twCoefLogitnorm = function (median, quant, perc = 0.975, ...) {
#   mu = logit(median)
#   upperLogit = logit(quant)
#   sigmaFac = qnorm(perc)
#   sigma = 1/sigmaFac * (upperLogit - mu)
#   tibble::tibble(mu = mu, sigma = sigma)
# }
# 
# # adapted from logitnorm package
# twCoefLogitnormCi = function (lower, upper, perc = 0.975, isTransScale = FALSE) {
#   sigmaFac = qnorm(perc)
#   if (!isTRUE(isTransScale)) {
#     lower <- logit(lower)
#     upper <- logit(upper)
#   }
#   halfWidth <- (upper - lower)/2
#   sigma <- halfWidth/sigmaFac
#   tibble::tibble(mu = upper - halfWidth, sigma = sigma)
# }

#' Generate mu and sigma parameters for a logitnormal distribution
#' 
#' The resulting logitnorm distribution will have a set median. The confidence intervals
#' will not match those provided as they are used as a inter-quartile range.
#'
#' @param median the median of the
#' @param lower the lower CI
#' @param upper the upper CI
#' @param ci the confidence limits
#' @param fix_median make the median of the logitnorm be the same as the median given.
#'   This can cause issues when very skewed distributions are used
#' @param ... not used
#'
#' @return a tibble with mu and sigma columns
#' @export
ci_to_logitnorm = function (median, lower, upper, ci = 0.95, fix_median = TRUE, ...) {
  pkgutils::recycle(median,lower,upper, fix_median)
  pkgutils::check_consistent(lower>0, median>=lower, upper>=median, upper<1)
  lower <- logit(lower)
  upper <- logit(upper)
  iqr = upper - lower
  mu = ifelse( fix_median, logit(median), upper-iqr/2 )
  sigmaFac = stats::qnorm(1-(1-ci)/2)
  sigma = 1/sigmaFac * iqr/2
  tibble::tibble(mu = mu, sigma = sigma)
}


beta_dist_to_logitnorm = function(beta_dist) {
  tibble::as_tibble(beta_dist) %>% 
    dplyr::mutate(
      lower = dplyr::if_else(lower==0, .Machine$double.eps, lower), 
      upper = dplyr::if_else(upper==1, 1-.Machine$double.eps, upper),
      median = dplyr::case_when(
        median<lower ~ lower,
        median>upper ~ upper,
        TRUE ~ median
      ), 
      ci_to_logitnorm(median,lower,upper, fix_median=(shape1>1 & shape2>1))
    ) %>%
    dplyr::select(mu,sigma)
}

.mean_logitnorm = function(mu, sigma) {
  pkgutils::recycle(mu,sigma)
  pkgutils::check_consistent(sigma>=0)
  sapply(seq_along(mu), function(i) {
    f = function(x) x*.d_logitnorm(x,mu[i],sigma[i])
    out = stats::integrate(f,0,1)
    out$value
  })
}

.logit_quantiles = function(mu, sigma, p=c(median=0.5,lower=0.025,upper=0.975)) {
  n = pkgutils::recycle(mu, sigma)
  if (is.null(names(p))) {
    names(p) = sprintf("Q.%1.3g",p)
  }
  dplyr::bind_cols(lapply(seq_along(p), function(i) {
    qs = .q_logitnorm(p[i],mu,sigma)
    tibble::as_tibble(as.list(qs))
  }))
}

.p_logitnorm = function(q, mu = 0, sigma = 1, ...) {
  pkgutils::recycle(q,mu,sigma)
  pkgutils::check_consistent(sigma>=0)
  ql <- logit(q)
  stats::pnorm(ql, mean = mu, sd = sigma, ...)
}

.d_logitnorm = function(x, mu = 0, sigma = 1, log = FALSE, ...) {
  pkgutils::recycle(x,mu,sigma)
  pkgutils::check_consistent(sigma>=0)
  ql <- logit(x)
  if (log) {
    ifelse(x <= 0 | x >= 1, 0, stats::dnorm(ql, mean = mu, sd = sigma, 
                                     log = TRUE, ...) - log(x) - log1p(-x))
  }
  else {
    ifelse(x <= 0 | x >= 1, 0, stats::dnorm(ql, mean = mu, sd = sigma, 
                                     ...)/x/(1 - x))
  }
}

.q_logitnorm = function(p, mu = 0, sigma = 1, ...) {
  pkgutils::recycle(p,mu,sigma)
  pkgutils::check_consistent(sigma>=0)
  qn <- stats::qnorm(p, mean = mu, sd = sigma, ...)
  inv_logit(qn)
}

.r_logitnorm = function(n, mu = 0, sigma = 1, ...) {
  pkgutils::recycle(n,mu,sigma)
  pkgutils::check_integer(n)
  pkgutils::check_consistent(sigma>=0)
  inv_logit(stats::rnorm(n = n, mean = mu, sd = sigma, ...))
}