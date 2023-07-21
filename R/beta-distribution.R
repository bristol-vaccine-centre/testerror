#' Generate concave beta distribution parameters from mean and confidence intervals
#'
#' @param mean the mean of the probability given
#' @param lower the lower ci of the probability given 
#' @param upper the upper ci of the probability given 
#' @param confint the ci limits 
#' @param widen widen the spread of the final beta by this factor.
#' @param limit the lowest posible value for the shape parameters of the resulting
#'   `beta_dist`
#'
#' @return a list with shape1, shape2 values, and d, p, q and r functions
#' @export
#'
#' @examples
#' beta = beta_params(0.25, 0.1, 0.3)
#' tmp = beta$r(n=1000)
#' diff = stats::quantile(tmp,0.975)-stats::quantile(tmp,0.025)
#' if(abs(diff-0.2) > 0.1) stop("confidence limits should be 0.2 apart")
#' if(abs(mean(tmp)-0.25) > 0.1) stop("mean value should be 0.25")
#' 
beta_params = function( mean, lower, upper, confint = 0.95, widen = 1, limit=1,  ...) {
  # find the effective size that gives an 95% CI 
  # k is a concentration parameter.
  
  zcrit = (1-confint)/2
  nsens = 
    tryCatch(
      stats::uniroot(f = function(k) {
        # shape1+shape2 = k
        # mean = shape1/(shape1+shape2)
        # mean * k = shape1
        # shape2 = k*(1-mean)
        # TODO: should mean be median in fact.
        # (median * (k-2/3)) - 1/3 = shape1
        # shape2 = k-((median * (k-2/3)) - 1/3)
        (
          stats::qbeta(1-zcrit, shape1 = mean*k, shape2 = (1-mean)*k)-
            stats::qbeta(zcrit, shape1 = mean*k, shape2 = (1-mean)*k)
        )-(
          upper-lower
        )}, 
        interval = c(2,10000000) 
      )$root,
      error = function(e) return(2))
  
  # widen CIs by reducing concentration parameter
  nsens = nsens * 1/widen
  ksens = nsens*mean
  
  if (ksens < limit) {
    nsens = nsens * limit/ksens
    ksens = limit
  }
  
  if ((nsens-ksens) < limit) {
    nsens = nsens * limit/(nsens-ksens)
    ksens = nsens-limit
  }
  
  return(beta_dist(
      p = ksens,
      q = nsens-ksens))
}



#' Fit a beta distribution to data using method of moments
#'
#' @param samples a set of probabilities
#' @param na.rm should we ignore NA values
#'
#' @return a `beta_dist` S3 object fitted to the data.
#' @export
#'
#' @examples
#' beta_fit(rbeta(10000,40,60))
#' beta_fit(rbeta(10000,1,99))
beta_fit = function(samples, na.rm=FALSE) {
  if (any(samples<0 | samples>1)) stop("samples out of range")
  if (na.rm) samples = na.omit(samples)
  # tmp = try(suppressWarnings(
  #   MASS::fitdistr(x = samples,densfun = "beta", start = list(shape1 = mean(samples),shape2 = 1-mean(samples)))
  # ),silent = TRUE)
  # if (isa(tmp,"try-error")) {
    v = sd(samples)^2
    e = mean(samples)
    return(beta_dist(
      ((e*(1-e))/v-1)*e,
      ((e*(1-e))/v-1)*(1-e)
    ))
  # } else {
  #   return(
  #     beta_dist(tmp$estimate["shape1"],tmp$estimate["shape2"])
  #   )
  # }
}

#' Generate a beta distribution out of probabilities, or positive and negative counts
#' 
#' @param p the first shape / the probability or count of success
#' @param q (optional) the second shape / the probability or count of failure
#' @param n (optional) the number of trials.
#' @param ... 
#'
#' @return either a single `beta_dist` object or a list of `beta_dist`s
#'
#' @export
#' @examples 
#' beta_dist(c(1,2,3),c(3,2,1))
beta_dist = function(p=n-q, q=n-p, n=p+q,...) {
  if (length(p) != length(q) || length(p) != length(n)) stop("all inputs must be the same length")
  if (length(p) > 1) {
    tmp = lapply(1:length(p), function(i) beta_dist(p=p[[i]], q=q[[i]], n=n[[i]], ...))
    return(structure(tmp, class=c("beta_dist_list",class(tmp))))
  }
  if (p<1 && q==floor(q) && q>1) message("beta_dist: second parameter could be interpreted as number of trials, did you mean to use `n=`")
  if (n != p+q) {
    shape1 = p/(p+q)*n
    shape2 = q/(p+q)*n
  } else {
    shape1 = p
    shape2 = q
  }
  return(structure(
    list(
      shape1 = shape1,
      shape2 = shape2,
      conc = shape1+shape2,
      d = function(x) stats::dbeta(x, shape1, shape2),
      p = function(q) stats::pbeta(q, shape1, shape2),
      q = function(p) stats::qbeta(p, shape1, shape2),
      r = function(n) stats::rbeta(n, shape1, shape2)
    ),
    class = "beta_dist"))
}

#' convert a beta distribution to a tibble
#' 
#' @param x the beta distribution
#' @param prefix name to output columns prefix.lower, prefix.upper etc
#' @param confint confidence intervals
#' @param ... not used
#' @importFrom tibble as_tibble
#' 
#'
#' @export
as_tibble.beta_dist = function(x, prefix=NULL, confint = 0.95, ...) {
  tmp = tibble(
    shape1 = x$shape1,
    shape2 = x$shape2,
    conc = x$conc,
    mean = x$shape1/x$conc,
    upper = x$q(1-(1-confint)/2),
    lower = x$q((1-confint)/2),
  ) 
  if (!is.null(prefix)) 
    tmp = tmp %>% rename_with(~ sprintf("%s.%s",prefix,.x))
  return(tmp)
}

#' Format a beta distribution
#' 
#' @param x the beta distribution
#' @param ... not used
#' 
#' @return nothing
#'
#' @export
format.beta_dist = function(x, ...) {
  do.call(sprintf, c(list(fmt = "%1.3f [%1.3f\u2013%1.3f] (N=%1.2f)"), x$shape1/ (x$shape1 + x$shape2), 
                     qbeta(c(0.025,0.975), shape1 = x$shape1, shape2 = x$shape2),x$shape1 + x$shape2
  ))
}

#' Print a beta distribution
#' 
#' @param x the beta distribution
#' @param ... not used
#' @return nothing
#' @export
print.beta_dist = function(x, ...) cat(format(x, ...), "\n")