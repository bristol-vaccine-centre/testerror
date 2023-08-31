#' Generate concave beta distribution parameters from mean and confidence intervals
#'
#' @param median the median of the probability given
#' @param lower the lower ci of the probability given 
#' @param upper the upper ci of the probability given 
#' @param confint the ci limits 
#' @param widen widen the spread of the final beta by this factor
#' @param limit the lowest possible value for the shape parameters of the resulting
#'   `beta_dist` (1 enforces that the distribution is convex)
#' @param ... not used
#'
#' @return a list with shape1, shape2 values, and d, p, q and r functions
#' @export
#'
#' @examples
#' beta = beta_params(0.25, 0.1, 0.3)
beta_params = function( median, lower, upper, confint = 0.95, widen = 1, limit=1,  ...) {
  # find the effective size that gives an 95% CI 
  # k is a concentration parameter.
  
  zcrit = (1-confint)/2
  
  data = tibble::tibble(
    x = c(0.5,zcrit,1-zcrit),
    y = c(median,lower,upper)
  )
  
  initc = max(1/c(median,1-median))
  
  tmp = stats::nls(y ~ stats::qbeta(x, p*conc, (1-p)*conc), data = data, 
                   start=list(p= median, conc=initc+1), 
                   lower=list(p=0, conc=initc), 
                   upper=list(p=1, conc=Inf),
                   algorithm="port")
  p = stats::coef(tmp)[["p"]]
  conc = stats::coef(tmp)[["conc"]]/widen
  
  if (p*conc < limit) {
    conc = limit/p
  }
  
  if ((1-p)*conc < limit) {
    conc = limit/(1-p)
  }
  
  return(beta_dist(p = p, n=conc))
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
#' beta_fit(stats::rbeta(10000,40,60))
#' beta_fit(stats::rbeta(10000,1,99))
beta_fit = function(samples, na.rm=FALSE) {
  if (any(samples<0 | samples>1)) stop("samples out of range")
  if (na.rm) samples = stats::na.omit(samples)
  
  v = stats::sd(samples)^2
  e = mean(samples)
  return(beta_dist(
    shape1=((e*(1-e))/v-1)*e,
    shape2=((e*(1-e))/v-1)*(1-e)
  ))

}

#' Generate a beta distribution out of probabilities, or positive and negative counts
#' 
#' @param p the first shape / the probability or count of success
#' @param q (optional) the second shape / the probability or count of failure
#' @param n (optional) the number of trials.
#' @param shape1 the first shape parameter (use this to force interpretation as shape)
#' @param shape2 the second shape parameter (use this to force interpretation as shape)
#' @param ... not used
#'
#' @return either a single `beta_dist` object or a list of `beta_dist`s
#'
#' @export
#' @examples 
#' beta_dist(shape1 = c(1,2,3),shape2 = c(3,2,1))
#' beta_dist(p = 0.7, n = 2)
beta_dist = function(..., p=NULL, q=NULL, n=NULL, shape1=NULL, shape2=NULL) {
  # TODO: change this to be explicit about shape versus probability
  if (length(rlang::list2(...))>0) stop("additional parameters detected to `beta_dist`: all parameters must be named")
  len = pkgutils::recycle(p,q,n,shape1,shape2)
  pkgutils::resolve_missing(
    n = shape1+shape2,
    p = shape1/n,
    q = shape2/n,
    shape1 = n*p,
    shape2 = n*q,
    p = 1-q,
    q = 1-p
  )
  pkgutils::check_consistent(
    n >= shape1,
    n >= shape2,
    shape1 >= 0,
    shape2 >= 0,
    p <= 1,
    p >= 0,
    q <= 1,
    q >= 0,
    n >= 0
  )
  
  # if (is.null(shape1) || is.null(shape2)) {
  #   if (is.null(p)+is.null(q)+is.null(n) > 1) stop("two of shape1 (p), shape2 (q) or concentration (n) must be given")
  #   if (is.null(n)) {
  #     # concentration parameter is not given
  #     # p & q are pure shape parameters
  #     shape1 = p
  #     shape2 = q
  #   } else {
  #     if (is.null(p)) {
  #       if (any(q > 1)) {
  #         # q is a shape param
  #         shape1 = n-q
  #         shape2 = q
  #       } else {
  #         # q is either a probability or a small shape number
  #         # we are assuming a probability, because shape2 is NULL
  #         shape1 = n*(1-q)
  #         shape2 = n*q
  #       }
  #     } else if (is.null(q)) {
  #       if (any(p > 1)) {
  #         # p is a shape param
  #         shape1 = p
  #         shape2 = n-p
  #       } else {
  #         # p is either a probability or a small shape number
  #         # we are assuming a probability, because shape2 is NULL
  #         shape1 = n*p
  #         shape2 = n*(1-p)
  #       }
  #     } else {
  #       # all 3 provided
  #       if (all(p+q == 1)) {
  #         shape1 = p*n
  #         shape2 = q*n
  #       } else if (all(p+q == n)) {
  #         shape1 = p
  #         shape2 = q
  #       } else {
  #         stop("shape parameters (p, q) must add up to 1 or n")
  #       }
  #     }
  #   }
  # }
  
  return(as.beta_dist(shape1, shape2))
}


as.beta_dist = function(shape1,shape2) {
  n=pkgutils::recycle(shape1,shape2)
  if (n>1) return(as.beta_dist_list(shape1,shape2))
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

as.beta_dist_list = function(x, ...) {
  UseMethod("as.beta_dist_list", x)
}

as.beta_dist_list.beta_dist_list = function(x, ...) {
  return(x)
}

as.beta_dist_list.beta_dist = function(x, ...) {
  return(as.beta_dist_list.list(list(x)))
}

as.beta_dist_list.list = function(x, ...) {
  tmp = class(x)
  if (!all(sapply(x, isa, "beta_dist"))) stop("list is not only of `beta_dist` objects")
  class(x)<-unique(c("beta_dist_list","beta_dist",tmp))
  return(x)
}

as.beta_dist_list.numeric = function(x, shape2, ...) {
  n=pkgutils::recycle(x,shape2)
  tmp = lapply(1:n, function(i) as.beta_dist(shape1=x[[i]], shape2=shape2[[i]]))
  return(as.beta_dist_list.list(tmp))
}


#' Update the posterior of a `beta_dist`
#'
#' @param x a `beta_dist` or `beta_dist_list` acting as the prior
#' @param pos positive observation(s)
#' @param neg negative observation(s)
#' @param n number observations
#' @param ... not used
#'
#' @return a new `beta_dist` o `beta_dist_list`
#' @export
#'
#' @examples
#' update_posterior(beta_dist(shape1=1,shape2=1), neg=10, n=30)
update_posterior = function(x, ..., pos=NULL, neg=NULL, n=NULL) {
  UseMethod("update_posterior", x)
}

#' @inherit update_posterior
#' @export
update_posterior.beta_dist = function(x, ..., pos=NULL, neg=NULL, n=NULL) {
  pkgutils::recycle(pos,neg,n)
  if (is.null(pos) && is.null(neg) && is.null(n)) return(x)
  pkgutils::resolve_missing(pos=n-neg, neg=n-pos, n=pos+neg)
  
  if (length(pos) > 1) {
    tmp = lapply(seq_along(pos), function(i) {update_posterior.beta_dist(x, pos=pos[[i]], neg=neg[[i]], n=n[[i]])})
    return(structure(tmp, class=c("beta_dist_list","beta_dist",class(tmp))))
  }
  
  return(beta_dist(shape1=x$shape1+pos, shape2=x$shape2+neg))
}

#' @inherit update_posterior
#' @export
update_posterior.beta_dist_list = function(x, ..., pos=NULL, neg=NULL, n=NULL) {
  pkgutils::recycle(x, pos, n)
  if (is.null(pos) && is.null(neg) && is.null(n)) return(x)
  pkgutils::resolve_missing(pos=n-neg, neg=n-pos, n=pos+neg)
  
  tmp = lapply(seq_along(x), function(i) {
    update_posterior(x[[i]], pos = pos[[i]], neg = neg[[i]], n = n[[i]])
  })
  return(structure(tmp, class=c("beta_dist_list","beta_dist",class(tmp))))
}

#' Repeat a `beta_dist`
#'
#' @param x a `beta_dist`
#' @param times n
#' @param ... not used
#'
#' @return a `beta_dist_list`
#' @export
rep.beta_dist = function(x, times, ...) {
  if (times==1) return(x)
  tmp = rep(list(x), times=times, ...)
  return(as.beta_dist_list.list(tmp))
}

#' convert a list of betas to a tibble
#'
#' @param x a beta dist list
#' @inheritDotParams as_tibble.beta_dist
#' @importFrom tibble as_tibble
#'
#' @return a tibble
#' @export
as_tibble.beta_dist_list = function(x, ...) {
  dplyr::bind_rows(lapply(x, as_tibble, ...))
}

#' convert a beta distribution to a tibble
#' 
#' @param x the beta distribution
#' @param prefix name to output columns prefix.lower, prefix.upper etc
#' @param confint confidence intervals
#' @param ... not used
#' @importFrom tibble as_tibble
#' 
#' @export
as_tibble.beta_dist = function(x, prefix=NULL, confint = 0.95, ...) {
  tmp = tibble::tibble(
    shape1 = x$shape1,
    shape2 = x$shape2,
    conc = x$conc,
    mean = x$shape1/x$conc,
    median = x$q(0.5),
    upper = x$q(1-(1-confint)/2),
    lower = x$q((1-confint)/2)
  ) 
  if (!is.null(prefix)) 
    tmp = tmp %>% dplyr::rename_with(~ sprintf("%s.%s",prefix,.x))
  return(tmp)
}

.default_beta_dist_format = function() {
  return(getOption("beta_dist.format",default = "{sprintf('%1.1f%% [%1.1f%%\u2014%1.1f%%] (N=%1.1f)',median*100,lower*100,upper*100,conc)}"))
}

#' Format a beta distribution
#' 
#' @param x the beta distribution
#' @param glue a glue spec taking any of `shape1`, `shape2`, `conc`, `mean`, `median`, `upper`, `lower`
#' @param ... not used
#' 
#' @return nothing
#' 
#' @export
#' @examples 
#' format(beta_dist(shape1=3,shape2=6), "{format(mean*100, digits=3)}%")
format.beta_dist = function(x, glue = .default_beta_dist_format(), ...) {
  tmp = as_tibble.beta_dist(x)
  glue::glue_data(tmp, glue)
  # suppressWarnings(do.call(sprintf, c(list(fmt = fmt), x$shape1/ (x$shape1 + x$shape2), 
  #                    stats::qbeta(c(0.025,0.975), shape1 = x$shape1, shape2 = x$shape2),x$conc
  # )))
}

#' Format a beta distribution list
#' 
#' @param x the beta distribution list
#' @inheritDotParams format.beta_dist
#' 
#' @return nothing
#'
#' @export
format.beta_dist_list = function(x, ...) {
  sapply(x, format.beta_dist, ...)
}

#' Detect the length of a beta distribution
#' 
#' @param x the beta distribution
#' @param ... not used
#' 
#' @return always 1
#'
#' @export
length.beta_dist = function(x, ...) {
  return(1)
}

#' Detect the length of a beta distribution list
#' 
#' @param x the beta distribution list
#' @param ... not used
#' 
#' @return the length of the list
#'
#' @export
length.beta_dist_list = function(x, ...) {
  return(sum(sapply(x, length)))
}

#' Print a beta distribution
#' 
#' @param x the beta distribution
#' @param ... not used
#' @return nothing
#' @export
print.beta_dist = function(x, ...) {
  cat(format(x, ...), "\n")
}

#' Print a beta distribution
#' 
#' @param x the beta distribution
#' @param ... not used
#' @return nothing
#' @export
print.beta_dist_list = function(x, ...) {
  cat(paste0(format(x, ...), collapse="\n"),"\n")
}


#' Get a parameter of the `beta_dist`
#'
#' @param x a `beta_dist` or `beta_dist_list` acting as the prior
#' @param type the parameter to extract one of `shape1` or `shape2` or `conc`
#'
#' @return a vector of doubles
#' @export
#'
#' @examples
#' get_beta_shape(beta_dist(shape1=1,shape2=1))
#' get_beta_shape(beta_dist(shape1=2:5,shape2=1:4))
get_beta_shape = function(x, type = c("shape1","shape2","conc")) {
  UseMethod("get_beta_shape", x)
}

#' @inherit get_beta_shape
#' @export
get_beta_shape.beta_dist = function(x, type = c("shape1","shape2","conc")) {
  type = match.arg(type)
  return(x[[type]])
}

#' @inherit get_beta_shape
#' @export
get_beta_shape.beta_dist_list = function(x, type = c("shape1","shape2","conc")) {
  type = match.arg(type)
  return(purrr::map_dbl(x, ~ .x[[type]]))
}


#' A uniform prior
#'
#' @return a `beta_dist`
#' @export
uniform_prior = function() {
  beta_dist(shape1 = 1, shape2 = 1)
}

#' Uninformative prior
#'
#' @return a `beta_dist`
#' @export
uninformed_prior = function() {
  beta_dist(shape1 = 0.0001, shape2 = 0.0001)
}

#' The default prior for specificity
#'
#' If undefined this is `r format(beta_dist(p=0.7, n=2), "{sprintf('%1.2f (%1.2f - %1.2f)',mean,lower,upper)}")`. This can be set with `options(testerror.sens_prior = beta_dist(p=??, n=??))`
#'
#' @return a `beta_dist`
#' @export
sens_prior = function() {
  getOption("testerror.sens_prior", beta_dist(p=0.7, n=2))
}

#' The default prior for specificity
#'
#' If undefined this is `r format(beta_dist(p=0.98, n=1), "{sprintf('%1.2f (%1.2f - %1.2f)',mean,lower,upper)}")`. This can be set with `options(testerror.spec_prior = beta_dist(p=??, n=??))`
#'
#' @return a `beta_dist`
#' @export
spec_prior = function() {
  getOption("testerror.spec_prior", beta_dist(p=0.98, n=1))
}