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
  # tmp = try(suppressWarnings(
  #   MASS::fitdistr(x = samples,densfun = "beta", start = list(shape1 = mean(samples),shape2 = 1-mean(samples)))
  # ),silent = TRUE)
  # if (isa(tmp,"try-error")) {
    v = stats::sd(samples)^2
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
#' @param ... not used
#'
#' @return either a single `beta_dist` object or a list of `beta_dist`s
#'
#' @export
#' @examples 
#' beta_dist(c(1,2,3),c(3,2,1))
beta_dist = function(p, q=NULL, n=NULL, ...) {
  # TODO: change this to be explicit about shape versus probability
  .recycle(p,q,n)
  if (is.null(q) & is.null(n)) stop("one of shape2 (q) or concentration (n) must be given")
  if (!is.null(q) & !is.null(n)) {
    if (all(p+q == 1)) {
      shape1 = p*n
      shape2 = q*n
    } else {
      if (!all(p+q == n)) stop("shape parameters (p, q) must add up to 1 or n")
      shape1 = p
      shape2 = q
    }
  } else if (!is.null(q)) {
    shape1 = p
    shape2 = q
  } else if (!is.null(n)) {
    if (all(p <= 1)) {
      shape1 = p*n
      shape2 = (1-p)*n
    } else {
      shape1 = p
      shape2 = n-p
    }
  }
  
  if (length(shape1) > 1) {
    tmp = lapply(seq_along(p), function(i) beta_dist(p=shape1[[i]], q=shape2[[i]], ...))
    return(structure(tmp, class=c("beta_dist_list",class(tmp))))
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

#' Update the posterior of a `beta_dist`
#'
#' @param x a `beta_dist` or `beta_dist_list` acting as the prior
#' @param pos positive observation(s)
#' @param n number observations
#' @param ... not used
#'
#' @return a new `beta_dist` o `beta_dist_list`
#' @export
#'
#' @examples
#' update_posterior(beta_dist(1,1), 10, 30)
update_posterior = function(x, pos, n, ... ) {
  UseMethod("update_posterior", x)
}

#' @inherit update_posterior
#' @export
update_posterior.beta_dist = function(x, pos, n, ...) {
  if (n<pos) stop("update_posterior called with more positives (pos) than observations (n)")
  beta_dist(p=x$shape1+pos, q=x$shape2+n-pos)
}

#' @inherit update_posterior
#' @export
update_posterior.beta_dist_list = function(x, pos, n, ...) {
  tmp = lapply(seq_along(x), function(i) {
    update_posterior(x[[i]],pos[[i]],n[[i]])
  })
  return(structure(tmp, class=c("beta_dist_list",class(tmp))))
}

rep.beta_dist = function(x, times) {
  tmp = rep(list(x), times)
  return(structure(tmp, class=c("beta_dist_list",class(tmp))))
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
#' format(beta_dist(3,6), "{format(mean*100, digits=3)}%")
format.beta_dist = function(x, glue = "{sprintf('%1.3f [%1.3f\u2013%1.3f] (N=%1.2f)',median,lower,upper,conc)}", ...) {
  tmp = as_tibble.beta_dist(x)
  glue::glue_data(tmp, glue)
  # suppressWarnings(do.call(sprintf, c(list(fmt = fmt), x$shape1/ (x$shape1 + x$shape2), 
  #                    stats::qbeta(c(0.025,0.975), shape1 = x$shape1, shape2 = x$shape2),x$conc
  # )))
}

#' Format a beta distribution list
#' 
#' @param x the beta distribution list
#' @param ... not used
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

#' Print a beta distribution
#' 
#' @param x the beta distribution
#' @param ... not used
#' @return nothing
#' @export
print.beta_dist = function(x, ...) {
  cat(format(x, ...), "\n")
}
