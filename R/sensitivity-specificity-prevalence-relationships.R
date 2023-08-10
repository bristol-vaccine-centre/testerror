
false_pos = function(n, spec) n*(1-spec)
true_neg = function(n, spec) n*spec
false_neg = function(n, sens) n*(1-sens)
true_pos = function(n, sens) n*sens

#' True prevalence from apparent prevalence
#' 
#' This estimator runs into problems with small AP as the Rogan-Gladen
#' conversion is really using expected apparent prevalence. Getting the expected
#' value of the AP distribution is complex and the expected value given a single
#' observation is not in general the ratio of positives / count. The expected
#' apparent prevalence is never less than the specificity but the observed value
#' often is. To deal with this the R-G estimator truncates at zero.
#'
#' @param ap the expected apparent prevalence.
#' @param sens the sensitivity of the test
#' @param spec the specificity of the test
#'
#' @return the estimate of 'true prevalence'
#' @export
#'
#' @examples
#' rogan_gladen(50/200, 0.75, 0.97)
rogan_gladen = function(ap, sens, spec) { scales::squish((ap-1+spec) / (sens-1+spec)) }
# est.prev = (obs.prev+obs.spec-1)/(obs.sens+obs.spec-1)


#' Test underestimation limit
#'
#' For a given sensitivity and specificity this give the critical threshold after which 
#' test error introduces underestimation rather than over estimation
#'
#' @param sens the sensitivity of the test
#' @param spec the specificity of the test
#'
#' @return the value where apparent prevalence equals true prevalence
#' @export
#'
#' @examples
#' tmp1 = underestimation_threshold(0.75, 0.97)
#' tmp2 = rogan_gladen(tmp1, 0.75, 0.97)
#' if (abs(tmp1-tmp2) > 0.0000001) stop("error")
#' 
underestimation_threshold = function(sens, spec) {
  return((1-spec) / (1-sens+1-spec))
}

#' Test underestimation critical limit
#'
#' For a given combination of prevalence, sensitivity and specificity this gives 
#' the critical threshold at which true prevalence equals apparent prevalence
#'
#' @param p the prevalence or apparent prevalence
#' @param sens the sensitivity of the test
#' @param spec the specificity of the test
#'
#' @return the combination of sensitivity and specificity where apparent prevalence equals true prevalence
#' @export
#'
#' @examples
#' optimal_performance(p=0.1, sens=0.75)
#' optimal_performance(p=0.005, spec=0.9975)
optimal_performance = function(p=NULL, sens=NULL, spec=NULL) {
  if (is.null(p)+is.null(sens)+is.null(spec) != 1) stop("exactly two of p, sens and spec must be given")
  if (length(p) > 1 || length(sens) > 1 || length(spec) > 1) stop("this is not vectorised. please use sapply.")
  if (is.null(sens)) {
    f = function(s) p-underestimation_threshold(s, spec)
    sens = stats::uniroot(f, interval = c(0,1))$root
  } else if (is.null(spec)) {
    f = function(s) p-underestimation_threshold(sens, s)
    spec = stats::uniroot(f, interval = c(0,1))$root
  } else {
    p = underestimation_threshold(sens, spec)
  }
  return(list(p=p,sens=sens,spec=spec))
}

#' Apparent prevalence from known prevalence
#' 
#' The observed counts of disease is going to be a binomial but 
#' with the apparent prevalence as a probability. This will never be less than 
#' `(1-specificity)` of the test (and never more than the sensitivity). When either
#' of those quantities are uncertain the shape of the distribution of observed counts 
#' is not clear cut.
#'
#' @param p the true value of the prevalence
#' @param sens the sensitivity of the test
#' @param spec the specificity of the test
#'
#' @return the expected value of apparent prevalence
#' @export
#'
#' @examples
#' apparent_prevalence(0, 0.75, 0.97)
#' apparent_prevalence(1, 0.75, 0.97)
apparent_prevalence = function(p, sens, spec) { p*sens+(1-p)*(1-spec) }
