#' Test optimal performance
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
  
  # N.B. this does not use pkgutils check_consistent for a good reason which is
  # they are not related quantities. however the idea of filling in missing
  # values is maybe equivalent
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

