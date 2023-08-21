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
apparent_prevalence = function(p, sens, spec) { 
  p*sens+(1-p)*(1-spec) 
}
