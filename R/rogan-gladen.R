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
rogan_gladen = function(ap, sens, spec) { 
  scales::squish((ap-1+spec) / (sens-1+spec)) 
}


