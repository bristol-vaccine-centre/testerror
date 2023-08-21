#' Expected test panel prevalence assuming independence
#' 
#' @param p a vector of prevalences of the component tests
#' @param na.rm remove NA values?
#'
#' @return a single value for the effective specificity of the combination of 
#' the tests
#' @export
#'
#' @examples
#' panel_prevalence(p = rep(0.01,24))
panel_prevalence = function(p, na.rm=FALSE) {
  return(1-prod(1-p, na.rm=na.rm))
}


