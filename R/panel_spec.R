#' Test panel combination specificity
#' 
#' Calculate the specificity of a combination of tests, where the tests are 
#' testing for different conditions and positive results are combined 
#' into a panel using a logical OR. Because false positives from each component
#' of a panel combine the false positive rate for the panel is higher than the
#' individual components (and hence the true negative rate a.k.a specificity is 
#' lower).
#'
#' @param spec a vector of specificity of the component tests
#' @param na.rm remove NA values?
#'
#' @return a single value for the effective specificity of the combination of 
#' the tests
#' @export
#'
#' @examples
#' panel_spec(spec = rep(0.9975,24))
panel_spec = function(spec, na.rm=FALSE) {
  return(prod(spec, na.rm=na.rm))
}

