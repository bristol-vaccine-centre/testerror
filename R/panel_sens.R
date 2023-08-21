#' Test panel combination sensitivity
#' 
#' Calculate the sensitivity of a combination of tests, where the tests are 
#' testing for different conditions and positive results are combined 
#' into a panel using a logical OR. Because false negatives from each component
#' of a panel can be cancelled out by true positives, or false positives from 
#' other components of the test depending on the prevalence of the underlying 
#' conditions, the combined false negative rate is lower the more cases there 
#' arecombine the false positive rate for the panel is higher than the
#' individual components (and hence the true negative rate a.k.a specificity is 
#' lower).
#'
#' @param p the true prevalence (one of p or ap must be given)
#' @param sens a vector of sensitivities of the component tests
#' @param spec a vector of specificity of the component tests
#' @param na.rm remove NA values?
#'
#' @return an effective specificity for the combination of the tests
#' @export
#'
#' @examples
#' #TODO
panel_sens = function(p, sens, spec, na.rm=FALSE) {
  
  if (length(sens) == 1) sens = rep(sens,length(p))
  if (length(spec) == 1) spec = rep(spec,length(p))
  
  if (na.rm) {
    p = p[!(is.na(p) | is.na(sens) | is.na(spec))]
    sens = sens[!(is.na(ap) | is.na(sens) | is.na(spec))]
    spec = spec[!(is.na(ap) | is.na(sens) | is.na(spec))]
  }
  
  sens_N = 1-(
    prod((1-sens) * p + spec * (1-p) ) - prod(spec * (1-p))
  )/(
    1 - prod(1-p)
  )
  
  return(sens_N)
}


