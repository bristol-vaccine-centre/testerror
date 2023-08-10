

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

#' Propagate component test specificity into panel specificity
#'
#' @inheritParams uncertain_panel_rogan_gladen
#' @param fit_beta return the result as a `beta_dist` object?
#' @param na.rm remove missing values
#'
#' @return a vector of possible specificities for the panel or a fitted `beta_dist`
#' @export
#'
#' @examples
#' uncertain_panel_spec(c(2,3,4,2,2), c(800,800,800,800,800), fit_beta=TRUE)
uncertain_panel_spec = function(
    false_pos_controls = NA, n_controls = NA, 
    ..., 
    spec = beta_dist(n_controls-false_pos_controls+1, false_pos_controls+1), 
    samples=1000, 
    na.rm=FALSE,
    fit_beta = FALSE) {
  
  if (inherits(spec,"beta_dist")) spec=list(spec)
  mat = sapply(1:length(spec), function(i) spec[[i]]$r(samples))
  
  specs = apply(mat, MARGIN = 1, prod, na.rm=na.rm)
  
  specs[specs < 0] = 0
  specs[specs > 1] = 1
  
  if (fit_beta) return(beta_fit(specs))
  return(specs)
}


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


#' Propagate component test sensitivity and specificity into panel specificity
#' assuming a known set of observations of component apparent prevalence
#'
#' @inheritParams uncertain_panel_rogan_gladen
#' @param fit_beta return the result as a `beta_dist` object?
#' 
#' 
#' @return a vector of possible sensitivity values
#'
#' @export
#' @examples 
#' uncertain_panel_sens_estimator(
#'   pos_obs = c(30,10,20,10,5), n_obs=1000,
#'   false_pos_controls = c(20,15,15,15,15), n_controls = c(800,800,800,800,800),
#'   false_neg_diseased = c(20,25,20,20,15), n_diseased = c(100,100,100,100,100),
#'   fit_beta = TRUE)
uncertain_panel_sens_estimator = function(
    pos_obs,
    n_obs,
    false_pos_controls = NA,
    n_controls = NA,
    false_neg_diseased = NA,
    n_diseased = NA,
    ...,
    sens = beta_dist(n_diseased-false_neg_diseased+1, false_neg_diseased+1), 
    spec = beta_dist(n_controls-false_pos_controls+1, false_pos_controls+1),
    samples = 1000,
    fit_beta = FALSE
  ) {
  
  if (length(n_obs) == 1) n_obs = rep(n_obs,length(pos_obs))
  if (length(pos_obs) != length(sens) || length(pos_obs) != length(sens) || length(pos_obs) != length(n_obs)) stop("all provided parameters should be the same length") 
  ap = rep(pos_obs/n_obs,samples)
  ap = matrix(ap,nrow = samples,byrow = TRUE)
  
  if (inherits(spec,"beta_dist")) spec=list(spec)
  if (inherits(sens,"beta_dist")) sens=list(sens)
  
  spec_mat = sapply(1:length(spec), function(i) spec[[i]]$r(samples))
  sens_mat = sapply(1:length(sens), function(i) sens[[i]]$r(samples))
  
  one_minus_prev = (sens_mat-ap)/(spec_mat+sens_mat-1)
  # one_minus_prev[ap < 1-spec_mat] = 1
  one_minus_prev[sens_mat < ap] = 0
  
  sens_N = 1-(
    apply(1-ap, MARGIN = 1, prod) - apply(spec_mat*one_minus_prev, MARGIN = 1, prod)
  )/(
    1 - apply( one_minus_prev, MARGIN=1, prod )
  )
  
  sens_N[sens_N < 0] = 0
  sens_N[sens_N > 1] = 1
  
  if (fit_beta) return(beta_fit(sens_N))
  
  return(sens_N)
}

#' Estimate test panel combination sensitivity
#' 
#' Estimate the sensitivity of a combination of tests, where the tests are 
#' testing for different conditions and positive results are combined 
#' into a panel using a logical OR. Because false negatives from each component
#' of a panel can be cancelled out by true positives, or false positives from 
#' other components of the test depending on the prevalence of the underlying 
#' conditions, the combined false negative rate is lower the more cases there 
#' arecombine the false positive rate for the panel is higher than the
#' individual components (and hence the true negative rate a.k.a specificity is 
#' lower).
#'
#' @param ap the apparent prevalence or test positivity (one of p or ap must be given)
#' @param sens a vector of sensitivities of the component tests
#' @param spec a vector of specificity of the component tests
#' @param na.rm remove NA values?
#'
#' @return an effective specificity for the combination of the tests
#' @export
#'
#' @examples
#' #TODO
panel_sens_estimator = function(ap, sens, spec, na.rm=FALSE) {
  
  if (length(sens) == 1) sens = rep(sens,length(ap))
  if (length(spec) == 1) spec = rep(spec,length(ap))
  
  
  if (na.rm) {
    ap = ap[!(is.na(ap) | is.na(sens) | is.na(spec))]
    sens = sens[!(is.na(ap) | is.na(sens) | is.na(spec))]
    spec = spec[!(is.na(ap) | is.na(sens) | is.na(spec))]
  }
  
  # this seems to cause an issue
  if (all(abs(ap - (1-spec)) <= .Machine$double.eps)) {
    return(NA_real_)
  }
  
  one_minus_prev = dplyr::case_when(
    # ap < 1-spec ~ 1,
    sens < ap ~ 0,
    TRUE ~ (sens-ap)/(spec+sens-1)
  )
  
  sens_N = 1-(
      prod(1-ap) - prod(spec*one_minus_prev)
    )/(
      1 - prod( one_minus_prev )
    )
    
  return(sens_N)
}


