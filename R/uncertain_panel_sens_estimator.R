#' Propagate component test sensitivity and specificity into panel specificity
#' assuming a known set of observations of component apparent prevalence
#'
#' @inheritParams uncertain_rogan_gladen
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
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    sens = sens_prior(), 
    spec = spec_prior(),
    samples = 1000,
    fit_beta = FALSE
) {
  
  
  n = pkgutils::recycle(
    pos_obs, n_obs, false_pos_controls, n_controls, 
    false_neg_diseased, n_diseased, sens, spec
  )
  
  if (n<2) stop("Panel results require more than one test")
  
  pkgutils::check_integer(pos_obs, n_obs)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  
  sens = update_posterior(sens, neg=false_neg_diseased, n=n_diseased)
  spec = update_posterior(spec, neg=false_pos_controls, n=n_controls)
  
  ap = rep(pos_obs/n_obs,samples)
  ap = matrix(ap,nrow = samples,byrow = TRUE)
  
  # spec=as.beta_dist_list(spec)
  # sens=as.beta_dist_list(sens)
  
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

