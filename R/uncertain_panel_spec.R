#' Propagate component test specificity into panel specificity
#'
#' @inheritParams uncertain_rogan_gladen
#' @param fit_beta return the result as a `beta_dist` object?
#' @param na.rm remove missing values
#'
#' @return a vector of possible specificities for the panel or a fitted `beta_dist`
#' @export
#'
#' @examples
#' uncertain_panel_spec(c(2,3,4,2,2), c(800,800,800,800,800), fit_beta=TRUE)
uncertain_panel_spec = function(
    false_pos_controls = NULL, 
    n_controls = NULL, 
    ..., 
    spec = spec_prior(), 
    samples=1000, 
    na.rm=FALSE,
    fit_beta = FALSE) {
  
  pkgutils::recycle(false_pos_controls,n_controls,spec)
  pkgutils::check_integer(false_pos_controls, n_controls)
  
  spec = update_posterior(spec, neg = false_pos_controls, n = n_controls)
  
  if (inherits(spec,"beta_dist")) spec=list(spec)
  mat = sapply(1:length(spec), function(i) spec[[i]]$r(samples))
  
  specs = apply(mat, MARGIN = 1, prod, na.rm=na.rm)
  
  specs[specs < 0] = 0
  specs[specs > 1] = 1
  
  if (fit_beta) return(beta_fit(specs))
  return(specs)
}

