#' True prevalence from apparent prevalence with uncertainty
#' 
#' Uses resampling to incorporate uncertainty of sensitivity and specificity into
#' an estimate of true prevalence from a given value of apparent prevalence.
#' 
#' @param pos_obs the number of positive observations for a given test
#' @param n_obs the number of observations for a given test
#' @param false_pos_controls the number of positives that appeared in the specificity
#'   disease-free control group. These are by definition false positives. This
#'   is (1-specificity)*n_controls
#' @param n_controls the number of controls in the specificity
#'   disease-free control group. 
#' @param false_neg_diseased the number of negatives that appeared in the sensitivity
#'   confirmed disease group. These are by definition false negatives. This
#'   is (1-sensitivity)*n_diseased
#' @param n_diseased the number of confirmed disease cases in the sensitivity
#'   control group.
#' @param ... not used
#' @param sens the prior sensitivity of the test as a `beta_dist`.
#' @param spec the prior specificity of the test as a `beta_dist`.
#' @param samples number fo random draws of sensitivity and specificity
#' @param confint confidence interval limits
#' @param fmt a `sprintf` formatting string accepting 3 numbers
#' @param seed set seed for reproducibility
#' @param prefix column name prefix
#' 
#' @return the expected value of apparent prevalence
#' @export
#'
#' @examples
#' uncertain_rogan_gladen(
#'   pos_obs = 20, n_obs = 1000, 
#'   false_pos_controls = 10, n_controls = 800, 
#'   false_neg_diseased = 20, n_diseased = 100)
#'   
#' uncertain_rogan_gladen(
#'   pos_obs = 5, n_obs = 1000, 
#'   sens = beta_dist(0.75,n=200), 
#'   spec = beta_dist(0.9975, n=800))
#' 
#' uncertain_rogan_gladen(
#'   pos_obs = c(5,10), n_obs = c(1000,1000), 
#'   false_pos_controls = c(2,1), n_controls = c(800,800), 
#'   false_neg_diseased = c(25,20),n_diseased = c(100,100))
#' 
uncertain_rogan_gladen = function(
    pos_obs,
    n_obs,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    # posterior beta distribution assuming a uniform prior 
    spec = spec_prior(),
    sens = sens_prior(),
    samples = 1000,
    confint = 0.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    seed = NA
) { 
  
  # This can handle vector inputs (but isn't vectorised per se)
  # N.B. this documentation is the core parameter documentation for all other functions
  
  if (!is.na(seed)) set.seed(seed)
  
  n = pkgutils::recycle(
    pos_obs, n_obs, false_pos_controls, n_controls, 
    false_neg_diseased, n_diseased, sens, spec
  )
  
  pkgutils::check_integer(pos_obs, n_obs)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  
  
  if (n>1) {
    # The sens and spec parameters will be a beta_dist_list
    return(dplyr::bind_rows(lapply(1:n, function(i) {
      return(uncertain_rogan_gladen(
        pos_obs = pos_obs[[i]], 
        n_obs = n_obs[[i]], 
        false_pos_controls = false_pos_controls[[i]],
        n_controls = n_controls[[i]],
        false_neg_diseased = false_neg_diseased[[i]],
        n_diseased = n_diseased[[i]],
        sens = sens[[i]], 
        spec = spec[[i]], 
        samples=samples, confint = confint, fmt = fmt, seed=seed, prefix = prefix))
    })))
  }
  
  # At this stage the inputs are of length 1.
  sens = update_posterior(sens, neg=false_neg_diseased, n=n_diseased)
  spec = update_posterior(spec, neg=false_pos_controls, n=n_controls)
  
  # Convert beta_dist inputs into a vector of samples
  rsens = sens$r(samples)
  rspec = spec$r(samples)
  
  # ap = stats::rbeta(samples, pos_obs, n_obs-pos_obs)
  ap = rep(pos_obs/n_obs,samples)
  p = rogan_gladen(ap,rsens,rspec)
  
  out = dplyr::bind_cols(
    .beta_label_3(p, "prevalence", ci = confint, fmt=fmt),
    .beta_label(spec, "spec", ci = confint, fmt=fmt),
    .beta_label(sens, "sens", ci = confint, fmt=fmt),
  ) %>% 
  dplyr::mutate(
    prevalence.method = "rogan-gladen (samples)"
  )
  
  return(out)
}
