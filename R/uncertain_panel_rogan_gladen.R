#' Rogan-Gladen true prevalence for panel with resampling
#' 
#' Uses resampling to incorporate uncertainty of sensitivity and specificity into
#' an estimate of true prevalence from a given value of apparent prevalence.
#' 
#' This is not vectorised
#' 
#' @param panel_pos_obs the number of positive observations for a given panel of tests
#' @param panel_n_obs the number of observations for each component test
#' @inheritParams uncertain_rogan_gladen
#' @inheritDotParams uncertain_rogan_gladen
#' 
#' @return the expected value of apparent prevalence
#' @export
uncertain_panel_rogan_gladen = function(
    panel_pos_obs,
    panel_n_obs,
    pos_obs,
    n_obs,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    sens = sens_prior(), 
    spec = spec_prior(),
    confint = 0.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    samples = 1000
) { 
  
  pkgutils::check_integer(panel_pos_obs, panel_n_obs, pos_obs, n_obs)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  
  # get sensitivity as set of samples
  panel_sens = uncertain_panel_sens_estimator(
    pos_obs = pos_obs, 
    n_obs = n_obs,
    false_pos_controls = false_pos_controls,
    n_controls = n_controls,
    false_neg_diseased = false_neg_diseased,
    n_diseased = n_diseased,
    sens = sens,
    spec = spec,
    samples = samples,
    fit_beta = FALSE
  )
  
  # get specificity as set of samples
  panel_spec = uncertain_panel_spec(
    false_pos_controls = false_pos_controls,
    n_controls = n_controls,
    spec = spec, 
    samples = samples,
    fit_beta = FALSE
  )
  
  # Use uncertain rogan gladen with sens and spec params being samples
  ap = rep(panel_pos_obs/panel_n_obs,samples)
  p = rogan_gladen(ap, panel_sens, panel_spec)
  
  panel = dplyr::bind_cols(
    .beta_label_3(p, "prevalence", ci = confint, fmt=fmt),
    .beta_label_3(panel_sens, "sens", ci = confint, fmt=fmt),
    .beta_label_3(panel_spec, "spec", ci = confint, fmt=fmt)
  ) %>% mutate(
    prevalence.method = "rogan-gladen (samples)"
  )
  
  return(panel)
}