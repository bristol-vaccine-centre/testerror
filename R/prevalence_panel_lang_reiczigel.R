#' Lang-Reiczigel true prevalence for panel
#' 
#' Uses resampling to incorporate uncertainty of sensitivity and specificity into
#' an estimate of true prevalence from a given value of apparent prevalence.
#' 
#' This is not vectorised
#' 
#' @param panel_pos_obs the number of positive observations for a given panel of tests
#' @param panel_n_obs a vector of the number of observations for each component test
#' @inheritParams uncertain_rogan_gladen
#' @inheritDotParams prevalence_lang_reiczigel
#' @param samples number of random draws of sensitivity and specificity (optional - default 1000)
#' 
#' @return the expected value of apparent prevalence
#' @export
#'
#' @examples
#' #TODO
prevalence_panel_lang_reiczigel = function(
    panel_pos_obs,
    panel_n_obs,
    pos_obs,
    n_obs,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    spec = spec_prior(),
    sens = sens_prior(),
    confint = 0.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    samples = 1000
) {
  
  # Fit beta distributions to panel sensitivity and specificity estimators:
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
    fit_beta = TRUE
  ) 
  
  panel_spec = uncertain_panel_spec(
    false_pos_controls = false_pos_controls,
    n_controls = n_controls,
    spec = spec, 
    samples = samples,
    fit_beta = TRUE
  )
  
  # Use beta distribution as parameters to normal Lang-Recizigel prevalence calc
  panel = prevalence_lang_reiczigel(
    pos_obs = panel_pos_obs,
    n_obs = panel_n_obs,
    spec = panel_spec,
    sens = panel_sens,
    confint = confint,
    fmt = fmt,
    ...
  )
  
  return(panel)
}