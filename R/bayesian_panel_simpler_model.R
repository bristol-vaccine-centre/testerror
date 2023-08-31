#' Bayesian simpler model true prevalence for panel
#'
#' @param panel_pos_obs the number of positive observations for a given panel of tests
#' @param panel_n_obs the number of observations for each component test
#' @inheritParams uncertain_rogan_gladen
#' @inheritParams rstan::sampling
#' @param test_names a vector of the component test names in desired order
#' @param panel_sens the prior sensitivity of the panel as a `beta_dist`
#'   (optional)
#' @param panel_spec the prior specificity of the panel as a `beta_dist`
#'   (optional)
#' @param panel_name the name of the panel for combined result
#' @param cache_result save the result of the sampling in memory for the current 
#'   session
#'
#' @return a list of dataframes containing the prevalence, sensitivity, and
#'   specificity estimates, and a `stanfit` object with the raw fit data
#' @export
bayesian_panel_simpler_model = function(
    panel_pos_obs,
    panel_n_obs,
    pos_obs,
    n_obs,
    test_names,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    sens = uniform_prior(), 
    spec = uniform_prior(),
    panel_sens = uniform_prior(),
    panel_spec = uniform_prior(),
    panel_name = "Panel",
    confint=.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    chains = 4,
    warmup = 1000,
    iter = 2000,
    cache_result = TRUE
) {
  
  if (is.null(sens)) sens = uniform_prior()
  if (is.null(spec)) spec = uniform_prior()
  
  n = pkgutils::recycle(pos_obs, n_obs, false_pos_controls, n_controls,
                        false_neg_diseased, n_diseased, sens, spec)
  pkgutils::check_integer(panel_pos_obs, panel_n_obs, pos_obs, n_obs)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  
  standata = list()
  
  standata$n_test = n
  standata$k_sample = n_obs
  standata$pos_sample = pos_obs
  
  standata$k_sample_combined = panel_n_obs
  standata$pos_sample_combined = panel_pos_obs
  
  #n_test = length(standata$pos_sample)
  #standata$n_test = n
  
  # Beta(1,1) = U(0,1)
  
  standata = c(
    standata,
    .standata_priors(
      n_test = n,
      false_pos_controls = false_pos_controls,
      n_controls = n_controls,
      false_neg_diseased = false_neg_diseased,
      n_diseased = n_diseased,
      sens = sens, 
      spec = spec,
      panel_sens = panel_sens,
      panel_spec = panel_spec
    ))
  
  if (!cache_result) {
    memoise::drop_cache(.sampling_cached)(
      .get_stan_model("panel-simpler"),
      data = standata,
      chains = chains,
      warmup = warmup,
      iter = iter,
      show_messages = getOption("testerror.debug",FALSE)
    )
  }
  
  fit_combined = .sampling_cached(
    .get_stan_model("panel-simpler"),
    data = standata,
    chains = chains,
    warmup = warmup,
    iter = iter,
    show_messages = getOption("testerror.debug",FALSE)
  )
  
  result = .summarise_stan_results(
    fit_combined, test_names, panel_name, confint, "bayes (binom)", fmt)
  
  return(list(
    summary = result,
    fit = fit_combined
  ))
}


