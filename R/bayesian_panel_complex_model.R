#' Bayesian models true prevalence for panel
#'
#' Uses resampling to incorporate uncertainty of sensitivity and specificity
#' into an estimate of true prevalence from a given value of apparent
#' prevalence.
#'
#' This is not vectorised
#'
#' @param test_results 
#'   `r interfacer::idocument(bayesian_panel_complex_model, test_results)`
#' @inheritParams uncertain_rogan_gladen
#' @inheritParams rstan::sampling
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
bayesian_panel_complex_model = function(
    test_results = testerror::.input_data,
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
  test_results = interfacer::ivalidate(test_results, ..., .prune=TRUE)
  
  if (is.null(sens)) sens = uniform_prior()
  if (is.null(spec)) spec = uniform_prior()
  
  n = pkgutils::recycle(false_pos_controls, n_controls,
                        false_neg_diseased, n_diseased, sens, spec)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  

  
  # Summarise by test to get a per patient panel result propagating NAs
  panel_results = test_results %>% 
    dplyr::group_by(id) %>% 
    dplyr::summarise(result = any(result))
  
  test_names = levels(test_results$test)
  
  standata = list()
  
  # The component tests
  # make input matrix for models that need line by line input
  standata$y_sample = test_results %>%
    dplyr::mutate(test = as.character(test), result = as.integer(result)) %>%
    tidyr::pivot_wider(names_from = test, values_from = result, values_fill = NA_integer_) %>%
    dplyr::select(-id) %>% 
    as.matrix()
  # The combined result
  standata$y_sample_combined = as.integer(panel_results$result)
  
  n_test = ncol(standata$y_sample)
  
  if (n_test != n) stop(
    "There are a different number of tests in the `test_results` data frame 
    compared to the sensitivity and specificity parameters (false_pos_controls 
    etc...)"
  )
  
  standata$n_test = n_test
  standata$k_sample = nrow(standata$y_sample)
  
  # Beta(1,1) = U(0,1)
  
  standata = c(
    standata,
    .standata_priors(
      n_test = n_test,
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
    .get_stan_model("panel-complex"),
    data = standata,
    chains = chains,
    warmup = warmup,
    iter = iter,
    show_messages = getOption("testerror.debug",FALSE)
  )
  
  comp = test_results %>% 
    dplyr::summarise(pos_obs = sum(result), n_obs=dplyr::n(), .by="test")
  pnl = panel_results %>% 
    dplyr::summarise(pos_obs = sum(result), n_obs=dplyr::n())
  
  result = .summarise_stan_results(
    fit_combined, test_names, panel_name, confint, "bayes (bernoulli)", fmt)
  
  return(list(
    summary = result,
    fit = fit_combined
  ))
}


