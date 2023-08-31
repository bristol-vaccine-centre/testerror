#' Bayesian simpler model true prevalence for component
#'
#' @inheritParams uncertain_rogan_gladen
#' @inheritParams rstan::sampling
#' @param cache_result save the result of the sampling in memory for the current 
#'   session
#'
#' @return a list of dataframes containing the prevalence, sensitivity, and
#'   specificity estimates, and a `stanfit` object with the raw fit data
#' @export
bayesian_component_simpler_model = function(
    pos_obs,
    n_obs,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    sens = uniform_prior(), 
    spec = uniform_prior(),
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
  pkgutils::check_integer(pos_obs, n_obs)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  
  standata = list()
  standata$n = n
  standata$k_sample = array(n_obs, dim=n)
  standata$pos_sample = array(pos_obs, dim=n)
  
  standata = c(
    standata,
    .standata_priors(
      n_test = n,
      false_pos_controls = false_pos_controls,
      n_controls = n_controls,
      false_neg_diseased = false_neg_diseased,
      n_diseased = n_diseased,
      sens = sens, 
      spec = spec
    ))
  
  # https://callr.r-lib.org/
  # https://www.jchau.org/2021/02/02/tracking-stan-sampling-progress-shiny/
  # https://github.com/r-lib/progress
  
  if (!cache_result) {
    memoise::drop_cache(.sampling_cached)(
      .get_stan_model("component-simpler"),
      data = standata,
      chains = chains,
      warmup = warmup,
      iter = iter,
      show_messages = getOption("testerror.debug",FALSE)
    )
  }
  
  fit_combined = .sampling_cached(
    .get_stan_model("component-simpler"),
    data = standata,
    chains = chains,
    warmup = warmup,
    iter = iter
    # refresh = -1,
    # show_messages = getOption("testerror.debug",FALSE)
  )
  
  result = .summarise_stan_results(
    fit_combined, 
    test_names = 1:n, 
    panel_name = NULL, 
    confint = confint, 
    model_name = "bayes (binom)", 
    fmt = fmt) %>% dplyr::arrange(test) %>% dplyr::select(-test)
  
  return(list(
    summary = result,
    fit = fit_combined
  ))
}





