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
bayesian_component_logit_model = function(
    pos_obs,
    n_obs,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    sens = sens_prior(), 
    spec = spec_prior(),
    confint=.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    chains = 4,
    warmup = 1000,
    iter = 2000,
    cache_result = TRUE
) {
  
  if (is.null(sens)) sens = sens_prior()
  if (is.null(spec)) spec = spec_prior()
  
  n = pkgutils::recycle(pos_obs, n_obs, false_pos_controls, n_controls,
                        false_neg_diseased, n_diseased, sens, spec)
  pkgutils::check_integer(pos_obs, n_obs)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  
  standata = list()
  standata$n = n
  standata$k_sample = array(n_obs, dim=n)
  standata$pos_sample = array(pos_obs, dim=n)
  
  # create the priors for prevalance
  sens_point_estimate = sens %>% update_posterior(neg = false_neg_diseased, n = n_diseased) %>% as_tibble() %>% dplyr::pull(mean)
  spec_point_estimate = spec %>% update_posterior(neg = false_pos_controls, n = n_controls) %>% as_tibble() %>% dplyr::pull(mean)
  ep = rogan_gladen(pos_obs/n_obs, sens = sens_point_estimate, spec=spec_point_estimate) %>% scales::squish(range=c(0.001, 0.999))
  
  standata$logit_mu_p_prior = array(logit(ep), dim=n)
  standata$logit_sigma_p_prior = array(abs(logit(ep))/2, dim=n)
  
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
      logit = TRUE
    ))
  
  # https://callr.r-lib.org/
  # https://www.jchau.org/2021/02/02/tracking-stan-sampling-progress-shiny/
  # https://github.com/r-lib/progress
  
  if (!cache_result) {
    memoise::drop_cache(.sampling_cached)(
      .get_stan_model("component-logit"),
      data = standata,
      chains = chains,
      warmup = warmup,
      iter = iter,
      show_messages = getOption("testerror.debug",FALSE)
    )
  }
  
  fit_combined = .sampling_cached(
    .get_stan_model("component-logit"),
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
    model_name = "bayes (logit)", 
    fmt = fmt) %>% dplyr::arrange(test) %>% dplyr::select(-test)
  
  return(list(
    summary = result,
    fit = fit_combined
  ))
}





