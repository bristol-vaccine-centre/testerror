#' Bayesian simpler model true prevalence for component
#'
#' @param pos_obs a vector of the number of positive observations for each component test
#' @inheritParams uncertain_rogan_gladen
#' @inheritParams rstan::sampling
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
    iter = 2000
) {
  
  n = pkgutils::recycle(pos_obs, n_obs, false_pos_controls, n_controls,
                        false_neg_diseased, n_diseased, sens, spec)
  pkgutils::check_integer(pos_obs, n_obs)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  
  standata = list()
  standata$n = n
  standata$k_sample = array(n_obs, dim=n)
  standata$pos_sample = array(pos_obs, dim=n)
  
  if (is.null(false_pos_controls)) {
    false_pos_controls = rep(0,n)
    n_controls = rep(0,n)
  }
  if (is.null(false_neg_diseased)) {
    false_neg_diseased = rep(0,n)
    n_diseased = rep(0,n)
  }
  
  standata$fp_spec = array(false_pos_controls, dim=n)
  standata$k_spec = array(n_controls, dim=n)
  standata$fn_sens = array(false_neg_diseased, dim=n)
  standata$k_sens = array(n_diseased, dim=n)
  
  standata$tp_sens_prior = array(get_beta_shape(sens, "shape1"), dim=n)
  standata$fn_sens_prior = array(get_beta_shape(sens, "shape2"), dim=n)
  standata$tn_spec_prior = array(get_beta_shape(spec, "shape1"), dim=n)
  standata$fp_spec_prior = array(get_beta_shape(spec, "shape2"), dim=n)
  
  # https://callr.r-lib.org/
  # https://www.jchau.org/2021/02/02/tracking-stan-sampling-progress-shiny/
  # https://github.com/r-lib/progress
  
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





