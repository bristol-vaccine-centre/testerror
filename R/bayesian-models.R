# adaptor code to execute the range of bayesian models included

# model caching ---- 



.model_cache = memoise::cache_filesystem(
  fs::path(rappdirs::user_cache_dir("testerror"),"models"))

#' @importFrom memoise memoise
#' @importFrom digest digest
#' @importFrom rappdirs user_cache_dir
NULL

.build_model = function(model, name, stan_version=rstan::stan_version()) {
  message("compiling stan model: ",name," on first use...")
  rstan::stan_model(model_code = model,model_name = name)
}


.build_model_cached = memoise::memoise(.build_model, cache = .model_cache)
.sampling_cached = memoise::memoise(rstan::sampling)


.get_stan_model = function(name) {
  file = fs::path_ext_set(fs::path(system.file("stan", package="testerror"),name),"stan")
  if (!fs::file_exists(file)) stop("Stan model ",name," is not defined")
  model = readr::read_file(file)
  tmp = .build_model_cached(model, name)
  return(tmp)
}

.get_stan_model_names = function() {
  paths = fs::dir_ls(system.file("stan", package="testerror"),glob = "*.stan")
  paths %>% fs::path_file() %>% fs::path_ext_remove()
}

# specific stan models ----

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
#' @param panel_results 
#'   `r interfacer::idocument(bayesian_panel_complex_model, panel_results)`
#' @param false_pos_controls the number of positives that appeared in the
#'   component tests of a specificity disease-free control group. These are by
#'   definition false positives.
#' @param n_controls the number of controls in the specificity (optional)
#'   disease-free control group (optional).
#' @param false_neg_diseased the number of negatives that appeared in the
#'   component tests of a sensitivity confirmed disease group. These are by
#'   definition false negatives. (optional)
#' @param n_diseased the number of confirmed disease cases in the sensitivity
#'   control group (potional)
#' @param ... not used
#' @param sens the prior sensitivity of the components as a `beta_dist_list`
#'   (optional)
#' @param spec the prior specificity of the components as a `beta_dist_list`
#'   (optional)
#' @param panel_sens the prior sensitivity of the panel as a `beta_dist`
#'   (optional)
#' @param panel_spec the prior specificity of the panel as a `beta_dist`
#'   (optional)
#' @param panel_name the name of the panel for combined result
#' @param confint confidence limit width
#' @param fmt a `sprintf` formatting string accepting 3 numbers
#' @inheritParams rstan::sampling
#'
#' @return a list of dataframes containing the prevalence, sensitivity, and
#'   specificity estimates, and a `stanfit` object with the raw fit data
#' @export
bayesian_panel_complex_model = function(
    test_results = testerror:::.input_data,
    panel_results = testerror:::.input_panel_data,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    sens = NULL, 
    spec = NULL,
    panel_sens = NULL,
    panel_spec = NULL,
    panel_name = "Panel",
    confint=.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    chains = 4,
    warmup = 1000,
    iter = 2000
) {
  test_results = interfacer::ivalidate(test_results, ..., .prune=TRUE)
  panel_results = interfacer::ivalidate(panel_results, ..., .prune=TRUE)
  
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
      spec = sens,
      panel_sens = panel_sens,
      panel_spec = panel_spec
    ))
  
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


#' Bayesian simpler model true prevalence for panel
#'
#' @param panel_pos_obs the number of positive observations for a given panel of tests
#' @param panel_n_obs the number of observations for each component test
#' @param pos_obs a vector of the number of positive observations for each component of a panel test
#' @param n_obs a vector of the number of observations for each component test
#' @param false_pos_controls the number of positives that appeared in the
#'   component tests of a specificity disease-free control group. These are by
#'   definition false positives.
#' @param n_controls the number of controls in the specificity (optional)
#'   disease-free control group (optional).
#' @param false_neg_diseased the number of negatives that appeared in the
#'   component tests of a sensitivity confirmed disease group. These are by
#'   definition false negatives. (optional)
#' @param n_diseased the number of confirmed disease cases in the sensitivity
#'   control group (potional)
#' @param ... not used
#' @param sens the prior sensitivity of the components as a `beta_dist_list`
#'   (optional)
#' @param spec the prior specificity of the components as a `beta_dist_list`
#'   (optional)
#' @param test_names a vector of the component test names in desired order
#'   (optional)
#' @param panel_sens the prior sensitivity of the panel as a `beta_dist`
#'   (optional)
#' @param panel_spec the prior specificity of the panel as a `beta_dist`
#'   (optional)
#' @param panel_name the name of the panel for combined result
#' @param confint confidence limit width
#' @param fmt a `sprintf` formatting string accepting 3 numbers
#' @inheritParams rstan::sampling
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
    sens = NULL, 
    spec = NULL,
    panel_sens = NULL,
    panel_spec = NULL,
    panel_name = "Panel",
    confint=.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    chains = 4,
    warmup = 1000,
    iter = 2000
) {
  
  standata = list()
  
  standata$k_sample = n_obs
  standata$pos_sample = pos_obs
  
  standata$k_sample_combined = panel_n_obs
  standata$pos_sample_combined = panel_pos_obs
  
  n_test = length(standata$pos_sample)
  standata$n_test = n_test
  
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
      spec = sens,
      panel_sens = panel_sens,
      panel_spec = panel_spec
    ))
  
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


#' Bayesian simpler model true prevalence for component
#'
#' @param pos_obs a vector of the number of positive observations for each component test
#' @param n_obs a vector of the number of observations for each component test
#' @param false_pos_controls the number of positives that appeared in the
#'   component tests of a specificity disease-free control group. These are by
#'   definition false positives.
#' @param n_controls the number of controls in the specificity (optional)
#'   disease-free control group (optional).
#' @param false_neg_diseased the number of negatives that appeared in the
#'   component tests of a sensitivity confirmed disease group. These are by
#'   definition false negatives. (optional)
#' @param n_diseased the number of confirmed disease cases in the sensitivity
#'   control group (potional)
#' @param ... not used
#' @param sens the prior sensitivity of the components as a `beta_dist_list`
#'   (optional)
#' @param spec the prior specificity of the components as a `beta_dist_list`
#'   (optional)
#' @param confint confidence limit width
#' @param fmt a `sprintf` formatting string accepting 3 numbers
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
    sens = NULL, 
    spec = NULL,
    confint=.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    chains = 4,
    warmup = 1000,
    iter = 2000
) {
  
  n = .recycle(pos_obs, n_obs, false_pos_controls, n_controls,
               false_neg_diseased, n_diseased, sens, spec)
  
  standata = list()
  standata$n = n
  standata$k_sample = array(n_obs, dim=n)
  standata$pos_sample = array(pos_obs, dim=n)
  
  if (is.null(sens)) sens = beta_dist(p = rep(1,n), q = rep(1,n))
  if (is.null(spec)) spec = beta_dist(p = rep(1,n), q = rep(1,n))
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
  standata$tp_sens_prior = array(as.double(purrr::map_dbl(sens,~ .x$shape1)), dim=n)
  standata$fn_sens_prior = array(as.double(purrr::map_dbl(sens,~ .x$shape2)), dim=n)
  standata$tn_spec_prior = array(as.double(purrr::map_dbl(spec,~ .x$shape1)), dim=n)
  standata$fp_spec_prior = array(as.double(purrr::map_dbl(spec,~ .x$shape2)), dim=n)
  
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




# helper functions ----

.summarise_stan_results = function(fit_combined, test_names, panel_name, confint,  model_name, fmt) {
  
  ci_prob = c(
    median = 0.5,
    lower = (1-confint)/2,
    upper = 1-(1-confint)/2
  )
  
  ci_label = sprintf("%s%%",as.character(ci_prob*100))
  names(ci_label) = names(ci_prob)
  
  summ = rstan::summary(
    fit_combined, 
    pars = c("p","sens","spec"), 
    probs = ci_prob
  )$summary
  
  summ2 = tibble::as_tibble(summ,rownames = "param") %>%
    dplyr::mutate(
      test = test_names[as.integer(stringr::str_extract(param,"[0-9]+"))],
      param = stringr::str_extract(param,"[a-zA-Z]+")
    ) %>% dplyr::rename( !!!ci_label)
  
  if (!is.null(panel_name)) {
    comb_summ = rstan::summary(
      fit_combined, 
      pars = c("p_combined","sens_combined","spec_combined"),
      probs = ci_prob
    )$summary
    comb_summ2 = tibble::as_tibble(comb_summ,rownames = "param") %>%
      dplyr::mutate(
        test = panel_name,
        param = stringr::str_extract(param,"[a-zA-Z]+")
      ) %>% dplyr::rename( !!!ci_label)
  } else {
    comb_summ2 = tibble::tibble()
  }
  
  out = dplyr::bind_rows(summ2,comb_summ2) %>%
    dplyr::select(param,median,lower,upper,test) %>%
    dplyr::mutate(
      label = sprintf(fmt,median*100,lower*100,upper*100),
      method = model_name
    )
  
  result = out %>% dplyr::filter(param=="p") %>% 
    dplyr::select(-param) %>% 
    dplyr::rename_with(.cols = c(-test), .fn = ~ paste0("prevalence.",.x)) %>%
    dplyr::inner_join(
      out %>% dplyr::filter(param=="sens") %>% dplyr::select(-param, -method) %>% 
        dplyr::rename_with(.cols = c(-test), .fn = ~ paste0("sens.",.x)),
      by = "test"
    ) %>% 
    dplyr::inner_join(
      out %>% dplyr::filter(param=="spec") %>% dplyr::select(-param, -method) %>% 
        dplyr::rename_with(.cols = c(-test), .fn = ~ paste0("spec.",.x)),
      by = "test"
    )
  
  return(result)
}

# population optional parts of the model.
.standata_priors = function(
    n_test,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    sens = NULL, 
    spec = NULL,
    panel_sens = NULL,
    panel_spec = NULL,
    ...
  ) {
  
  standata = list()
  
  # Spec priors
  if (!is.null(spec))  {
    # Use spec as beta prior
    standata$tn_spec_prior = as.double(purrr::map_dbl(spec,~ .x$shape1))
    standata$fp_spec_prior = as.double(purrr::map_dbl(spec,~ .x$shape2))
  } else {
    standata$tn_spec_prior = rep(1,n_test)
    standata$fp_spec_prior = rep(1,n_test)
  }
  
  # Sens priors
  if (!is.null(sens)) {
    standata$tp_sens_prior = as.double(purrr::map_dbl(sens,~ .x$shape1))
    standata$fn_sens_prior = as.double(purrr::map_dbl(sens,~ .x$shape2))
  } else {
    standata$tp_sens_prior = rep(1,n_test)
    standata$fn_sens_prior = rep(1,n_test)
  }
  
  # negative controls (spec)
  if (!is.null(false_pos_controls)) {
    # Use fp_spec as binomial outcome
    standata$k_spec = as.integer(n_controls)
    standata$fp_spec = as.integer(false_pos_controls)
  } else {
    standata$k_spec = rep(0,n_test)
    standata$fp_spec = rep(0,n_test)
  }
  
  # positive controls (sens)  
  if (!is.null(false_neg_diseased)) {
    standata$k_sens = as.integer(n_diseased)
    standata$fn_sens = as.integer(false_neg_diseased)
  } else {
    standata$k_sens = rep(0,n_test)
    standata$fp_sens = rep(0,n_test)
  }
  
  # Panel sens priors
  if (!is.null(panel_sens)) {
    standata$tp_panel_sens_prior = as.double(panel_sens$shape1)
    standata$fn_panel_sens_prior = as.double(panel_sens$shape2)
  } else {
    standata$tp_panel_sens_prior = 1
    standata$fn_panel_sens_prior = 1
  }
  
  # Panel spec priors
  if (!is.null(panel_spec)) {
    standata$tn_panel_spec_prior = as.double(panel_spec$shape1)
    standata$fp_panel_spec_prior = as.double(panel_spec$shape2)
  } else {
    standata$tn_panel_spec_prior = 1
    standata$fp_panel_spec_prior = 1
  } 
  
  return(standata)
}
