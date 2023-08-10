.input_data = interfacer::iface(
  id = character ~ "the patient identifier",
  test = factor ~ "the test type",
  result = logical ~ "the test result",
  .groups = FALSE
)

.input_panel_data = interfacer::iface(
  id = character ~ "the patient identifier",
  result = logical ~ "the panel result",
  .groups = FALSE
)

.output_data = interfacer::iface(
  test = character ~ "the name of the test or panel",
  prevalence.lower = numeric ~ "the lower estimate",
  prevalence.median = numeric ~ "the median estimate",
  prevalence.upper = numeric ~ "the upper estimate",
  prevalence.method = character ~ "the method of estimation",
  prevalence.label = character ~ "a fomatted label of the true prevalence estimate with CI",
  .groups = FALSE
)

#' Calculate an estimate of true prevalence for a single panel and components
#' 
#' Uses apparent prevalence, and uncertain
#' estimates of test sensitivity and test specificity for the 3 methods
#' described in Supplementary 2. This function works for a single panel per
#' dataframe, multiple panels will need to call this function multiple times
#' in a `group_modify`. 
#'
#' @param test_results `r interfacer::idocument(true_panel_prevalence, test_results)`
#' @param false_pos_controls the number of positives that appeared in the specificity
#'   disease-free control group. These are by definition false positives. This
#'   is (1-specificity)*n_controls (optional if `spec` is given)
#' @param n_controls the number of controls in the specificity
#'   disease-free control group. (optional if `spec` is given) 
#' @param false_neg_diseased the number of negatives that appeared in the sensitivity
#'   confirmed disease group. These are by definition false negatives. This
#'   is (1-sensitivity)*n_diseased  (optional if `sens` is given)
#' @param n_diseased the number of confirmed disease cases in the sensitivity
#'   control group.  (optional if `sens` is given)
#' @param ... passed onto methods
#' @param sens the sensitivity of the component tests of test as a `beta_dist_list`. 
#'   (optional if `false_pos_controls` and `n_controls` is given)
#' @param spec the specificity of the test as a `beta_dist_list`
#'   (optional if `false_neg_diseased` and `n_diseased` is given)
#' @param panel_name the name of the panel for combined result
#' @param confint confidence limit width
#' @param method one of:
#' * "lang-reiczigel": Frequentist confidence limits
#' * "bayes": Bayesian analysis
#' * "rogan-gladen": Rogan gladen with uncertainty
#' @param na.rm exclude patients with missing results
#'
#' @return `r format(testerror:::.output_data)`
#' @export
#'
#' @examples
#' tmp = testerror:::panel_example()
#' true_panel_prevalence(
#'   test_results = tmp$samples %>% dplyr::select(id,test,result = observed),
#'   false_pos_controls = tmp$performance$false_pos_controls,
#'   n_controls = tmp$performance$n_controls,
#'   false_neg_diseased = tmp$performance$false_neg_diseased,
#'   n_diseased = tmp$performance$n_diseased,
#'   method = "rogan-gladen"
#' )
true_panel_prevalence = function(
    test_results = testerror:::.input_data,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    sens = NULL, 
    spec = NULL,
    panel_name = "Panel",
    confint=.95,
    method = c("rogan-gladen", "lang-reiczigel", "bayes"),
    na.rm = TRUE
) {
  
  test_results = interfacer::ivalidate(test_results, ..., .prune=TRUE)
  test_names = levels(test_results$test)
  n_test = length(test_names)
  if (inherits(sens, "beta_dist")) sens = list(sens) 
  if (inherits(spec, "beta_dist")) spec = list(spec)
  
  
  # Fill in missing test results
  test_results = test_results %>% 
    tidyr::complete(test, fill = list(result=NA))
  
  if (any(is.na(test_results$result))) {
    if (na.rm) {
      message("Excluding patients with some missing results.")
      test_results = test_results %>% 
        dplyr::group_by(id) %>%
        dplyr::filter(!any(is.na(result))) %>% 
        dplyr::ungroup()
    } else {
      stop("Missing values in test results.")
    }
  }
  
  # Summarise by test to get a per patient panel result propagating NAs
  panel_results = test_results %>% 
    dplyr::group_by(id) %>% 
    dplyr::summarise(result = any(result))
  
  # Summarise by test to get population counts of positivity per test
  # make sure they are ordered correctly
  test_counts = test_results %>% dplyr::group_by(test) %>% 
    dplyr::summarise(pos_obs = sum(result,na.rm = TRUE), n_obs=sum(!is.na(result))) 
  
  panel_counts = panel_results %>% dplyr::summarise(
    panel_pos_obs = sum(result,na.rm = TRUE), 
    panel_n_obs=sum(!is.na(result))
  )
  
  method = match.arg(method)
  if (stringr::str_starts(method,"l")) {
    
    # fill in binomial from shape paramenters if sens is given as a beta
    if ((is.null(n_diseased) || is.null(false_neg_diseased)) && !is.null(sens)) {
      if (length(sens) == 1) sens = rep(sens, length(n_test))
      n_diseased = purrr::map_dbl(sens, ~ .x$conc)
      false_neg_diseased = purrr::map_dbl(sens, ~ .x$shape2)
    }
    # fill in binomial from shape paramenters if spec is given as a beta
    if ((is.null(n_controls) || is.null(false_pos_controls)) && !is.null(spec)) {
      if (length(spec) == 1) spec = rep(spec, length(n_test))
      n_controls = purrr::map_dbl(spec, ~ .x$conc)
      false_pos_controls = purrr::map_dbl(spec, ~ .x$shape2)
    }
    
    # Lang-Reiczigel (frequentist)
    panel = prevalence_panel_lang_reiczigel(
      panel_pos_obs = panel_counts$panel_pos_obs,
      panel_n_obs = panel_counts$panel_n_obs,
      pos_obs = test_counts$pos_obs,
      n_obs = test_counts$n_obs,
      false_pos_controls = false_pos_controls,
      n_controls = n_controls,
      false_neg_diseased = false_neg_diseased,
      n_diseased = n_diseased,
      ci=confint,
      ...
    ) %>% dplyr::mutate(
      test = panel_name
    )
    
    comp = prevalence_lang_reiczigel(
      pos_obs = test_counts$pos_obs,
      n_obs = test_counts$n_obs,
      false_pos_controls = false_pos_controls,
      n_controls = n_controls,
      false_neg_diseased = false_neg_diseased,
      n_diseased = n_diseased,
      ci=confint
    ) %>% dplyr::mutate(
      test = test_names
    ) %>% dplyr::bind_cols(
      .beta_label_2(n_controls-false_pos_controls+1, false_pos_controls+1, "spec"),
      .beta_label_2(n_diseased-false_neg_diseased+1, false_neg_diseased+1, "sens"),
    )
    
    return(
      interfacer::ireturn(dplyr::bind_rows(comp,panel) %>% dplyr::mutate(
        pos_obs = c(test_counts$pos_obs,panel_counts$panel_pos_obs),
        n_obs = c(test_counts$n_obs,panel_counts$panel_n_obs)
      ), testerror:::.output_data)
    )
    
  } else if (stringr::str_starts(method,"r")) {
    
    if (is.null(sens)) sens = beta_dist(n_diseased-false_neg_diseased+1, false_neg_diseased+1)
    if (is.null(spec)) spec = beta_dist(n_controls-false_pos_controls+1, false_pos_controls+1)
    
    # Rogan Gladen
    comp = uncertain_rogan_gladen(
        pos_obs = test_counts$pos_obs,
        n_obs = test_counts$n_obs,
        sens = sens,
        spec = spec
      ) %>% 
      dplyr::mutate(
        test = test_names
      ) %>% 
      dplyr::bind_cols(
        dplyr::bind_rows(purrr::map(spec, ~ .beta_label(.x, "spec"))),
        dplyr::bind_rows(purrr::map(sens, ~ .beta_label(.x, "sens")))
      )
    
    
    panel = uncertain_panel_rogan_gladen(
      panel_pos_obs = panel_counts$panel_pos_obs,
      panel_n_obs = panel_counts$panel_n_obs,
      pos_obs = test_counts$pos_obs,
      n_obs = test_counts$n_obs,
      sens = sens,
      spec = spec,
      ...
    ) %>% dplyr::mutate(
      test = panel_name
    )
    
    return(
      interfacer::ireturn(dplyr::bind_rows(comp,panel) %>% dplyr::mutate(
        pos_obs = c(test_counts$pos_obs,panel_counts$panel_pos_obs),
        n_obs = c(test_counts$n_obs,panel_counts$panel_n_obs)
      ), testerror:::.output_data)
    )
  } else if (stringr::str_starts(method,"b")) {
    
    # Some bayesian models can have both data and priors specified for
    
    # Bayesian
    bayes = bayesian_panel_true_prevalence_model(
      
      # Some models use count data
      panel_pos_obs = panel_counts$panel_pos_obs,
      panel_n_obs = panel_counts$panel_n_obs,
      pos_obs = test_counts$pos_obs,
      n_obs = test_counts$n_obs,
      
      # other models use raw_data
      test_results = test_results,
      panel_results = panel_results,
      
      false_pos_controls = false_pos_controls,
      n_controls = n_controls,
      false_neg_diseased = false_neg_diseased,
      n_diseased = n_diseased,
      ...,
      sens = sens,
      spec = spec,
      
      ci=confint,
      test_names = levels(test_results$test),
      panel_name = panel_name
    )
    
    return(
      interfacer::ireturn(
        bayes$summary %>% dplyr::mutate(
          pos_obs = c(test_counts$pos_obs,panel_counts$panel_pos_obs),
          n_obs = c(test_counts$n_obs,panel_counts$panel_n_obs)
        ), 
        testerror:::.output_data)
    )
    
  } 
  
}


#' Rogan-Gladen true prevalence for panel with resampling
#' 
#' Uses resampling to incorporate uncertainty of sensitivity and specificity into
#' an estimate of true prevalence from a given value of apparent prevalence.
#' 
#' This is not vectorised
#' 
#' @param panel_pos_obs the number of positive observations for a given panel of tests
#' @param panel_n_obs the number of observations for each component test
#' @param pos_obs a vector of the number of positive observations for each component of a panel test
#' @param n_obs a vector of the number of observations for each component test
#' @param false_pos_controls the number of positives that appeared in the component tests of a specificity
#'   disease-free control group. These are by definition false positives. This
#'   is (1-specificity)*n_controls (ignored if `spec` is given)
#' @param n_controls the number of controls in the specificity
#'   disease-free control group (ignored if `spec` is given). 
#' @param false_neg_diseased the number of negatives that appeared in the component tests of a sensitivity
#'   confirmed disease group. These are by definition false negatives. This
#'   is (1-sensitivity)*n_diseased (ignored if `sens` is given)
#' @param n_diseased the number of confirmed disease cases in the sensitivity
#'   control group (ignored if `sens` is given)
#' @param ... not used
#' @param sens the sensitivity of the component tests of test as a `beta_dist_list` (optional)
#' @param spec the specificity of the test as a `beta_dist_list` (optional)
#' @param panel_sens the sensitivity of the component tests of test as a set of samples (optional)
#' @param panel_spec the specificity of the test as a set of samples (optional)
#' @param samples number of random draws of sensitivity and specificity (optional - default 1000)
#' 
#' @return the expected value of apparent prevalence
#' @export
uncertain_panel_rogan_gladen = function(
    panel_pos_obs,
    panel_n_obs,
    pos_obs,
    n_obs,
    false_pos_controls = NA,
    n_controls = NA,
    false_neg_diseased = NA,
    n_diseased = NA,
    ...,
    # posterior beta distribution assuming a uniform prior 
    spec = beta_dist(p=n_controls-false_pos_controls+1, q=false_pos_controls+1),
    sens = beta_dist(p=n_diseased-false_neg_diseased+1, q=false_neg_diseased+1),
    panel_sens = uncertain_panel_sens_estimator(
      pos_obs = pos_obs, 
      n_obs = n_obs,
      sens = sens,
      spec = spec,
      samples = samples,
      fit_beta = FALSE
    ), 
    panel_spec = uncertain_panel_spec(
      spec = spec, 
      samples = samples,
      fit_beta = FALSE
    ),
    samples = 1000
) { 
  
  if (inherits(panel_sens,"beta_dist")) {
    panel_sens = panel_sens$r(samples)
  }
  if (inherits(panel_spec,"beta_dist")) {
    panel_spec = panel_spec$r(samples)
  }
  
  panel = uncertain_rogan_gladen(
      panel_pos_obs, 
      panel_n_obs, 
      sens = panel_sens, 
      spec = panel_spec, 
      ...
    ) %>% 
    dplyr::bind_cols(
      .beta_label_3(panel_sens, "sens"),
      .beta_label_3(panel_spec, "spec")
    )
  
  return(panel)
}

#' Lang-Reiczigel true prevalence for panel
#' 
#' Uses resampling to incorporate uncertainty of sensitivity and specificity into
#' an estimate of true prevalence from a given value of apparent prevalence.
#' 
#' This is not vectorised
#' 
#' @param panel_pos_obs the number of positive observations for a given panel of tests
#' @param panel_n_obs a vector of the number of observations for each component test
#' @param pos_obs a vector of the number of positive observations for each component of a panel test
#' @param n_obs a vector of the number of observations for each component test
#' @param false_pos_controls the number of positives that appeared in the component tests of a specificity
#'   disease-free control group. These are by definition false positives. This
#'   is (1-specificity)*n_controls (ignored if `spec` is given)
#' @param n_controls the number of controls in the specificity
#'   disease-free control group (ignored if `spec` is given). 
#' @param false_neg_diseased the number of negatives that appeared in the component tests of a sensitivity
#'   confirmed disease group. These are by definition false negatives. This
#'   is (1-sensitivity)*n_diseased (ignored if `sens` is given)
#' @param n_diseased the number of confirmed disease cases in the sensitivity
#'   control group (ignored if `sens` is given)
#' @param ... not used
#' @param sens the sensitivity of the component tests of test as a `beta_dist_list` (optional)
#' @param spec the specificity of the test as a `beta_dist_list` (optional)
#' @param panel_sens the sensitivity of the component tests of test as a `beta_dist` (optional)
#' @param panel_spec the specificity of the test as a `beta_dist` (optional)
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
    false_pos_controls = NA,
    n_controls = NA,
    false_neg_diseased = NA,
    n_diseased = NA,
    ...,
    # posterior beta distribution assuming a uniform prior 
    spec = beta_dist(p=n_controls-false_pos_controls+1, q=false_pos_controls+1),
    sens = beta_dist(p=n_diseased-false_neg_diseased+1, q=false_neg_diseased+1),
    panel_sens = uncertain_panel_sens_estimator(
                    pos_obs = pos_obs, 
                    n_obs = n_obs,
                    sens = sens,
                    spec = spec,
                    samples = samples,
                    fit_beta = TRUE
                  ), 
    panel_spec = testerror::uncertain_panel_spec(
                    spec = spec, 
                    samples = samples,
                    fit_beta = TRUE
                  ),
    samples = 1000
) {
  
  
  # Fit sensitivity and specificity to beta distribution with moments
  if (!inherits(panel_sens,"beta_dist")) {
    panel_sens = beta_fit(panel_sens)
  }
  if (!inherits(panel_spec,"beta_dist")) {
    panel_spec = beta_fit(panel_spec)
  }
  
  # Use beta distribution as parameters to normal Lang-Recizigel
  panel = prevalence_lang_reiczigel(
    pos_obs = panel_pos_obs,
    n_obs = panel_n_obs,
    false_pos_controls = panel_spec$shape2,
    n_controls = panel_spec$conc,
    false_neg_diseased = panel_sens$shape2,
    n_diseased = panel_sens$conc,
    ...
  ) %>% dplyr::bind_cols(
    .beta_label(panel_sens, "sens"),
    .beta_label(panel_spec, "spec")
  )
  
  return(panel)
}

# decide which model to displatch to.
bayesian_panel_true_prevalence_model = function(..., model_type = c("simpler","complex")) {
  model_type = match.arg(model_type)
  if (stringr::str_starts(model_type,"c")) {
    bayesian_panel_complex_model(...)
  } else if (stringr::str_starts(model_type,"s")) {
    bayesian_panel_simpler_model(...)
  } else {
    stop("unknown stan model: ",model_type)
  }
}

# formatting functions for sens and spec ----

.beta_label = function(beta_dist, prefix, fmt = "%1.1f%% [%1.1f%% \u2013 %1.1f%%]") {
  tibble::as_tibble(beta_dist) %>%
    dplyr::select(median,lower,upper) %>%
    dplyr::mutate(label = sprintf(fmt, median*100, lower*100, upper*100)) %>%
    dplyr::rename_with(~ sprintf("%s.%s", prefix, .x))
}

.beta_label_2 = function(shape1, shape2, prefix, ci = 0.95, fmt = "%1.1f%% [%1.1f%% \u2013 %1.1f%%]") {
  tibble::tibble(
    median = stats::qbeta(0.5, shape1, shape2),
    lower = stats::qbeta((1-ci)/2, shape1, shape2),
    upper = stats::qbeta(1-(1-ci)/2, shape1, shape2)
  ) %>%
  dplyr::mutate(label = sprintf(fmt, median*100, lower*100, upper*100)) %>%
  dplyr::rename_with(~ sprintf("%s.%s", prefix, .x))
}

.beta_label_3 = function(samples, prefix, ci = 0.95, fmt = "%1.1f%% [%1.1f%% \u2013 %1.1f%%]") {
  tibble::tibble(
    median = unname(stats::quantile(samples,0.5)),
    lower = unname(stats::quantile(samples,(1-ci)/2)),
    upper = unname(stats::quantile(samples,1-(1-ci)/2))
  ) %>%
    dplyr::mutate(label = sprintf(fmt, median*100, lower*100, upper*100)) %>%
    dplyr::rename_with(~ sprintf("%s.%s", prefix, .x))
}


