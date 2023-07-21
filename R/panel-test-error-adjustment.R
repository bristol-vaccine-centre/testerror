#' Calculate an estimate of true prevalence for a single panel and components
#' 
#' Uses apparent prevalence, and uncertain
#' estimates of test sensitivity and test specificity for the 3 methods
#' described in Supplementary 2. This function works for a single panel per
#' dataframe, mulitple panels will need to call this function multiple times
#' in a `group_modify`. 
#'
#' @param test_results `r interfacer::idocument(true_panel_prevalence, test_results)`
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
#' @param ... passed onto methods
#' @param sens the sensitivity of the component tests of test as a `beta_dist_list`
#' @param spec the specificity of the test as a `beta_dist_list`
#' @param test_names a vector of the component test names in desired order (optional)
#' @param panel_name the name of the panel for combined result
#' @param confint confidence limit width
#' @param method one of:
#' * "lang-reiczigel": Frequentist confidence limits
#' * "bayes": Bayesian analysis
#' * "rogan-gladen": Rogan gladen with uncertainty
#'
#' @return a dataframe with prevalence.lower, prevalence.median and prevalence.upper columns
#' @export
#'
#' @examples
#' tmp = random_example()
#' true_panel_prevalence(
#'   test_results = tmp$samples %>% select(id,test,result = observed),
#'   false_pos_controls = tmp$performance$false_pos_controls,
#'   n_controls = tmp$performance$n_controls,
#'   false_neg_diseased = tmp$performance$false_neg_diseased,
#'   n_diseased = tmp$performance$n_diseased,
#'   method = "rogan-gladen"
#' )
true_panel_prevalence = function(
    test_results = interfacer::iface(
          id = character ~ "the patient identifier",
          test = factor ~ "the test type",
          result = logical ~ "the test result",
          .groups = FALSE
    ),
    false_pos_controls = NA,
    n_controls = NA,
    false_neg_diseased = NA,
    n_diseased = NA,
    ...,
    sens = beta_dist(n_diseased-false_neg_diseased, false_neg_diseased), 
    spec = beta_dist(n_controls-false_pos_controls, false_pos_controls),
    test_names = levels(test_results$test),
    panel_name = "Panel",
    confint=.95,
    method = c("rogan-gladen", "lang-reiczigel", "bayes")
) {
  
  test_results = interfacer::ivalidate(test_results, ..., .prune=TRUE)
  n_test = length(test_names)
  if (length(levels(test_results$test)) != n_test) stop("The test results dataframe does not have the same length as the test names.")
  
  # Fill in missing test results
  test_results = test_results %>% tidyr::complete(test, fill = list(result=NA))
  
  # Summarise by test to get a per patient panel result propagating NAs
  panel_results = test_results %>% group_by(id) %>% 
    summarise(result = any(result))
  
  # Summarise by test to get population counts of positivity per test
  # make sure they are ordered correctly
  test_counts = test_results %>% group_by(test) %>% 
    summarise(pos_obs = sum(result,na.rm = TRUE), n_obs=sum(!is.na(result))) 
  
  panel_counts = panel_results %>% summarise(
    panel_pos_obs = sum(result,na.rm = TRUE), 
    panel_n_obs=sum(!is.na(result))
  )
  
  if (length(sens) == 1) sens = rep(sens, length(n_test))
  if (length(spec) == 1) spec = rep(spec, length(n_test))
  
  if (!(
    length(sens) == n_test && length(spec) == n_test
  )) stop("Mismatch between number of tests and lengths of sensitivity or specificity parameters")
  
    method = match.arg(method)
  if (stringr::str_starts(method,"l")) {
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
      ci=confint
    ) %>% mutate(
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
    ) %>%
    mutate(
      test = test_names
    )
    return(bind_rows(comp,panel))
    
  } else if (stringr::str_starts(method,"b")) {
    # Bayesian
    stop("not yet implemented")
    
    
    
  } else if (stringr::str_starts(method,"r")) {
    # Rogan Gladen
    comp = uncertain_rogan_gladen(
      pos_obs = test_counts$pos_obs,
      n_obs = test_counts$n_obs,
      sens = sens,
      spec = spec
    ) %>% 
      bind_rows() %>%
      mutate(
        test = test_names
      )
    panel = uncertain_panel_rogan_gladen(
      panel_pos_obs = panel_counts$panel_pos_obs,
      panel_n_obs = panel_counts$panel_n_obs,
      pos_obs = test_counts$pos_obs,
      n_obs = test_counts$n_obs,
      sens = sens,
      spec = spec
    ) %>% mutate(
      test = panel_name
    )
    return(bind_rows(comp,panel))
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
#' @param ... 
#' @param sens the sensitivity of the component tests of test as a `beta_dist_list` (optional)
#' @param spec the specificity of the test as a `beta_dist_list` (optional)
#' @param samples number of random draws of sensitivity and specificity (optional - default 1000)
#' 
#' @return the expected value of apparent prevalence
#' @export
#'
#' @examples
#' #TODO uncertain_panel_rogan_gladen(pos_obs = 20, n_obs = 1000, false_pos_controls = 10, n_controls = 800, false_neg_diseased = 20, n_diseased = 100)
#' #TODO uncertain_rogan_gladen(pos_obs = 5, n_obs = 1000, sens = beta_dist(0.75,n=200), spec = beta_dist(0.9975, n=800))
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
    sens = beta_dist(n_diseased-false_neg_diseased, false_neg_diseased), 
    spec = beta_dist(n_controls-false_pos_controls, false_pos_controls),
    samples = 1000
) { 
  
  panel_spec = uncertain_panel_spec(spec=spec, samples = samples)
  panel_sens = uncertain_panel_sens_estimator(
    pos_obs=pos_obs, n_obs=n_obs,
    sens = sens, spec=spec, samples = samples)
  
  uncertain_rogan_gladen(panel_pos_obs, panel_n_obs, sens = panel_sens, spec = panel_spec, ...)
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
#' @param ... 
#' @param sens the sensitivity of the component tests of test as a `beta_dist_list` (optional)
#' @param spec the specificity of the test as a `beta_dist_list` (optional)
#' @param samples number of random draws of sensitivity and specificity (optional - default 1000)
#' 
#' @return the expected value of apparent prevalence
#' @export
#'
#' @examples
#' #TODO uncertain_panel_rogan_gladen(pos_obs = 20, n_obs = 1000, false_pos_controls = 10, n_controls = 800, false_neg_diseased = 20, n_diseased = 100)
#' #TODO uncertain_rogan_gladen(pos_obs = 5, n_obs = 1000, sens = beta_dist(0.75,n=200), spec = beta_dist(0.9975, n=800))
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
    sens = beta_dist(n_diseased-false_neg_diseased, false_neg_diseased), 
    spec = beta_dist(n_controls-false_pos_controls, false_pos_controls),
    samples = 1000
) {
  # Get the panel sensitivity as a set of samples
  panel_sens_samples = testerror::uncertain_panel_sens_estimator(
    pos_obs = pos_obs, 
    n_obs = n_obs,
    sens = sens,
    spec = spec,
    samples = samples
  )
  
  # Get the panel specificity as a set of samples
  panel_spec_samples = testerror::uncertain_panel_spec(
    spec = spec, 
    samples = samples
  )
  
  # Fit sensitivity and specificity to beta distribution with moments
  panel_sens_beta = beta_fit(panel_sens_samples)
  panel_spec_beta = beta_fit(panel_spec_samples)
  
  # Use beta distribution as parameters to normal Lang-Recizigel
  panel = prevalence_lang_reiczigel(
    pos_obs = panel_pos_obs,
    n_obs = panel_n_obs,
    false_pos_controls = panel_spec_beta$shape2,
    n_controls = panel_spec_beta$conc,
    false_neg_diseased = panel_sens_beta$shape2,
    n_diseased = panel_sens_beta$conc,
    ...
  ) 
  
  return(panel)
}