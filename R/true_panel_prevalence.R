#' Calculate an estimate of true prevalence for a single panel and components
#' 
#' Uses apparent prevalence, and uncertain
#' estimates of test sensitivity and test specificity for the 3 methods
#' described in Supplementary 2. This function works for a single panel per
#' dataframe, multiple panels will need to call this function multiple times
#' in a `group_modify`. 
#'
#' @param test_results `r interfacer::idocument(true_panel_prevalence, test_results)`
#' @inheritParams uncertain_rogan_gladen
#' @inheritDotParams uncertain_panel_rogan_gladen -pos_obs -n_obs -panel_pos_obs -panel_n_obs
#' @inheritDotParams prevalence_panel_lang_reiczigel -pos_obs -n_obs -panel_pos_obs -panel_n_obs
#' @inheritDotParams bayesian_panel_complex_model
#' @inheritDotParams bayesian_panel_true_prevalence_model 
#' @inheritDotParams bayesian_panel_simpler_model -pos_obs -n_obs -test_names -panel_pos_obs -panel_n_obs
#' @inheritDotParams bayesian_panel_logit_model -pos_obs -n_obs -test_names -panel_pos_obs -panel_n_obs
#' @param panel_name the name of the panel for combined result
#' @param method one of:
#' * "lang-reiczigel": Frequentist confidence limits
#' * "bayes": Bayesian analysis
#' * "rogan-gladen": Rogan gladen with uncertainty
#' @param na.rm exclude patients with missing results
#'
#' @return `r .output_data`
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
    test_results = testerror::.input_data,
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
    
    # if (is.null(sens)) sens = sens_prior()
    # if (is.null(spec)) spec = spec_prior()
    if (is.null(sens)) sens = uninformed_prior()
    if (is.null(spec)) spec = uninformed_prior()
    
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
      sens = sens,
      spec = spec,
      confint=confint,
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
      sens = sens,
      spec = spec,
      confint=confint
    ) %>% dplyr::mutate(
      test = test_names
    )
    
    return(
      interfacer::ireturn(dplyr::bind_rows(comp,panel) %>% dplyr::mutate(
        pos_obs = c(test_counts$pos_obs,panel_counts$panel_pos_obs),
        n_obs = c(test_counts$n_obs,panel_counts$panel_n_obs)
      ), testerror::.output_data)
    )
    
  } else if (stringr::str_starts(method,"r")) {
    
    if (is.null(sens)) sens = sens_prior()
    if (is.null(spec)) spec = spec_prior()
    # if (is.null(sens)) sens = uninformed_prior()
    # if (is.null(spec)) spec = uninformed_prior()
    
    
    # Rogan Gladen
    comp = uncertain_rogan_gladen(
        pos_obs = test_counts$pos_obs,
        n_obs = test_counts$n_obs,
        false_pos_controls = false_pos_controls,
        n_controls = n_controls,
        false_neg_diseased = false_neg_diseased,
        n_diseased = n_diseased,
        sens = sens,
        spec = spec,
        confint=confint,
        ...
      ) %>% 
      dplyr::mutate(
        test = test_names
      )
    
    panel = uncertain_panel_rogan_gladen(
        panel_pos_obs = panel_counts$panel_pos_obs,
        panel_n_obs = panel_counts$panel_n_obs,
        pos_obs = test_counts$pos_obs,
        n_obs = test_counts$n_obs,
        false_pos_controls = false_pos_controls,
        n_controls = n_controls,
        false_neg_diseased = false_neg_diseased,
        n_diseased = n_diseased,
        sens = sens,
        spec = spec,
        confint=confint,
        ...
      ) %>% dplyr::mutate(
        test = panel_name
      )
    
    return(
      interfacer::ireturn(dplyr::bind_rows(comp,panel) %>% dplyr::mutate(
        pos_obs = c(test_counts$pos_obs,panel_counts$panel_pos_obs),
        n_obs = c(test_counts$n_obs,panel_counts$panel_n_obs)
      ), testerror::.output_data)
    )
  } else if (stringr::str_starts(method,"b")) {
    
    # Some bayesian models can have both data and priors specified for
    # if (is.null(sens)) sens = uniform_prior()
    # if (is.null(spec)) spec = uniform_prior()
    # if (is.null(sens)) sens = uninformed_prior()
    # if (is.null(spec)) spec = uninformed_prior()
    # if (is.null(sens)) sens = sens_prior()
    # if (is.null(spec)) spec = spec_prior()
    
    # Bayesian
    bayes = bayesian_panel_true_prevalence_model(
      
      # Some models use count data (but will ignore additional)
      panel_pos_obs = panel_counts$panel_pos_obs,
      panel_n_obs = panel_counts$panel_n_obs,
      pos_obs = test_counts$pos_obs,
      n_obs = test_counts$n_obs,
      
      # other models use raw_data (but will ignore additional)
      test_results = test_results,
      
      false_pos_controls = false_pos_controls,
      n_controls = n_controls,
      false_neg_diseased = false_neg_diseased,
      n_diseased = n_diseased,
      ...,
      sens = sens,
      spec = spec,
      
      confint = confint,
      test_names = levels(test_results$test),
      panel_name = panel_name
    )
    
    return(
      interfacer::ireturn(
        bayes$summary %>% dplyr::mutate(
          pos_obs = c(test_counts$pos_obs,panel_counts$panel_pos_obs),
          n_obs = c(test_counts$n_obs,panel_counts$panel_n_obs)
        ), 
        testerror::.output_data)
    )
    
  } 
  
}


#' Execute one of a set of bayesian models
#'
#' @inheritDotParams bayesian_panel_complex_model
#' @inheritDotParams bayesian_panel_simpler_model
#' @inheritDotParams bayesian_panel_logit_model
#' @param model_type The bayesian model used - `r pkgutils::doc_formals(bayesian_panel_true_prevalence_model, model_type)`
#'
#' @return `r .output_data`
#' @export
bayesian_panel_true_prevalence_model = function(..., model_type = c("logit","simpler","complex")) {
  model_type = match.arg(model_type)
  if (stringr::str_starts(model_type,"c")) {
    bayesian_panel_complex_model(...)
  } else if (stringr::str_starts(model_type,"s")) {
    bayesian_panel_simpler_model(...)
  } else if (stringr::str_starts(model_type,"l")) {
    bayesian_panel_logit_model(...)
  } else {
    stop("unknown stan model: ",model_type)
  }
}



