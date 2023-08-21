#' Vectorised true prevalence estimates
#' 
#' Calculate an estimate of true prevalence from apparent prevalence, and uncertain
#' estimates of test sensitivity and test specificity, using one of 3 methods.
#'
#' @inheritParams uncertain_rogan_gladen
#' @inheritDotParams uncertain_rogan_gladen
#' @param method one of:
#' * "lang-reiczigel": Frequentist confidence limits: see `prevalence_lang_reiczigel()`
#' * "rogan-gladen": Rogan gladen incorporating uncertainty with resampling: see `uncertain_rogan_gladen()`
#' * "bayes": Bayesian analysis: see `bayesian_component_simpler_model()`
#'
#' @return `r .output_data`
#' @export
#'
#' @examples
#' true_prevalence(c(1:50), 200, 2, 800, 25, 75)
#' true_prevalence(c(1:10)*2, 200, 25, 800, 1, 6, method="rogan-gladen")
#' true_prevalence(c(1:10)*2, 200, 5, 800, 1, 6, method="bayes")
true_prevalence = function(
    pos_obs,
    n_obs,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    confint=.95,
    method = c("lang-reiczigel", "rogan-gladen", "bayes"),
    ...,
    spec = NULL,
    sens = NULL
) {
  
  n = pkgutils::recycle(pos_obs,n_obs,false_pos_controls,n_controls,
               false_neg_diseased, n_diseased, sens, spec)
  method = match.arg(method)
  
  if (!stringr::str_starts(method,"b")) {
    
    if (is.null(sens)) sens = sens_prior()
    if (is.null(spec)) spec = spec_prior()
    
    if (stringr::str_starts(method,"l")) {  
      # Lang-Reiczigel (frequentist)
      return(prevalence_lang_reiczigel(
        pos_obs,n_obs,
        false_pos_controls = false_pos_controls,
        n_controls = n_controls,
        false_neg_diseased = false_neg_diseased,
        n_diseased = n_diseased,
        sens = sens,
        spec = spec,
        confint=confint, 
        ...) %>%
        dplyr::mutate(
          pos_obs = pos_obs,
          n_obs = n_obs
        ))
      
    } else if (stringr::str_starts(method,"r")) {
      
      return(uncertain_rogan_gladen(
        pos_obs,n_obs,
        false_pos_controls = false_pos_controls,
        n_controls = n_controls,
        false_neg_diseased = false_neg_diseased,
        n_diseased = n_diseased,
        sens = sens,
        spec = spec,
        confint=confint, 
        ...) %>% 
        dplyr::mutate(
          pos_obs = pos_obs,
          n_obs = n_obs
        ))
    }
    
  } else {
    
    if (is.null(sens)) sens = uniform_prior()
    if (is.null(spec)) spec = uniform_prior()
    
    tmp = bayesian_component_simpler_model(
      pos_obs,n_obs,
      false_pos_controls = false_pos_controls,
      n_controls = n_controls,
      false_neg_diseased = false_neg_diseased,
      n_diseased = n_diseased,
      sens = sens,
      spec = spec,
      confint=confint, 
      ...)
    return(
      tmp$summary %>% dplyr::mutate(
        pos_obs = pos_obs,
        n_obs = n_obs
      )
    )
  }
  
}

