# adaptor code to execute the range of bayesian models included

# model caching ---- 

#' @importFrom memoise memoise
#' @importFrom digest digest
#' @importFrom rappdirs user_cache_dir
#' @import BH 
#' @importFrom Rcpp getRcppVersion
#' @import RcppEigen
#' @import RcppParallel
#' @import StanHeaders
NULL

.sampling_cached = memoise::memoise(rstan::sampling)

.get_stan_model = function(name, .nocache = FALSE) {
  
  file = fs::path_ext_set(fs::path(system.file("stan", package="testerror"),name),"stan")
  path = fs::path_ext_set(fs::path(rappdirs::user_cache_dir("testerror"),"models",rstan::stan_version(),name),"stan")
  fs::dir_create(fs::path_dir(path))
  if (!fs::file_exists(path) || .nocache) fs::file_copy(file, path,overwrite = TRUE)
  
  # must call the stan_model in the global environment to allow caching to occur
  tmp = do.call(rstan::stan_model, args = list(
    file = path,
    model_name = name,
    auto_write = TRUE
  ), envir = globalenv())
  
  return(tmp)
}

.get_stan_model_names = function() {
  paths = fs::dir_ls(system.file("stan", package="testerror"),glob = "*.stan")
  paths %>% fs::path_file() %>% fs::path_ext_remove()
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
    false_pos_controls,
    n_controls,
    false_neg_diseased,
    n_diseased,
    sens = NULL, 
    spec = NULL,
    panel_sens = NULL,
    panel_spec = NULL,
    logit = FALSE,
    ...
) {
  
  standata = list()
  
  # TODO remove redundant parts due to sens and spec being always present
  
  # Spec priors
  if (!is.null(spec)) {
    if (!logit) {
      standata$tn_spec_prior = array(as.double(get_beta_shape(spec, "shape1")), dim=n_test)
      standata$fp_spec_prior = array(as.double(get_beta_shape(spec, "shape2")), dim=n_test)
    } else {
      tmp = beta_dist_to_logitnorm(spec)
      standata$logit_mu_spec_prior = array(tmp$mu, dim=n_test)
      standata$logit_sigma_spec_prior = array(tmp$sigma, dim=n_test)
    }
  }
  
  # Sens priors
  if (!is.null(sens)) {
    if (!logit) {
      standata$tp_sens_prior = array(as.double(get_beta_shape(sens, "shape1")), dim=n_test)
      standata$fn_sens_prior = array(as.double(get_beta_shape(sens, "shape2")), dim=n_test)
    } else {
      tmp = beta_dist_to_logitnorm(sens)
      standata$logit_mu_sens_prior = array(tmp$mu, dim=n_test)
      standata$logit_sigma_sens_prior = array(tmp$sigma, dim=n_test)
    }
  }
  
  # negative controls (spec)
  if (!is.null(false_pos_controls)) {
    # Use fp_spec as binomial outcome
    standata$k_spec = array(as.integer(n_controls), dim=n_test)
    standata$fp_spec = array(as.integer(false_pos_controls), dim=n_test)
  } else {
    standata$k_spec = array(rep(0L,n_test), dim=n_test)
    standata$fp_spec = array(rep(0L,n_test), dim=n_test)
  }
  
  # positive controls (sens)  
  if (!is.null(false_neg_diseased)) {
    standata$k_sens = array(as.integer(n_diseased), dim=n_test)
    standata$fn_sens = array(as.integer(false_neg_diseased), dim=n_test)
  } else {
    standata$k_sens = array(rep(0L,n_test), dim=n_test)
    standata$fn_sens = array(rep(0L,n_test), dim=n_test)
  }
  
  # Panel sens priors
  if (!is.null(panel_sens)) {
    if (!logit) {
      standata$tp_panel_sens_prior = as.double(panel_sens$shape1)
      standata$fn_panel_sens_prior = as.double(panel_sens$shape2)
    } else {
      tmp = beta_dist_to_logitnorm(panel_sens)
      standata$logit_mu_panel_sens_prior = tmp$mu
      standata$logit_sigma_panel_sens_prior = tmp$sigma
    }
  }
  
  # Panel spec priors
  if (!is.null(panel_spec)) {
    if (!logit) {
      standata$tn_panel_spec_prior = as.double(panel_spec$shape1)
      standata$fp_panel_spec_prior = as.double(panel_spec$shape2)
    } else {
      tmp = beta_dist_to_logitnorm(panel_spec)
      standata$logit_mu_panel_spec_prior = tmp$mu
      standata$logit_sigma_panel_spec_prior = tmp$sigma
    }
  } 
  
  return(standata)
}