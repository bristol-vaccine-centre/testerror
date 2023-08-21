#' True prevalence from apparent prevalence with uncertainty
#' 
#' Uses lang-reiczigel estimators to incorporate uncertainty of sensitivity and specificity into
#' an estimate of true prevalence from a given value of apparent prevalence.
#' 
# Reference: Lang Zs, Reiczigel J (2014) Confidence limits for prevalence of disease 
#            adjusted for estimated sensitivity and specificity, Preventive Veterinary 
#            Medicine 113, 13-22.
# 
# Function adapted by Robert Challen from script at
# Function adapted by Matthias Flor from script at
# http://www2.univet.hu/users/jreiczig/CI4prevSeSp/CI4TruePrevalence_15_03_2013.r
#' 
#' @inheritParams uncertain_rogan_gladen
#' 
#' @return the expected value of apparent prevalence
#' @export
prevalence_lang_reiczigel = function(
    pos_obs,
    n_obs,
    false_pos_controls = NULL,
    n_controls = NULL,
    false_neg_diseased = NULL,
    n_diseased = NULL,
    ...,
    spec = spec_prior(),
    sens = sens_prior(),
    confint = 0.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]"
) {
  
  # N.b. this is vectorised
  
  n = pkgutils::recycle(pos_obs,n_obs,false_pos_controls,n_controls,
                        false_neg_diseased, n_diseased, sens, spec)
  
  pkgutils::check_integer(pos_obs, n_obs)
  pkgutils::check_integer(false_pos_controls, n_controls, false_neg_diseased, n_diseased)
  
  sens = update_posterior(sens, neg=false_neg_diseased, n=n_diseased)
  spec = update_posterior(spec, neg=false_pos_controls, n=n_controls)
  
  nprev = n_obs                               # Sample size for prevalence
  kprev = pos_obs                             # Frequency of positive diagnoses in sample of size nprev
  nspec = get_beta_shape(spec, "conc")        # Sample size for specificity
  kspec = get_beta_shape(spec, "shape1")      # Frequency of negative diagnoses in sample of size nspec
  nsens = get_beta_shape(sens, "conc")        # Sample size for sensitivity
  ksens = get_beta_shape(sens, "shape1")      # Frequency of positive diagnoses in sample of size nsens
  
  # TODO: check if this is part of original. 
  # Maybe not required. maybe inserted to make sure that test with zero
  # observations do not cause NaN. 
  nprev = dplyr::if_else(nprev==0, 1, nprev)
  nsens = dplyr::if_else(nsens==0, 1, nsens)
  nspec = dplyr::if_else(nspec==0, 1, nspec)
  
  # Observed relative frequencies
  obs.prev = kprev/nprev
  obs.sens = ksens/nsens
  obs.spec = kspec/nspec
  
  # MF inserted the following line:
  # if (obs.sens + obs.spec <= 1) return(tibble::tibble(lower=NA_real_,median=NA_real_, upper=NA_real_))
  
  # Rogan-Gladen point estimate of true prevalence
  est.prev = (obs.prev+obs.spec-1)/(obs.sens+obs.spec-1)
  est.prev = scales::squish(est.prev)
  
  # Adjustments
  zcrit = stats::qnorm((1+confint)/2)
  plus  = 2
  
  nprev. = nprev+zcrit^2
  kprev. = kprev+zcrit^2/2
  
  nsens. = nsens+plus
  ksens. = ksens+plus/2
  
  nspec. = nspec+plus
  kspec. = kspec+plus/2
  
  obs.prev. = kprev./nprev.
  obs.sens. = ksens./nsens.
  obs.spec. = kspec./nspec.
  
  
  est.prev. = (obs.prev.+obs.spec.-1)/(obs.sens.+obs.spec.-1)
  
  # Youden index
  Youden. = obs.sens.+obs.spec.-1
  
  # Standard error of est.prev.
  se.est.prev. = sqrt(
    obs.prev.*(1-obs.prev.)/nprev. +
      obs.sens.*(1-obs.sens.)/nsens. * est.prev.^2 +
      obs.spec.*(1-obs.spec.)/nspec. * (1-est.prev.)^2
  )/abs(Youden.)
  
  # Shift parameter
  dprev = 2*zcrit^2*
    (est.prev.*obs.sens.*(1-obs.sens.)/nsens. - (1-est.prev.)*obs.spec.*(1-obs.spec.)/nspec.)
  
  # Adjusted confidence limits
  LCL = est.prev.+dprev - zcrit*se.est.prev.
  UCL = est.prev.+dprev + zcrit*se.est.prev.
  LCL = scales::squish(LCL)
  UCL = scales::squish(UCL)
  
  #     cat("\nSensitivity: ", round(obs.sens, 4), ", adjusted: ", round(obs.sens., 4),
  #         "\nSpecificity: ", round(obs.spec, 4), ", adjusted: ", round(obs.spec., 4),
  #         "\nObserved prevalence: ", round(obs.prev, 8), ", adjusted: ", round(obs.prev., 8),
  #         "\nRogan-Gladen true prevalence: ", round(est.prev,4), "\n",
  #         "Adjusted ", round(100*conflevel), "% CI: ", round(LCL,4), " - ", round(UCL,4), "\n",
  #         sep="")
  tmp = tibble::tibble(
    # apparent = obs.prev,
    prevalence.lower = LCL, 
    prevalence.median = est.prev,
    prevalence.upper = UCL,
    prevalence.method = "lang-reiczigel"
  ) %>% 
  dplyr::mutate(
    prevalence.label = sprintf(fmt, 
      prevalence.median*100,
      prevalence.lower*100,
      prevalence.upper*100)
  ) %>% 
  dplyr::bind_cols(
    .beta_label(spec, "spec", ci = confint, fmt=fmt),
    .beta_label(sens, "sens", ci = confint, fmt=fmt),
  )
  
  return(tmp)
}




