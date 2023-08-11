# The count of positives of an individual serotype is a function of test
# specificity, prevalence, and sensitivity. For any given serotype prevalence is
# generally lower than 1-specificity and we are in the universe of the "false
# positive paradox" where there are more false positives than true positives, and
# the PPV of the test is low. The degree of false positives depends critically on the precise
# value of sensitivity and prevalence both of which are unknown.
# 
# How can we tell if an observation of a test positive count of x/n is significant?
#   
# * We assume a null hypothesis that the prevalence of this serotype is zero.
#
# * We can assume an uninformative prior on the sensitivity as a beta distribution (e.g
# beta(0.001, 0.001))
#
# * We can update that prior based on the initial cutoff determination (e.g. 2 false
# pos, 398 true neg). This will give us a ~ Beta(2.001, 398.001) posterior
# for test results.
#
# * We can use the posterior predictive distribution (a beta-binomial) to
# predict the probability of observing at least the number of observations P(X
# >= x) X ~ BetaBin(n, a, b). If this is low then we can reject null hypothesis
# that prevalence is zero. 
# 
# * The limit of the probability must be adjusted for
# multiple testing with a boneferoni adjustment for multiple serotypes.


#' Identify the minimum number of positive test result observations needed
#' to be confident the disease has a non-zero prevalence.
#'
#' @param n_obs the number of tests performed.
#' @param false_pos_controls the number of positives that appeared in the specificity
#'   disease-free control group. These are by definition false positives. This
#'   is `(1-specificity)*n_controls`
#' @param n_controls the number of controls in the specificity
#'   disease-free control group.
#' @param bonferroni the number of simultaneous tests considered.
#' @param ... not used
#'
#' @return a vector of test positive counts which are the lowest significant value
#' that could be regarded as not due to chance.
#' @export
#'
#' @examples
#' # lowest significant count of positives in 1000 tests 
#' fp_signif_level(1000, false_pos_controls = 0:5, n_controls=800)
#' fp_signif_level(c(1000,800,600,400), false_pos_controls = 1:4, n_controls=800)
fp_signif_level = function(n_obs, false_pos_controls, n_controls, bonferroni = NULL, ...) {
  
  .check_length(false_pos_controls,n_controls)
  true_neg_controls = n_controls - false_pos_controls
  
  .check_length(false_pos_controls,n_obs)
  if (length(n_obs) == 1) n_obs = rep(n_obs, length(false_pos_controls))
  if (length(false_pos_controls) == 1) {
    false_pos_controls = rep(false_pos_controls, length(n_obs))
    true_neg_controls = rep(true_neg_controls, length(n_obs))
    n_controls = rep(n_controls, length(n_obs))
  }
  
  if (is.null(bonferroni)) bonferroni = 1
  
  # assume prevalence of pneumo in controls is 0
  # assuming somewhere between 2 and 5 false positives out of 400+400 controls
  # 2 controls always positive in UAD design.
  corrected = sapply(1:length(false_pos_controls),
    function(i) {
      # A fully non-informative prior fails if pos_controls = 0 here so we use a very close to zero prior of beta(0.0001,0.0001)
      cdf = extraDistr::pbbinom(q=0:(n_obs[[i]]%/%10), n_obs[[i]], alpha = false_pos_controls[[i]]+0.0001, beta = true_neg_controls[[i]]+0.0001)
      names(cdf) = 0:(n_obs[[i]]%/%10)
      min(which(cdf > 1-(0.05/bonferroni)))-1
    })
  return(corrected)
}


.check_length = function(x,y) {
  if (length(x) != 1 && length(y) != 1 && length(x) != length(y))
    stop(deparse(substitute(x))," is incompatible length to ",deparse(substitute(y)))
}

.recycle = function(..., .min=1) {
  names = sapply(rlang::ensyms(...), rlang::as_label)
  dots = rlang::list2(...)
  lengths = sapply(dots, length)
  ml = max(c(lengths,.min))
  
  
  if (!all(lengths %in% c(0,1,ml))) {
    names = names[lengths != 1 & lengths != ml]
    stop(sprintf("%s is/are the wrong lengths. They should be 1 or %d",paste0(names,collapse=",") ,ml))
  }
  
  env = rlang::caller_env()
  
  for (i in seq_along(dots)) {
    x = dots[[i]]
    name = names[[i]]
    if (length(x) == 1) 
      env[[name]] = rep(x,ml)
  }
  
  return(ml)
  
}

#' Significance of an uncertain test result
#'
#' Calculates a p-value for a count of positive test results based
#' on false positive (specificity) controls. The null
#' hypothesis is that the prevalence of the disease is zero.
#'
#' This p_value does not tell you whether this count can be trusted only if the
#' prevalence of this disease is significantly more than zero after this
#' observation.
#' 
#' @param pos_obs the number of positive observations for a given test
#' @param n_obs the number of observations for a given test
#' @param false_pos_controls the number of positives that appeared in the specificity
#'   disease-free control group. These are by definition false positives. This
#'   is `(1-specificity)*n_controls`
#' @param n_controls the number of controls in the specificity
#'   disease-free control group.
#' @param format a sprintf fmt string for the p-value
#' @param lim a lower value to display
#' @param bonferroni the number of simultaneous hypotheses that are being tested
#' @param ... not used
#' 
#' @return a vector of p-values for the count
#' @export
#'
#' @examples
#' # calculate p-values for counts derived from 300 samples
#' # 10 observations is within noise of test
#' # 20 observations is unlikely on 1200 observations
#' fp_p_value(c(10,2,4,3,10,20), 1200, c(0,0,2,0,2,0)+2, 800)
#' 
#' # if the same observations are made against a smaller group then we get 
#' # a positive result for 10
#' fp_p_value( c(10,2,4,3,10,20), 1000, c(2,2,4,2,4,2), 800)
#' 
#' tibble::tibble(
#'   x = c(1,2,5,10,20,40,20,20,20,20,20),
#'   n = 1000,
#'   fp_controls = c(0,0,0,0,0,0,0,1,2,3,4)+2,
#'   n_controls = 800
#' ) %>% dplyr::mutate(
#'   p_value = fp_p_value(x, n, fp_controls, n_controls)
#' ) %>% dplyr::glimpse()
fp_p_value = function(pos_obs, n_obs, false_pos_controls, n_controls, format = "%1.3g", lim = 0.0001, bonferroni = NULL, ...) {
  # probabilities = extraDistr::pbbinom(q=1:80, samples, alpha = x, beta=800-x, lower.tail = FALSE)
  
  n = .recycle(pos_obs, n_obs, false_pos_controls, n_controls)
  if (is.null(bonferroni)) bonferroni = n
  true_neg_controls = n_controls - false_pos_controls
  neg_obs = n_obs - pos_obs
  
  if (length(false_pos_controls) == 1) false_pos_controls = rep(false_pos_controls, length(pos_obs))
  if (length(true_neg_controls) == 1) true_neg_controls = rep(true_neg_controls, length(pos_obs))
  
  x = extraDistr::pbbinom(q=pos_obs, pos_obs+neg_obs, alpha = false_pos_controls+1, beta = true_neg_controls+1, lower.tail = FALSE)
  
  tmp = dplyr::if_else(x<lim, sprintf(paste0("<",format),lim), sprintf(format, x))
  dplyr::if_else(x<0.05/bonferroni, sprintf("%s \u2020",tmp), tmp)
  
}

#' Vectorised true prevalence estimates
#' 
#' Calculate an estimate of true prevalence from apparent prevalence, and uncertain
#' estimates of test sensitivity and test specificity, using one of 3 methods.
#'
#' @param pos_obs the number of positive observations for a given test
#' @param n_obs the number of observations for a given test
#' @param false_pos_controls the number of positives that appeared in the specificity
#'   disease-free control group. These are by definition false positives. This
#'   is (1-specificity)*n_controls (ignored if `spec` given)
#' @param n_controls the number of controls in the specificity
#'   disease-free control group. 
#' @param false_neg_diseased the number of negatives that appeared in the sensitivity
#'   confirmed disease group. These are by definition false negatives. This
#'   is (1-sensitivity)*n_diseased (ignored if `sens` given)
#' @param n_diseased the number of confirmed disease cases in the sensitivity
#'   control group.
#' @param confint confidence limit width
#' @param method one of:
#' * "lang-reiczigel": Frequentist confidence limits
#' * "rogan-gladen": Rogan gladen incorporating uncertainty with resampling
#' * "bayes": Bayesian analysis
#' @param ... passed onto methods
#' @param sens the sensitivity of the test as a `beta_dist` (optional)
#' @param spec the specificity of the test as a `beta_dist` (optional)
#'
#' @return `r format(testerror:::.output_data)`
#' @export
#'
#' @examples
#' true_prevalence(c(1:50), 200, 2, 800, 25, 75)
#' true_prevalence(c(1:10)*2, 200, 25, 800, 4.5, 5.5)
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
    sens = NULL,
    spec = NULL
) {
  
  if (inherits(sens, "beta_dist")) sens = rep(sens,1)
  if (inherits(spec, "beta_dist")) spec = rep(spec,1)
  
  n = .recycle(pos_obs,n_obs,false_pos_controls,n_controls,
               false_neg_diseased, n_diseased, sens, spec)
  
  method = match.arg(method)
  if (!stringr::str_starts(method,"b")) {
    
    # Make sure false_pos_controls etc populated  
    if (is.null(false_pos_controls)) {
      if (is.null(spec)) stop("one of `false_pos_controls` or `spec` must be given")
      false_pos_controls = purrr::map_dbl(spec, ~ .x$shape2)
      n_controls = purrr::map_dbl(spec, ~ .x$conc)
    }
    
    if (is.null(false_neg_diseased)) {
      if(is.null(sens)) stop("one of `false_neg_diseased` or `sens` must be given")
      false_neg_diseased = purrr::map_dbl(sens, ~ .x$shape2)
      n_diseased = purrr::map_dbl(sens, ~ .x$conc)
    }
    
    if (stringr::str_starts(method,"l")) {  
      # Lang-Reiczigel (frequentist)
      return(prevalence_lang_reiczigel(
        pos_obs,n_obs,
        false_pos_controls,n_controls,
        false_neg_diseased,n_diseased,
        confint = confint, ...) %>%
        dplyr::bind_cols(
          .beta_label_2(n_controls-false_pos_controls+1, false_pos_controls+1, "spec"),
          .beta_label_2(n_diseased-false_neg_diseased+1, false_neg_diseased+1, "sens")
        ) %>% dplyr::mutate(
          pos_obs = pos_obs,
          n_obs = n_obs
        ))
    
    } else if (stringr::str_starts(method,"r")) {
      
      return(uncertain_rogan_gladen(
        pos_obs,n_obs,
        false_pos_controls,n_controls,
        false_neg_diseased,n_diseased,
        confint=confint, ...) %>%
        dplyr::bind_cols(
          .beta_label_2(n_controls-false_pos_controls+1, false_pos_controls+1, "spec"),
          .beta_label_2(n_diseased-false_neg_diseased+1, false_neg_diseased+1, "sens")
        ) %>% dplyr::mutate(
          pos_obs = pos_obs,
          n_obs = n_obs
        ))
    }
      
  } else {
    
    tmp = bayesian_component_simpler_model(
      pos_obs,n_obs,
      false_pos_controls,n_controls,
      false_neg_diseased,n_diseased,
      sens = sens,
      spec = spec,
      confint = confint, ...)
    return(
      tmp$summary %>% dplyr::mutate(
        pos_obs = pos_obs,
        n_obs = n_obs
      )
    )
  }

}

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
#' @param pos_obs the number of positive observations for a given test
#' @param n_obs the number of observations for a given test
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
#' @param ... not used
#' @param confint confidence interval limits
#' @param fmt a `sprintf` formatting string accepting 3 numbers
#' @param prefix column name prefix
#' 
#' @return the expected value of apparent prevalence
#' @export
prevalence_lang_reiczigel = function(
    pos_obs,
    n_obs,
    false_pos_controls,
    n_controls,
    false_neg_diseased,
    n_diseased,
    ...,
    confint=0.95,
    prefix = "prevalence",
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]"
    ) # Confidence level
{
  
  nprev = n_obs                               # Sample size for prevalence
  kprev = pos_obs                             # Frequency of positive diagnoses in sample of size nprev
  nspec = n_controls                          # Sample size for specificity
  kspec = n_controls-false_pos_controls       # Frequency of negative diagnoses in sample of size nspec
  nsens = n_diseased                          # Sample size for sensitivity
  ksens = n_diseased-false_neg_diseased       # Frequency of positive diagnoses in sample of size nsens
  
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
    lower = LCL, 
    median = est.prev,
    upper = UCL,
    method = "lang-reiczigel"
  ) %>% dplyr::mutate(
    label = sprintf(fmt, median*100,lower*100,upper*100)
  )
  
  if (!is.null(prefix)) tmp = tmp %>% dplyr::rename_with(~ sprintf("%s.%s",prefix,.x))
  
  return(tmp)
}




#' True prevalence from apparent prevalence with uncertainty
#' 
#' Uses resampling to incorporate uncertainty of sensitivity and specificity into
#' an estimate of true prevalence from a given value of apparent prevalence.
#' 
#' This is not vectorised
#' 
#' @param pos_obs the number of positive observations for a given test
#' @param n_obs the number of observations for a given test
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
#' @param ... not used
#' @param sens the sensitivity of the test as a `beta_dist` or as a vector of samples
#' @param spec the specificity of the test as a `beta_dist` or as a vector of samples
#' @param samples number fo random draws of sensitivity and specificity
#' @param confint confidence interval limits
#' @param fmt a `sprintf` formatting string accepting 3 numbers
#' @param seed set seed for reproducibility
#' @param prefix column name prefix
#' 
#' @return the expected value of apparent prevalence
#' @export
#'
#' @examples
#' uncertain_rogan_gladen(
#'   pos_obs = 20, n_obs = 1000, 
#'   false_pos_controls = 10, n_controls = 800, 
#'   false_neg_diseased = 20, n_diseased = 100)
#'   
#' uncertain_rogan_gladen(
#'   pos_obs = 5, n_obs = 1000, 
#'   sens = beta_dist(0.75,n=200), 
#'   spec = beta_dist(0.9975, n=800))
#' 
#' uncertain_rogan_gladen(
#'   pos_obs = c(5,10), n_obs = c(1000,1000), 
#'   false_pos_controls = c(2,1), n_controls = c(800,800), 
#'   false_neg_diseased = c(25,20),n_diseased = c(100,100))
#' 
uncertain_rogan_gladen = function(
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
    samples = 1000,
    confint = 0.95,
    fmt = "%1.2f%% [%1.2f%% \u2014 %1.2f%%]",
    seed = NA,
    prefix = "prevalence"
) { 
  
  if (!is.na(seed)) set.seed(seed)
  
  if (inherits(sens,"beta_dist_list")) {
    return(dplyr::bind_rows(lapply(1:length(sens), function(i) {
      sens2 = sens[[i]]
      spec2 = spec[[i]]
      return(uncertain_rogan_gladen(pos_obs = pos_obs[[i]], n_obs = n_obs[[i]], sens = sens2, spec = spec2))
    })))
  }
  
  if (inherits(sens,"beta_dist")) {
    sens = sens$r(samples)
  } else {
    sens = unlist(sens)
    samples = length(sens)
  }
  
  if (inherits(spec,"beta_dist")) {
    spec = spec$r(samples)
  } else {
    spec = unlist(spec)
    if (length(spec) != samples) stop("mismatch in length of samples in sens and spec")
  }
  
  ap = stats::rbeta(samples, pos_obs+1, n_obs-pos_obs+1)  
  # ap = rep(pos_obs/n_obs,samples)
  p = rogan_gladen(ap,sens,spec)
  zcrit = (1-confint)/2
  if (is.na(confint)) return(p)
  median = stats::quantile(p, 0.5)
  lower = stats::quantile(p, zcrit)
  upper = stats::quantile(p, 1-zcrit)
  label = sprintf(fmt, median*100,lower*100,upper*100)
  out = tibble::tibble(
    median=median, lower = lower, upper=upper, label=label, method = "rogan-gladen (samples)"
  )
  if (!is.null(prefix)) out = out %>% dplyr::rename_with(~ sprintf("%s.%s",prefix,.x))
  return(out)
}
