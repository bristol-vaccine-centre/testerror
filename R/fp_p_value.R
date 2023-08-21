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
  
  n = pkgutils::recycle(pos_obs, n_obs, false_pos_controls, n_controls)
  if (is.null(bonferroni)) bonferroni = n
  true_neg_controls = n_controls - false_pos_controls
  neg_obs = n_obs - pos_obs
  
  if (length(false_pos_controls) == 1) false_pos_controls = rep(false_pos_controls, length(pos_obs))
  if (length(true_neg_controls) == 1) true_neg_controls = rep(true_neg_controls, length(pos_obs))
  
  x = extraDistr::pbbinom(q=pos_obs, pos_obs+neg_obs, alpha = false_pos_controls+1, beta = true_neg_controls+1, lower.tail = FALSE)
  
  tmp = dplyr::if_else(x<lim, sprintf(paste0("<",format),lim), sprintf(format, x))
  dplyr::if_else(x<0.05/bonferroni, sprintf("%s \u2020",tmp), tmp)
  
}

