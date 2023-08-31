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
#' @param spec a prior value for specificity as a `beta`
#'
#' @return a vector of test positive counts which are the lowest significant value
#' that could be regarded as not due to chance.
#' @export
#'
#' @examples
#' # lowest significant count of positives in 1000 tests 
#' fp_signif_level(1000, false_pos_controls = 0:5, n_controls=800)
#' fp_signif_level(c(1000,800,600,400), false_pos_controls = 1:4, n_controls=800)
fp_signif_level = function(n_obs, false_pos_controls, n_controls, bonferroni = NULL, ..., spec = NULL) {
  
  if (is.null(spec)) spec = uninformed_prior()
    
  pkgutils::recycle(n_obs, false_pos_controls, n_controls, spec)
  true_neg_controls = n_controls - false_pos_controls
  
  
  
  if (is.null(bonferroni)) bonferroni = 1
  
  # assume prevalence of pneumo in controls is 0
  # assuming somewhere between 2 and 5 false positives out of 400+400 controls
  # 2 controls always positive in UAD design.
  corrected = sapply(1:length(false_pos_controls),
    function(i) {
      # A fully non-informative prior fails if pos_controls = 0 here so we use a very close to zero prior of beta(0.0001,0.0001)
      cdf = extraDistr::pbbinom(
        q=0:(n_obs[[i]]%/%10), n_obs[[i]], 
        alpha = false_pos_controls[[i]]+get_beta_shape(spec,type = "shape1"), 
        beta = true_neg_controls[[i]]+get_beta_shape(spec,type = "shape2")
      )
      names(cdf) = 0:(n_obs[[i]]%/%10)
      min(which(cdf > 1-(0.05/bonferroni)))-1
    })
  return(corrected)
}



