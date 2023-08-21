#' Test underestimation limit
#'
#' For a given sensitivity and specificity this give the critical threshold after which 
#' test error introduces underestimation rather than over estimation
#'
#' @param sens the sensitivity of the test
#' @param spec the specificity of the test
#'
#' @return the value where apparent prevalence equals true prevalence
#' @export
#'
#' @examples
#' tmp1 = underestimation_threshold(0.75, 0.97)
#' tmp2 = rogan_gladen(tmp1, 0.75, 0.97)
#' if (abs(tmp1-tmp2) > 0.0000001) stop("error")
#' 
underestimation_threshold = function(sens, spec) {
  return((1-spec) / (1-sens+1-spec))
}

