#' Dataframe format for component test results
#' 
#' 
#' `r .input_data`
#'
#' @name .input_data
#' @docType data
#' @keywords data
#' @export
.input_data = interfacer::iface(
  id = character ~ "the patient identifier",
  test = factor ~ "the test type",
  result = logical ~ "the test result",
  .groups = FALSE
)

#' Dataframe format for panel test results
#' 
#'  
#' `r .input_panel_data`
#'
#' @name .input_panel_data
#' @docType data
#' @keywords data
#' @export
.input_panel_data = interfacer::iface(
  id = character ~ "the patient identifier",
  result = logical ~ "the panel result",
  .groups = FALSE
)

#' Dataframe format for true prevalence results
#' 
#'  
#' `r .output_data`
#'
#' @name .output_data
#' @docType data
#' @keywords data
#' @export
.output_data = interfacer::iface(
  test = character ~ "the name of the test or panel",
  prevalence.lower = numeric ~ "the lower estimate",
  prevalence.median = numeric ~ "the median estimate",
  prevalence.upper = numeric ~ "the upper estimate",
  prevalence.method = character ~ "the method of estimation",
  prevalence.label = character ~ "a fomatted label of the true prevalence estimate with CI",
  .groups = FALSE
)