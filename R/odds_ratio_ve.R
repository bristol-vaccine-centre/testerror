#' Calculate a vaccine effectiveness estimate based on an odds ratio
#' 
#' This assumes that OR ~ RR which is only true if controls >> cases
#' The OR method can be used in test negative designs where disease positive
#' relates to vaccine treatable disease and disease negative relates to non 
#' vaccine treatable disease
#'
#' @param vaccinatedCase count of disease positive vaccine positive
#' @param unvaccinatedCase count of disease positive vaccine negative
#' @param vaccinatedControl count of disease negative vaccine positive
#' @param unvaccinatedControl count of disease negative vaccine positive
#' @param confint the confidence intervals
#'
#' @return a dataframe 
#' @export
#'
#' @examples
#' tibble::tibble(
#'   N_vacc = 42240,
#'   N_unvacc = 42256,
#'   N_vacc_pn_pos = 49,
#'   N_unvacc_pn_pos = 90
#' ) %>% dplyr::mutate(
#'   odds_ratio_ve(N_vacc_pn_pos, N_unvacc_pn_pos, N_vacc-N_vacc_pn_pos, N_unvacc-N_unvacc_pn_pos)
#' )
#' 
#' 
odds_ratio_ve = function(vaccinatedCase, unvaccinatedCase, vaccinatedControl, unvaccinatedControl, confint=c(0.025,0.975)) {
  p = confint
  oddsR = dplyr::case_when(
    vaccinatedCase == 0 | unvaccinatedCase == 0 | unvaccinatedControl == 0 | vaccinatedControl == 0 ~ NA_real_,
    TRUE ~ (vaccinatedCase/vaccinatedControl) / (unvaccinatedCase/unvaccinatedControl)
  )
  logOdds = log(oddsR)
  logOddsSD = sqrt(1/vaccinatedCase+1/unvaccinatedCase+1/vaccinatedControl+1/unvaccinatedControl)
  oddsQ = purrr::map(p, ~ exp(stats::qnorm(.x,logOdds, logOddsSD)))
  names(oddsQ) = paste0("OR.q.",p)
  veQ = purrr::map(rev(p), ~ 1 - exp(stats::qnorm(.x,logOdds, logOddsSD)))
  names(veQ) = paste0("VE.OR.q.",p)
  cbind(tibble::tibble(OR = oddsR),as.data.frame(oddsQ),tibble::tibble(VE.OR = 1-oddsR),as.data.frame(veQ))
  
}