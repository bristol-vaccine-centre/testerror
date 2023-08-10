either_or = function(.data, condition, if_true, if_false, ...) {
  if (condition) 
    return(purrr::as_mapper(if_true)(.data, ...))
  else 
    return(purrr::as_mapper(if_false)(.data, ...))
}
 
# flag = FALSE
# diamonds %>% either_or(flag,
#   ~ .x %>% dplyr::group_by(cut),
#   ~ .x %>% dplyr::group_by(color)
# ) %>% dplyr::count()
