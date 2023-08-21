# Internal functions used in some vignettes and examples

# distribution functions ----

# generate a theoretical distribution of component prevalences that combine to produce
# a panel prevalence at a given level. Based on a relative frequency of observation.
distribute = function(dist, p) {
  dist = dist/sum(dist)
  p = unique(p)
  if (length(p) > 1) stop("p must be unique")
  u = stats::uniroot(f = function(k) 1-prod(1-dist^k*p)-p,interval = c(0,.Machine$double.max.exp))
  tmp = u$root
  out = dist^tmp*p
  if (abs(1-prod(1-out)-p) > .Machine$double.eps^0.25 ) stop("did not find a distribution")
  return(out)
}

# # rfixed2(900,c(1,0,0.5))
# rfixed2 = function(n, prob) {
#   s = n/length(prob)
#   if (s != round(s)) stop("`n` must be a whole multiple of the number of probabilities (",length(prob),")")
#   pos = round(s*prob)
#   tmp2 = lapply(pos, function(p) {
#     neg = round(s-p)
#     tmp = c(rep(0,neg),rep(1,p))
#     sample(tmp)
#   })
#   browser()
# }

# Create a sample with exactly n*prev positives.
rfixed = function(boots, n, prev) {
  prev=unique(prev)
  if (length(prev) != 1) stop("prev must be a unique value over the group (maybe use group by)")
  # if (n == 1) stop("sample size must be large enough to make a meaningful group")
  pos = round(n*prev)
  neg = round(n-pos)
  s = c(rep(0,neg),rep(1,pos))
  lapply(1:boots,  function(...) sample(s)) %>% unlist()
}


# pipeline function ----

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

# formatting functions for sens and spec ----

.beta_label = function(beta_dist, prefix, ci=0.95, fmt = "%1.1f%% [%1.1f%% \u2013 %1.1f%%]") {
  tibble::as_tibble(beta_dist, confint=ci) %>%
    dplyr::select(median,lower,upper) %>%
    dplyr::mutate(label = sprintf(fmt, median*100, lower*100, upper*100)) %>%
    dplyr::rename_with(~ sprintf("%s.%s", prefix, .x))
}

# .beta_label_2 = function(shape1, shape2, prefix, ci = 0.95, fmt = "%1.1f%% [%1.1f%% \u2013 %1.1f%%]") {
#   tibble::tibble(
#     median = stats::qbeta(0.5, shape1, shape2),
#     lower = stats::qbeta((1-ci)/2, shape1, shape2),
#     upper = stats::qbeta(1-(1-ci)/2, shape1, shape2)
#   ) %>%
#     dplyr::mutate(label = sprintf(fmt, median*100, lower*100, upper*100)) %>%
#     dplyr::rename_with(~ sprintf("%s.%s", prefix, .x))
# }

.beta_label_3 = function(samples, prefix, ci = 0.95, fmt = "%1.1f%% [%1.1f%% \u2013 %1.1f%%]") {
  tibble::tibble(
    median = unname(stats::quantile(samples,0.5)),
    lower = unname(stats::quantile(samples,(1-ci)/2)),
    upper = unname(stats::quantile(samples,1-(1-ci)/2))
  ) %>%
    dplyr::mutate(label = sprintf(fmt, median*100, lower*100, upper*100)) %>%
    dplyr::rename_with(~ sprintf("%s.%s", prefix, .x))
}


# fp to spec etc ----

false_pos = function(n, spec) {
  n*(1-spec)
}

true_neg = function(n, spec) {
  n*spec
}

false_neg = function(n, sens) {
  n*(1-sens)
}

true_pos = function(n, sens) {
  n*sens
}
