distribute = function(dist, p) {
  dist = dist/sum(dist)
  p = unique(p)
  if (length(p) > 1) stop("p must be unique")
  u = uniroot(f = function(k) 1-prod(1-dist^k*p)-p,interval = c(0,.Machine$double.max.exp))
  tmp = u$root
  out = dist^tmp*p
  if (abs(1-prod(1-out)-p) > .Machine$double.eps^0.25 ) stop("did not find a distribution")
  return(out)
}

random_example = function(
    panel_prev = 0.1,
    n_comp = floor(runif(1, 5, 11)), 
    comp_dist = rpois(n_comp,10) * c(1,rpois(n_comp-1,2)),
    comp_spec = rbeta(n_comp, 995, 5),
    comp_sens = rbeta(n_comp, 40, 10),
    n_samples = 1000,
    n_controls = 800,
    n_diseased = 30
) {
  
  comp_test = factor(sprintf("Component %s",LETTERS[1:n_comp]))
  comp_prev = distribute(comp_dist, panel_prev)
  
  comp_design = tibble::tibble(
    test = comp_test,
    design_prev = comp_prev,
    design_spec = comp_spec,
    design_sens = comp_sens
  ) 
  
  panel_design = comp_design %>%
    summarise(
      panel_design_prev = panel_prevalence(design_prev),
      panel_design_spec = panel_spec(design_spec),
      panel_design_sens = panel_sens(design_prev,design_sens,design_spec)
    ) %>% mutate(test = "Panel") %>%
    rename_with(~ stringr::str_remove(.x,"panel_"))
  
  design_summary = bind_rows(
    comp_design,
    panel_design
  )
  
  performance = comp_design %>% mutate(
    false_pos_controls = rbinom(n_comp, n_controls, (1-design_spec)),
    n_controls = n_controls,
    false_neg_diseased = rbinom(n_comp, n_diseased, (1-design_sens)),
    n_diseased = n_diseased,
    spec = beta_dist(q=false_pos_controls, n=n_controls),
    sens = beta_dist(q=false_neg_diseased, n=n_diseased),
  ) %>% select(-design_prev) 
  
  samples = tibble::tibble(
      id = 1:n_samples
    ) %>% 
    dplyr::cross_join(comp_design) %>%
    dplyr::mutate(
      actual = rbinom(n_samples * n_comp, 1, design_prev),
      observed = rbinom(n_samples * n_comp, 1, actual * design_sens + (1-actual) * (1-design_spec))
    ) %>%
    dplyr::select(-starts_with("design"))
    
  panel = samples %>%
    group_by(id) %>%
    summarise( 
      actual = 1-prod(1-actual),
      observed = 1-prod(1-observed),
      .groups = "drop"
    )
  
  comp_summary = samples %>% 
    group_by(test) %>%
    summarise(
      actual_pos = sum(actual),
      observed_pos = sum(observed),
      n_samples = n(),
      .groups = "drop"
    )
  
  panel_summary = panel %>% 
    summarise(
      actual_pos = sum(actual),
      observed_pos = sum(observed),
      n_samples = n(),
      .groups = "drop"
    ) %>% mutate(
      test = "Panel"
    )
  
  summary = bind_rows(
    comp_summary, 
    panel_summary
  ) %>% mutate(
    true_prev = actual_pos/n_samples,
    apparent_prev = observed_pos/n_samples,
  )
  
  return(list(
    design = design_summary,
    samples = samples,
    summary = summary,
    performance = performance
  ))
}