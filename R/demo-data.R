# Internal functions used in some vignettes and examples

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
  pos = round(n*prev)
  neg = round(n-pos)
  s = c(rep(0,neg),rep(1,pos))
  lapply(1:boots,  function(...) sample(s)) %>% unlist()
}

rfixed_mnom = function(boots, n, dist, prev) {
  dist = dist/sum(dist)
  dist = dist*prev
  tmp = unlist(sapply(1:length(dist), function(x) rep(x,round(dist[x]*n))))
  tmp = c(tmp, rep(0,n-length(tmp)))
  lapply(1:boots,  function(...) sample(tmp)) %>% unlist()
}

## Panel test example ----

# generate a simple dataset and summarise it.
panel_example = function(
    panel_prev = 0.1,
    n_comp = floor(stats::runif(1, 5, 11)), 
    comp_dist = stats::rpois(n_comp,10) * c(1,stats::rpois(n_comp-1,2)),
    comp_spec = stats::rbeta(n_comp, 995, 5),
    comp_sens = stats::rbeta(n_comp, 40, 10),
    n_samples = 1000,
    n_controls = 800,
    n_diseased = 30,
    seed = 1001,
    exact = TRUE
) {
  
  set.seed(seed)
  comp_test = factor(sprintf("Component %s",LETTERS[1:n_comp]))
  comp_prev = distribute(comp_dist, panel_prev)
  
  # the component design parameters includiong test sens, spec and prev
  comp_design = tibble::tibble(
    test = comp_test,
    design_prev = comp_prev,
    design_spec = comp_spec,
    design_sens = comp_sens
  ) 
  
  # the theoretical combination of components into panel sens, spec and prev
  panel_design = comp_design %>%
    dplyr::summarise(
      panel_design_prev = panel_prevalence(design_prev),
      panel_design_spec = panel_spec(design_spec),
      panel_design_sens = panel_sens(design_prev,design_sens,design_spec)
    ) %>% dplyr::mutate(test = "Panel") %>%
    dplyr::rename_with(~ stringr::str_remove(.x,"panel_"))
  
  # Panel and component design
  design_summary = dplyr::bind_rows(
    comp_design,
    panel_design
  )
  
  # sample from the design to create a disease positive and a disease
  # negative control group sample. This simultates the outcome of controls
  # for each component
  performance = comp_design %>% either_or(
    exact,
    if_true = ~ .x %>% dplyr::mutate(
      false_pos_controls = round(n_controls * (1-design_spec)),
      n_controls = n_controls,
      false_neg_diseased = round(n_diseased * (1-design_sens)),
      n_diseased = n_diseased
    ), 
    if_false = ~ .x %>% dplyr::mutate(
      false_pos_controls = stats::rbinom(n_comp, n_controls, (1-design_spec)),
      n_controls = n_controls,
      false_neg_diseased = stats::rbinom(n_comp, n_diseased, (1-design_sens)),
      n_diseased = n_diseased
    )
  ) %>% dplyr::mutate(
    # posterior beta distribtuion assuming a uniform prior 
    spec = beta_dist(p=n_controls-false_pos_controls+1, q=false_pos_controls+1),
    sens = beta_dist(p=n_diseased-false_neg_diseased+1, q=false_neg_diseased+1),
  ) %>% dplyr::select(-design_prev) 
  
  # sample from the design to create a set of test results for each component
  # test with an actual value (defined by prevalence) and an observed (defined 
  # by actual, component sens and spec)
  samples = tibble::tibble(
      id = 1:n_samples
    ) %>% 
    dplyr::cross_join(comp_design) %>%
    dplyr::mutate(
      actual = stats::rbinom(n_samples * n_comp, 1, design_prev),
      observed = stats::rbinom(n_samples * n_comp, 1, actual * design_sens + (1-actual) * (1-design_spec))
    ) %>%
    dplyr::select(-tidyselect::starts_with("design"))
    
  # calculate a panel test result (at level of individual)
  panel = samples %>%
    dplyr::group_by(id) %>%
    dplyr::summarise( 
      actual = 1-prod(1-actual),
      observed = 1-prod(1-observed),
      .groups = "drop"
    )
  
  # summarise the samples into per component test counts of positives
  # this is simulation actual and obseved counts
  comp_summary = samples %>% 
    dplyr::group_by(test) %>%
    dplyr::summarise(
      actual_pos = sum(actual),
      observed_pos = sum(observed),
      n_samples = dplyr::n(),
      .groups = "drop"
    )
  
  # summarise the panel test results into panel test counts of positives
  panel_summary = panel %>% 
    dplyr::summarise(
      actual_pos = sum(actual),
      observed_pos = sum(observed),
      n_samples = dplyr::n(),
      .groups = "drop"
    ) %>% dplyr::mutate(
      test = "Panel"
    )
  
  # combine simulation actual and observed counts to calculate
  # true and apparent prevalence.
  summary = dplyr::bind_rows(
    comp_summary, 
    panel_summary
  ) %>% dplyr::mutate(
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

# Simple test example ----

# generate a simple dataset and summarise it.
test_example = function(
    prev = seq(0,0.2,length.out=11),
    spec = 0.95,
    sens = 0.8,
    n_samples = 1000,
    n_controls = 800,
    n_diseased = 30,
    n_boots = 1,
    name = sprintf("%1.1f%%",prev*100),
    exact = FALSE,
    seed = 1001
) {
  
  set.seed(seed)
  name = rlang::enexpr(name)
  
  n = .recycle(prev,spec,sens,n_controls,n_diseased)
  
  # the component design parameters including test sens, spec and prev
  comp_design = tibble::tibble(
    design_prev = prev,
    design_spec = spec,
    design_sens = sens
  ) %>% dplyr::mutate(
    test = forcats::as_factor(!!name)
  )
  
  if (anyDuplicated(comp_design$test)) stop("test names must be unique")
  
  # sample from the design to create a disease positive and a disease
  # negative control group sample. This simultates the outcome of controls
  # for each component
  performance = comp_design %>% either_or(
    exact,
    if_true = ~ .x %>% dplyr::mutate(
      false_pos_controls = round(n_controls * (1-design_spec)),
      false_neg_diseased = round(n_diseased * (1-design_sens)),
    ), 
    if_false = ~ .x %>% dplyr::mutate(
      false_pos_controls = stats::rbinom(n, n_controls, (1-design_spec)),
      false_neg_diseased = stats::rbinom(n, n_diseased, (1-design_sens)),
    )
  ) %>% dplyr::mutate(
    n_controls = n_controls,
    n_diseased = n_diseased,
    spec = beta_dist(p=n_controls-false_pos_controls+1, q=false_pos_controls+1),
    sens = beta_dist(p=n_diseased-false_neg_diseased+1, q=false_neg_diseased+1),
  ) %>% dplyr::select(-design_prev) 
  
  # sample from the design to create a set of test results for each component
  # test with an actual value (defined by prevalence) and an observed (defined 
  # by actual, component sens and spec)
  samples = tidyr::crossing(
    boot = 1:n_boots,
    id = 1:n_samples
  ) %>% 
    dplyr::cross_join(comp_design) %>%
    
    dplyr::mutate(
      actual = ( #dplyr::if_else(exact,
          #rfixed2(n_boots, n_samples * n, design_prev),
          stats::rbinom(n_boots * n_samples * n, 1, design_prev)
      ),
      observed = stats::rbinom(n_boots * n_samples * n, 1, actual * design_sens + (1-actual) * (1-design_spec))
    ) %>%
    dplyr::select(-tidyselect::starts_with("design"))
  
  # summarise the samples into per component test counts of positives
  # this is simulation actual and obseved counts
  comp_summary = samples %>% 
    dplyr::group_by(boot,test) %>%
    dplyr::summarise(
      actual_pos = sum(actual),
      observed_pos = sum(observed),
      n_samples = dplyr::n(),
      .groups = "drop"
    ) %>% dplyr::mutate(
      true_prev = actual_pos/n_samples,
      apparent_prev = observed_pos/n_samples,
    )
  
  return(list(
    design = comp_design,
    samples = samples,
    summary = comp_summary,
    performance = performance
  ))
}



# Plots ----

# plot true and predicted test results for a single panel.
demo_bar_plot = function(
    prediction = interfacer::iface(
      apparent_prev = double ~ "observed test positive rate",
      true_prev = double ~ "true positive rate",
      n_samples = integer ~ "the overall number of patients tested",
      testerror:::.output_data
    ), 
    ...
) {
  prediction = interfacer::ivalidate(prediction, ...)
  prediction = prediction %>% dplyr::mutate(
    test = forcats::as_factor(test)
  )
  tmp = prediction %>% dplyr::select(test,apparent_prev) %>% dplyr::distinct()
  n_samples = unique(prediction$n_samples)
  ggplot2::ggplot(prediction)+
    ggplot2::geom_bar(ggplot2::aes(x=test,y=apparent_prev), data = tmp, stat="identity", fill="grey80", colour=NA,width=0.8)+
    ggplot2::geom_errorbar(ggplot2::aes(x=test,y=apparent_prev,ymin=apparent_prev,ymax=apparent_prev), colour="red",width=0.8)+
    ggplot2::geom_errorbar(ggplot2::aes(x=test,y=true_prev,ymin=true_prev,ymax=true_prev), colour="blue",width=0.8)+
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = ~ .*n_samples, name="counts"),name = "proportion")+
    ggplot2::geom_point(ggplot2::aes(x=test, y=prevalence.median, colour = prevalence.method), size=1,position = ggplot2::position_dodge(width=0.2))+
    ggplot2::geom_errorbar(ggplot2::aes(x=test, y=prevalence.median, ymin=prevalence.lower, ymax=prevalence.upper, colour = prevalence.method), width=0.15, position = ggplot2::position_dodge(width=0.2))+
    ggplot2::scale_colour_grey(start=0,end = 0.4,name="")+
    ggplot2::xlab(NULL)
}

demo_qq_plot = function(
    prediction = interfacer::iface(
      apparent_prev = double ~ "observed test positive rate",
      true_prev = double ~ "true positive rate",
      n_samples = integer ~ "the overall number of patients tested",
      testerror:::.output_data
    ), 
    ...
) {
  
  prediction = interfacer::ivalidate(prediction, ...)
  
  show = unique(prediction$prevalence.method)
  mx = ceiling(max(prediction$prevalence.upper)*20)/20
  sep = mx/100
  
  offset_df = tibble::tibble(
    prevalence.method = show,
    offset = seq(-(length(show)-1)*sep/2,(length(show)-1)*sep/2,length.out=length(show))
  )
  
  prediction = prediction %>% 
    dplyr::inner_join(offset_df, by="prevalence.method") %>%
    dplyr::mutate(x = true_prev+offset)
    
    p1 = ggplot2::ggplot(prediction, ggplot2::aes(x=x, y= prevalence.median, ymin=prevalence.lower, ymax=prevalence.upper,colour=prevalence.method))+
      ggplot2::geom_point(size=0.5)+
      ggplot2::geom_errorbar(width = sep, size=0.25)+
      ggplot2::coord_fixed(xlim=c(0,mx),ylim=c(0,mx))+
      #ggplot2::geom_point(ggplot2::aes(y=prev.est), colour="magenta", size=1)+
      #ggplot2::geom_errorbar(ggplot2::aes(ymin=apparent.0.025, ymax=apparent.0.975), colour="red",position=ggplot2::position_nudge(x=-0.001/2.5), width = 0, alpha=0.4)+
      # ggplot2::geom_point(ggplot2::aes(y=prevalence.median),position=ggplot2::position_nudge(x=-sep/2.5),  size=0.5)+
      # ggplot2::geom_errorbar(ggplot2::aes(ymin=prevalence.lower, ymax=prevalence.upper,linetype="lang-reiczigel"),position=ggplot2::position_nudge(x=-sep/2.5), width = 0, alpha=0.4)+
      # ggplot2::annotate("text", x=0.005, y=0.095, label = combined$params, vjust="inward", hjust="inward")+
      # ggplot2::annotate("text", x=0.095, y=0.005, label = combined$priors, vjust="inward", hjust="inward")+
      ggplot2::geom_point(ggplot2::aes(x=true_prev, y=apparent_prev), colour="red",shape=4, size=2, inherit.aes = FALSE)+
      ggplot2::geom_abline(colour="#8080FF")+
      ggplot2::xlab("true prevalence")+
      ggplot2::ylab("estimated prevalence")+
      ggplot2::facet_wrap(~"components")+
      #ggplot2::guides(colour= ggplot2::guide_legend(title=ggplot2::element_blank()))+
      ggplot2::scale_color_grey(start = 0, end = 0.6, name="")+
      ggplot2::theme(legend.position = "bottom",legend.justification = "center")
    
    return(p1)
}
  

