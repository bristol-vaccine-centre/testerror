
## Panel test example ----

# a set of panels with different prevalences.
# configurable with component distribution, overall prevalence, and number of 
# components
multi_panel_example = function(
    panel_prev = seq(0.025,0.2, length.out=8), 
    n_comp = 20,
    comp_dist = stats::rpois(n_comp,10) * c(1,stats::rpois(n_comp-1,2)),
    ...) {
  tmp1 = tibble(
    scenario_name = sprintf("P: %1.2f%%",panel_prev*100),
    scenario_prev = panel_prev
  ) 
  tmp2=NULL
  for (p in panel_prev) {
    tmp3 = panel_example(
      panel_prev = p, 
      n_comp = n_comp,
      comp_dist = comp_dist, 
      comp_prev = distribute(comp_dist, panel_prev),
      ...)
    tmp4 = as_tibble(purrr:::map(tmp3, ~ list(.x)))
    tmp2 = bind_rows(tmp2, tmp4)
  }
  return(bind_cols(tmp1,tmp2))
}

# generate a simple dataset and summarise it.
panel_example = function(
    panel_prev = 0.1,
    n_comp = floor(stats::runif(1, 5, 11)), 
    comp_dist = stats::rpois(n_comp,10) * c(1,stats::rpois(n_comp-1,2)),
    comp_spec = stats::rbeta(n_comp, 995, 5),
    comp_sens = stats::rbeta(n_comp, 40, 10),
    comp_test = factor(sprintf("Component %s",LETTERS[1:n_comp])),
    comp_prev = distribute(comp_dist, panel_prev),
    panel_name = "Panel",
    panel_spec = NULL,
    n_samples = 1000,
    n_controls = 800,
    n_diseased = 30,
    n_boots = 1,
    seed = 1001,
    exact = FALSE,
    exact_controls = exact,
    exact_samples = exact
) {
  
  set.seed(seed)
  
  if (is.null(panel_spec)) {
    panel_spec = tibble::tibble(panel_name = panel_name, comp_test = comp_test)
  } else if (is.data.frame(panel_spec)) {
    # leave as is
  } else {
    panel_spec = dplyr::bind_rows(lapply(seq_along(panel_spec), 
                     function(i) tibble::tibble(
                        panel_name = names(panel_spec)[[i]],
                        comp_test = panel_spec[[i]]
                      )
    ))
  }
  panel_spec = panel_spec %>% group_by(panel_name) %>% mutate(n_components=n())
  
  tmp = test_example(
    prev = comp_prev,
    spec = comp_spec,
    sens = comp_sens,
    name = comp_test,
    n_samples = n_samples,
    n_controls = n_controls,
    n_diseased = n_diseased,
    n_boots = n_boots,
    seed = seed,
    exact_controls = exact_controls,
    exact_samples = exact_samples
  )
  
  # the theoretical combination of components into panel sens, spec and prev
  # this is not grouped by boot
  panel_design = tmp$design %>%
    dplyr::mutate(test = as.character(test)) %>%
    dplyr::inner_join(panel_spec %>% mutate(comp_test = as.character(comp_test)), by=c("test"="comp_test")) %>%
    dplyr::group_by(panel_name, n_components) %>%
    dplyr::summarise(
      tmp_design_prev = panel_prevalence(design_prev),
      tmp_design_spec = panel_spec(design_spec),
      tmp_design_sens = panel_sens(design_prev,design_sens,design_spec)
    ) %>% 
    dplyr::rename(test = panel_name) %>%
    dplyr::rename_with(.cols = tidyselect::starts_with("tmp_"), .fn = ~ stringr::str_remove(.x,"tmp_"))
  
  # Panel and component design
  tmp$design = dplyr::bind_rows(
    tmp$design %>% mutate(n_components=1), 
    panel_design
  )
  
  # calculate a panel test result (at level of individual)
  panel = tmp$samples %>%
    dplyr::mutate(test = as.character(test)) %>%
    dplyr::inner_join(panel_spec %>% mutate(comp_test = as.character(comp_test)), by=c("test"="comp_test")) %>%
    dplyr::group_by(panel_name, n_components,boot,id) %>%
    dplyr::summarise( 
      actual = 1-prod(1-actual),
      observed = 1-prod(1-observed),
      .groups = "drop"
    )
  
  # summarise the panel test results into panel test counts of positives
  panel_summary = panel %>% 
    dplyr::group_by(test = panel_name, n_components,boot) %>%
    dplyr::summarise(
      actual_pos = sum(actual),
      observed_pos = sum(observed),
      n_samples = dplyr::n(),
      .groups = "drop"
    ) %>% dplyr::mutate(
      true_prev = actual_pos/n_samples,
      apparent_prev = observed_pos/n_samples,
    )
  
  # combine simulation actual and observed counts to calculate
  # true and apparent prevalence.
  tmp$summary = dplyr::bind_rows(
    tmp$summary %>% mutate(n_components=1), 
    panel_summary
  )
  
  return(tmp)
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
    exact_controls = exact,
    exact_samples = exact,
    seed = 1001
) {
  
  if (n_boots > 1 && exact_controls && exact_samples) 
    stop("more than one bootstrap with no randomness requested")
  
  set.seed(seed)
  
  n = pkgutils::recycle(prev,spec,sens,n_controls,n_diseased)
  
  # the component design parameters including test sens, spec and prev
  comp_design = tibble::tibble(
    design_prev = prev,
    design_spec = spec,
    design_sens = sens
  ) %>% dplyr::mutate(
    test = forcats::as_factor(name)
  )
  
  if (anyDuplicated(comp_design$test)) stop("test names must be unique")
  
  # sample from the design to create a disease positive and a disease
  # negative control group sample. This simultates the outcome of controls
  # for each component
  performance = comp_design %>% 
    cross_join(tibble(boot = 1:n_boots)) %>%
    either_or(
    exact_controls,
    if_true = ~ .x %>% dplyr::mutate(
      false_pos_controls = round(n_controls * (1-design_spec)),
      false_neg_diseased = round(n_diseased * (1-design_sens)),
    ), 
    if_false = ~ .x %>% dplyr::mutate(
      false_pos_controls = stats::rbinom(dplyr::n(), n_controls, (1-design_spec)),
      false_neg_diseased = stats::rbinom(dplyr::n(), n_diseased, (1-design_sens)),
    )
  ) %>% dplyr::mutate(
    n_controls = rep(n_controls,n_boots),
    n_diseased = rep(n_diseased,n_boots),
    spec = beta_dist(shape1=n_controls-false_pos_controls+1, shape2=false_pos_controls+1),
    sens = beta_dist(shape1=n_diseased-false_neg_diseased+1, shape2=false_neg_diseased+1),
  ) %>% dplyr::select(-design_prev) 
  
  # sample from the design to create a set of test results for each component
  # test with an actual value (defined by prevalence) and an observed (defined 
  # by actual, component sens and spec)
  samples = tidyr::crossing(
    id = 1:n_samples,
    boot = 1:n_boots
  ) %>% 
    dplyr::cross_join(comp_design) %>%
    either_or(
      exact_samples,
      if_true = ~ .x %>%
        group_by(boot,test) %>%
        mutate(actual = rfixed(1, dplyr::n(), design_prev)) %>%
        group_by(boot,test,actual) %>%
        mutate(
          observed = rfixed(1, dplyr::n(), dplyr::if_else(actual==1, design_sens, 1-design_spec))
        ) %>% 
        ungroup(), 
      if_false = ~ .x %>% dplyr::mutate(
        actual = stats::rbinom(n_boots * n_samples * n, 1, design_prev),
        observed = stats::rbinom(n_boots * n_samples * n, 1, actual * design_sens + (1-actual) * (1-design_spec))
      )
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
demo_bar_plot_base = function(
    summary = interfacer::iface(
      test = character ~ "the test",
      apparent_prev = double ~ "observed test positive rate",
      true_prev = double ~ "true positive rate",
      n_samples = integer ~ "the overall number of patients tested"
    ), 
    ...
) {
  prediction = interfacer::ivalidate(summary, ...)
  prediction = prediction %>% dplyr::mutate(
    test = forcats::as_factor(test)
  )
  tmp = prediction %>% dplyr::select(test,apparent_prev) %>% dplyr::distinct()
  n_samples = unique(prediction$n_samples)
  ggplot2::ggplot(prediction)+
    ggplot2::geom_bar(ggplot2::aes(x=test,y=apparent_prev*100), data = tmp, stat="identity", fill="grey80", colour=NA,width=0.8)+
    ggplot2::geom_errorbar(ggplot2::aes(x=test,y=apparent_prev*100,ymin=apparent_prev*100,ymax=apparent_prev*100), colour="red",width=0.8)+
    ggplot2::geom_errorbar(ggplot2::aes(x=test,y=true_prev*100,ymin=true_prev*100,ymax=true_prev*100), colour="blue",width=0.8)+
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(trans = ~ ./100*n_samples, name="counts"),name = "prevalence (%)")
}

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
  

