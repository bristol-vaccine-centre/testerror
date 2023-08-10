# title: "Uncertainty propagation in VE studies"
# date: "2023-03-19"
# saving files to vignettes/latex/s2/fig

library(tidyverse)
library(patchwork)
library(rstan)
devtools::load_all("~/Git/ggrrr")
devtools::load_all("~/Git/avoncap")
devtools::load_all()
here::i_am("vignettes/bayesian-model-figures.R")
source(here::here("vignettes/vignette-utils.R"))
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
ggrrr::gg_pedantic()

app=apparent_prevalence_plot(p = c(0.1,0.3),top_left = "")

tmp = crossing(
  prev=seq(0,1,length.out=1001),
  test=0:100
) %>% mutate(
  sens = 0.8,
  spec = 0.95,
  ap = prev*sens+(1-prev)*(1-spec),
  p_test = dbinom(test, prob = ap, size=max(test)),
  test_n = test/max(test)
) %>%
glimpse()

app$layers <- c(geom_tile(aes(x=prev, y=test_n, fill = p_test),data = tmp, inherit.aes = FALSE), app$layers)
app$layers[[2]]$aes_params$colour = "grey40"
app2 = app+
  scale_fill_gradient(high = "#0A0AFF" , low = "white", name="probability", oob=scales::squish)+
  scale_y_continuous(sec.axis = sec_axis( trans=~.*100, name="Binomial count"), breaks=c(0,1),expand = c(0, 0))

save_as(app2,  here::here("vignettes/latex/s2/fig/rogan-gladen"))


# ggplot(crossing(
#   sens=seq(0.5,1,length.out=1001),
#   spec=seq(0.95,1,length.out=1001)
# ) %>% mutate(crit = (1-spec)/(2-spec-sens)),
# aes(x=spec, y=sens))+
#   geom_tile(aes(fill=crit))+
#   scale_fill_gradient(limits=c(0,1), high = "orange" , low = "white", name="critical threshold", oob=scales::squish, guide = "none") +
#   metR::geom_contour2(aes(z=crit,label=sprintf("%1.1f%%",after_stat(level)*100)), 
#                       breaks=c(0.99,0.5,0.2,0.1,0.05,0.01), 
#                       label.placer = metR::label_placer_fraction(frac = 0.3),
#                       skip=0,
#                       label_size = 6/ggplot2:::.pt,
#                       margin = grid::unit(c(2, 2, 2, 2), "pt")
#   )
#   #scale_x_continuous(trans="logit", breaks=c(0.5, 0.75,0.9, 0.95, 0.98, 0.99, 0.9975, 0.9999))

scenario = do_scenario(n_controls=800, n_diseased=26)

# Just use scenario setup to describe sens uncertainty ----
# TODO: need to figure out what te fit looks like.

scenario_config = scenario %>% select(-id,-boot,-actual,-test) %>% unnest(test_pcv_group) %>% distinct()

scenario_config_2 = scenario_config %>% 
  glimpse() %>%
  group_by(across(c(starts_with("group"),starts_with("panel")))) %>%
  summarise(
    panel_size = n(),
    panel_prevalence = testerror::panel_prevalence(test_prevalence),
    panel_sens = testerror::panel_sens(test_prevalence, test_sens, test_spec),
    panel_expected_apparent_prevalence = testerror::panel_prevalence(test_expected_apparent_prevalence),
    panel_spec = testerror::panel_spec(test_spec),
    panel_sens_estimator = testerror::panel_sens_estimator(ap = test_expected_apparent_prevalence,sens = test_sens, spec = test_spec),
    panel_sens_samples = list(testerror::uncertain_panel_sens_estimator(
      pos_obs = test_expected_apparent_prevalence*1000, n_obs = 1000,
      false_pos_controls = test_false_pos_controls, n_controls = test_n_controls,
      false_neg_diseased = test_false_neg_diseased, n_diseased = test_n_diseased
    )),
    panel_spec_samples = list(testerror::uncertain_panel_spec(
      false_pos_controls = test_false_pos_controls, n_controls = test_n_controls
    ))
  ) %>%
  group_by(across(c(starts_with("group"),starts_with("panel")))) %>%
  mutate(
    panel_sens_beta = purrr::map(panel_sens_samples, ~ as_tibble(beta_fit(.x))),
    panel_spec_beta = purrr::map(panel_spec_samples, ~ as_tibble(beta_fit(.x))),
    testerror::uncertain_rogan_gladen(
      pos_obs = panel_expected_apparent_prevalence*1000, n_obs = 1000,
      sens = panel_sens_samples,
      spec = panel_spec_samples,
      prefix = "rogan_gladen"
    )
  ) %>%
  unnest(cols = c(panel_sens_beta,panel_spec_beta),names_sep = ".") %>%
  mutate(
    testerror::true_prevalence(
      pos_obs = panel_expected_apparent_prevalence*1000, n_obs = 1000,
      false_pos_controls = panel_spec_beta.shape2, n_controls = panel_spec_beta.conc,
      false_neg_diseased = panel_sens_beta.shape2, n_diseased = panel_sens_beta.conc,
      prefix = "lang_reiczigel"
    )
  ) %>%
  glimpse()

# tmp = tibble(
#   type = "sens",
#   sample = scenario_config_2$panel_sens_samples[[1]]
# ) %>% bind_rows(tibble(
#   type = "spec",
#   sample = scenario_config_2$panel_spec_samples[[1]]
# ))
# 
# ggplot(tmp, aes(x=sample))+geom_density()+facet_wrap(~type)

# ggplot(serotype_prevalence_2, aes(x=false_neg_rate))+geom_density()
# ggplot(serotype_prevalence_2, aes(x=false_pos_rate))+geom_density()

# Compile the multiple tests model ----

stan_model_combined = rstan::stan_model(file = here::here("vignettes/multiple-tests.stan"))

sampling = memoise::memoise(rstan::sampling,cache = cd)

# The simulated data (which is grouped) needs to be fed to stan.
# the data has a patient id (id) and a test id (test_id)
do_bayesian_model = function(d, g = tibble::tibble(), idCol = id, testIdCol = test_id, resultCol = test, sens = 0.8, spec = 0.9975, 
                             #false_pos_controls = (n_controls-2/3)*(1-spec)+1/3,
                             false_pos_controls = (n_controls)*(1-spec),
                             n_controls = 1/(1-spec),
                             #false_neg_diseased = (n_diseased-2/3)*(1-sens)+1/3,
                             false_neg_diseased = (n_diseased)*(1-sens),
                             n_diseased = 1/(1-sens), 
                             controls = NULL,
                             model = stan_model_combined,
                             ...) {
  idCol = rlang::enexpr(idCol)
  testIdCol = rlang::enexpr(testIdCol)
  resultCol = rlang::enexpr(resultCol)
  
  tmp = d %>%
    transmute(
      id = !!idCol,
      test_id = !!testIdCol,
      result = !!resultCol
    ) 
  if (is.factor(tmp$test_id)) {
    tests = levels(forcats::fct_drop(tmp$test_id))
  } else {
    tests = unique(tmp$test_id)
  }
  n_test = length(tests)
  
  tmp = tmp %>%
    pivot_wider(names_from = test_id, values_from = result,values_fill = NA_integer_) %>%
    select(-id) %>% 
    as.matrix()
  
  if (!is.null(controls)) {
    
    y_control = controls %>%
      transmute(
        id = !!idCol,
        test_id = !!testIdCol,
        result = !!resultCol
      ) %>%
      pivot_wider(names_from = test_id, values_from = result,values_fill = NA_integer_) %>%
      select(-id) %>% as.matrix()
    
  } else {
    y_control =  matrix(numeric(), ncol = n_test, nrow = 0, byrow = TRUE)
  }
  
  est_data_combined = list(
    n_test = n_test,
    k_sample = nrow(tmp),
    y_sample = tmp,
    k_control = nrow(y_control),
    y_control = y_control,
    # n_disease = nrow(y_disease),
    # y_disease = y_disease,
    # Priors on sensitivity / specificity
    a_spec = n_controls-false_pos_controls,
    b_spec = false_pos_controls,
    a_sens = n_diseased-false_neg_diseased,
    b_sens = false_neg_diseased
  )
  
  sens_prior = do.call(sprintf, c(list(fmt = "%1.4f [%1.4f \u2013 %1.4f]"), 1-false_neg_diseased/n_diseased, 
                                  qbeta(c(0.025,0.975), shape1 = n_diseased-false_neg_diseased, shape2 = false_neg_diseased)))
  
  spec_prior = do.call(sprintf, c(list(fmt = "%1.4f [%1.4f \u2013 %1.4f]"), 1-false_pos_controls/n_controls, 
                                  qbeta(c(0.025,0.975), shape1 = n_controls-false_pos_controls, shape2 = false_pos_controls)))
  
  
  fit_combined = sampling(
      model,
      data = est_data_combined,
      chains = 4,
      warmup = 1000,          # number of warmup iterations per chain
      iter = 2000,            # total number of iterations per chain
      show_messages = TRUE
    )
  
  summ = summary(fit_combined, pars = c("p","sens","spec"))$summary
  
  summ2 = as_tibble(summ,rownames = "param") %>% 
    mutate(
      test_id = tests[as.integer(stringr::str_extract(param,"[0-9]+"))],
      param = stringr::str_extract(param,"[a-zA-Z]+")
    )
  
  comb_summ = summary(fit_combined, pars = c("p_combined","sens_combined","spec_combined"))$summary
  comb_summ2 = as_tibble(comb_summ,rownames = "param") %>% 
    mutate(
      param = stringr::str_extract(param,"[a-zA-Z]+")
    )
  
  return(
    tibble(
      data = list(d),
      standata = list(est_data_combined),
      stanfit = list(fit_combined),
      component = list(summ2),
      panel = list(comb_summ2),
      tests = list(tests),
      sens_prior = sens_prior,
      spec_prior = spec_prior
    )
  )
}


# Test the multiple tests model ----

bayesian_result = scenario %>% 
  unnest(test_pcv_group) %>%
  #glimpse() %>%
  #filter(panel_id == "PCV20") %>%
  group_by(group, panel_id) %>%
  group_modify(do_bayesian_model, sens=0.8, spec = 0.9975, n_controls=800, n_diseased=26)


do_combine_results = function(bayesian_result) {
  
  grp = bayesian_result %>% groups()
  serotype_tests = bayesian_result %>% 
    mutate(
      prior_false_pos_controls = purrr::map_dbl(standata, ~ .x$b_spec),
      prior_n_controls = purrr::map_dbl(standata, ~ .x$a_spec+.x$b_spec),
      prior_false_neg_diseased = purrr::map_dbl(standata, ~ .x$b_sens),
      prior_n_diseased = purrr::map_dbl(standata, ~ .x$a_sens+.x$b_sens),
      prior_spec = purrr::map_dbl(standata, ~ .x$a_spec/(.x$a_spec+.x$b_spec)),
      prior_sens = purrr::map_dbl(standata, ~ .x$a_sens/(.x$a_sens+.x$b_sens)),
    ) %>% 
    select(!!!grp, data, starts_with("prior")) %>% 
    unnest(data)
  
  disagg_prev = serotype_tests %>%
    group_by(!!!grp, across(c(starts_with("prior"), starts_with("group"),starts_with("test"),-test))) %>%
    summarise(
      actual = sum(actual),
      test = sum(test),
      total = n()
    ) %>%
    group_by(!!!grp, across(c(starts_with("prior"), starts_with("group"),starts_with("test"),-test))) %>%
    mutate(
      actual_prevalence = actual/total,
      prev.est = testerror::rogan_gladen(test/total, test_sens, test_spec),
      testerror::true_prevalence(pos_obs = test, n_obs = total, false_pos_controls = prior_false_pos_controls, n_controls = prior_n_controls, false_neg_diseased = prior_false_neg_diseased, n_diseased = prior_n_diseased,),
      testerror::true_prevalence(pos_obs = test, n_obs = total, false_pos_controls = prior_false_pos_controls, n_controls = prior_n_controls, false_neg_diseased = prior_false_neg_diseased, n_diseased = prior_n_diseased,method = "rogan-gladen", prefix="rogan_gladen")
    )
  
  combined_prev = serotype_tests %>%
    mutate(prevalence = test_prevalence) %>%
    group_by(!!!grp, across(starts_with("group"))) %>%
    estimate_panel_performance_uncertain(test_id) %>%
    mutate(
      prev_est.0.5 = testerror::rogan_gladen(panel_apparent_prevalence, panel_sens_est, panel_spec),
      panel_actual_prevalence = actual/total
    )
  
  comp_tmp = bayesian_result %>% ungroup() %>% select(!!!grp,group,component) %>% unnest(component) %>% filter(param=="p") %>%
    left_join(disagg_prev %>% mutate(test_id = as.character(test_id)), by=join_by(!!!grp, test_id)) %>%
    mutate(binom_ci_2(test,total,"apparent"))
  
  panel_tmp = bayesian_result %>% ungroup() %>% select(!!!grp,group,panel) %>% unnest(panel) %>%  filter(param=="p") %>%
    inner_join(combined_prev,by=join_by(!!!grp)) %>%
    mutate(binom_ci_2(test,total,"apparent"))
  
  return(list(
    components = comp_tmp,
    panels = panel_tmp,
    priors = sprintf("prior:\nsens: %s\nspec: %s", unique(bayesian_result$sens_prior), unique(bayesian_result$spec_prior)),
    params = sprintf("simulation:\nsens: %1.4f\nspec: %1.4f", unique(serotype_tests$test_sens),unique(serotype_tests$test_spec)),
    raw_result = bayesian_result
  ))
}

combined = do_combine_results(bayesian_result)

do_plots = function(combined, show=c("bayes","lang-reiczigel","rogan-gladen")) {

  comp_tmp = combined$components %>% filter(panel_id == "PCV20")
  panel_tmp = combined$panels %>% filter(panel_id == "PCV20")
  
  sep = 0.002
  
  # b = "bayes" %in% show
  # l = "lang-reiczigel" %in% show
  
  offset = seq(-(length(show)-1)*sep/2,(length(show)-1)*sep/2,length.out=length(show))
  
  b = c(offset[which(show=="bayes")],0)[[1]]
  l = c(offset[which(show=="lang-reiczigel")],0)[[1]]
  r = c(offset[which(show=="rogan-gladen")],0)[[1]]
  
  points = bind_rows(
    comp_tmp %>% transmute(x=actual_prevalence+b/2.5, test_id = test_id, type="bayes", lower = `2.5%`, mid = `50%`, upper=`97.5%`),
    comp_tmp %>% transmute(x=actual_prevalence+r/2.5, test_id = test_id, type="rogan-gladen", lower = rogan_gladen.lower, mid = rogan_gladen.median, upper=rogan_gladen.upper),
    comp_tmp %>% transmute(x=actual_prevalence+l/2.5, test_id = test_id, type="lang-reiczigel", lower = prevalence.lower, mid = prevalence.median, upper=prevalence.upper)
  )
  
  points2 = bind_rows(
    panel_tmp %>% transmute(x=panel_actual_prevalence+b, panel_id = panel_id, type="bayes", lower = `2.5%`, mid = `50%`, upper=`97.5%`),
    panel_tmp %>% transmute(x=panel_actual_prevalence+r, panel_id = panel_id, type="rogan-gladen", lower = rogan_gladen.lower, mid = rogan_gladen.median, upper=rogan_gladen.upper),
    panel_tmp %>% transmute(x=panel_actual_prevalence+l, panel_id = panel_id, type="lang-reiczigel", lower = prevalence.lower, mid = prevalence.median, upper=prevalence.upper)
  )
  
  points = points %>% filter(type %in% show)
  points2 = points2 %>% filter(type %in% show)
  
  cilim = binom::binom.confint(seq(0,0.25,length.out=51)*1000,rep(1000,51),methods="wilson")
  
  p1 = ggplot(points, aes(x=x, y= mid, ymin=lower, ymax=upper,colour=type))+
    geom_point(size=0.5)+
    geom_errorbar(width = sep/2.5, size=0.25)+
    coord_fixed(xlim=c(0,0.1),ylim=c(0,0.1))+
    #geom_point(aes(y=prev.est), colour="magenta", size=1)+
    #geom_errorbar(aes(ymin=apparent.0.025, ymax=apparent.0.975), colour="red",position=position_nudge(x=-0.001/2.5), width = 0, alpha=0.4)+
    # geom_point(aes(y=prevalence.median),position=position_nudge(x=-sep/2.5),  size=0.5)+
    # geom_errorbar(aes(ymin=prevalence.lower, ymax=prevalence.upper,linetype="lang-reiczigel"),position=position_nudge(x=-sep/2.5), width = 0, alpha=0.4)+
    annotate("text", x=0.005, y=0.095, label = combined$params, vjust="inward", hjust="inward")+
    annotate("text", x=0.095, y=0.005, label = combined$priors, vjust="inward", hjust="inward")+
    geom_point(aes(x=actual_prevalence, y=apparent.0.5), data=comp_tmp, colour="red",shape=4, size=2, inherit.aes = FALSE)+
    geom_abline(colour="#8080FF")+
    # geom_line(aes(x=x/n,y=upper),data = cilim, colour="#8080FF",linetype="dotted", inherit.aes = FALSE)+
    # geom_line(aes(x=x/n,y=lower),data = cilim, colour="#8080FF",linetype="dotted", inherit.aes = FALSE)+
    xlab("simulation prevalence")+
    ylab("estimated prevalence")+
    facet_wrap(~"components")+
    guides(colour= guide_legend(title=element_blank()))+
    scale_color_grey(start = 0, end = 0.6)+
    theme(legend.position = "bottom")
  
  p2 = ggplot(points2, aes(x=x, y= mid, ymin=lower, ymax=upper,colour=type))+
    geom_point(size=0.5)+
    geom_errorbar(width = sep, size=0.25)+
    coord_fixed(xlim=c(0,0.25),ylim=c(0,0.25))+
    #geom_point(aes(y=prev_est.0.5), colour="magenta", size=1)+
    #geom_errorbar(aes(ymin=apparent.0.025, ymax=apparent.0.975), colour="red",position=position_nudge(x=-0.001), width = 0, alpha=0.4)+
    # geom_point(aes(y=prevalence.median),position=position_nudge(x=-sep),  size=0.5)+
    # geom_errorbar(aes(ymin=prevalence.lower, ymax=prevalence.upper, linetype="lang-reiczigel"),position=position_nudge(x=-sep), width = 0, alpha=0.4)+
    annotate("text", x=0.01, y=0.24, label = combined$params, vjust="inward", hjust="inward")+
    annotate("text", x=0.24, y=0.01, label = combined$priors, vjust="inward", hjust="inward")+
    geom_point(aes(x=panel_actual_prevalence, y=apparent.0.5), data=panel_tmp, colour="red",shape=4,  size=2, inherit.aes = FALSE)+
    geom_abline(colour="#8080FF")+
    # geom_line(aes(x=x/n,y=upper),data = cilim, colour="#8080FF",linetype="dotted", inherit.aes = FALSE)+
    # geom_line(aes(x=x/n,y=lower),data = cilim, colour="#8080FF",linetype="dotted", inherit.aes = FALSE)+
    xlab("simulation prevalence")+
    ylab("estimated prevalence")+
    facet_wrap(~"panel")+
    guides(colour= guide_legend(title=element_blank()))+
    scale_color_grey(start = 0, end = 0.6)+
    theme(legend.position = "bottom")
  
  p3 = p1+p2+patchwork::plot_layout(guides="collect") & theme(legend.position = 'bottom',legend.justification = "center")
  # browser()
  return(list(p1,p2,p3))
}

p = do_plots(combined) #, show="bayes")
p

save_as(p[[3]], here::here("vignettes/latex/s2/fig/simulation_result_sens_80_80"),size = std_size$third)

do_bar_plot = function(combined, prev_filter = 0.1, show=c("bayes","lang-reiczigel","rogan-gladen")) {
  
  tests = # combined$raw_result %>% filter(panel_id == "PCV20") %>% pull(tests) %>% `[[`(1)
    avoncap::serotype_data$xr %>% filter(order==1) %>% pull(serotype)
  
  comp_tmp = combined$components %>% 
    filter(group_prevalence == prev_filter) %>%
    mutate(test_id = factor(test_id, tests)) %>% filter(panel_id == "PCV20" & !is.na(test_id))
  panel_tmp = combined$panels %>% filter(group_prevalence == prev_filter) %>%
    mutate(panel_id = factor(panel_id, c("PCV7","PCV13","PCV15","PCV20")))
  
  points = bind_rows(
    comp_tmp %>% transmute(test_id = test_id, type="bayes", lower = `2.5%`*100, mid = `50%`*100, upper=`97.5%`*100),
    comp_tmp %>% transmute(test_id = test_id, type="lang-reiczigel", lower = prevalence.lower*100, mid = prevalence.median*100, upper=prevalence.upper*100),
    comp_tmp %>% transmute(test_id = test_id, type="rogan-gladen", lower = rogan_gladen.lower*100, mid = rogan_gladen.median*100, upper=rogan_gladen.upper*100)
  )
  
  points2 = bind_rows(
    panel_tmp %>% transmute(panel_id = panel_id, type="bayes", lower = `2.5%`*100, mid = `50%`*100, upper=`97.5%`*100),
    panel_tmp %>% transmute(panel_id = panel_id, type="lang-reiczigel", lower = prevalence.lower*100, mid = prevalence.median*100, upper=prevalence.upper*100),
    panel_tmp %>% transmute(panel_id = panel_id, type="rogan-gladen", lower = rogan_gladen.lower*100, mid = rogan_gladen.median*100, upper=rogan_gladen.upper*100)
  )
  
  points = points %>% filter(type %in% show)
  points2 = points2 %>% filter(type %in% show)
  
  p3 = ggplot(comp_tmp, aes(x=test_id, y=test/total*100))+
    geom_bar(stat="identity",fill="grey90", colour="black",size=0.05)+
    geom_errorbar(aes(ymin=test/total*100,ymax=test/total*100),colour="red")+
    geom_errorbar(aes(ymin=actual/total*100,ymax=actual/total*100),colour="#8080FF")+
    geom_point(aes(x=test_id, y=mid, colour=type),data=points,inherit.aes = FALSE,size=0.5,position = position_dodge(width=0.4))+
    geom_errorbar(aes(x=test_id, ymin=lower, ymax=upper, colour=type),data=points,inherit.aes = FALSE,width=0,position = position_dodge(width=0.4))+
    
    ggpp::annotate("text_npc", npcx = 0.05, npcy = 0.95, label = sprintf("prevalence: %1.2g\n\n%s \n\n%s", prev_filter, combined$params, combined$priors), size=6/ggplot2::.pt )+
    xlab(NULL)+ylab("prevalence (%)")+coord_cartesian(ylim=c(0,15))+
    guides(colour= guide_legend(title=element_blank()))+
    scale_color_grey(start = 0, end = 0.6)
  
  p4 = ggplot(panel_tmp, 
              aes(x=panel_id, y=test/total*100))+
    geom_bar(stat="identity",fill="grey90", colour="black",size=0.05)+
    geom_errorbar(aes(ymin=test/total*100,ymax=test/total*100),colour="red")+
    geom_errorbar(aes(ymin=actual/total*100,ymax=actual/total*100),colour="#8080FF")+
    geom_point(aes(x=panel_id, y=mid, colour=type),data=points2,inherit.aes = FALSE,size=0.5,position = position_dodge(width=0.4))+
    geom_errorbar(aes(x=panel_id, ymin=lower, ymax=upper, colour=type),data=points2,inherit.aes = FALSE,width=0,position = position_dodge(width=0.4))+
    
    xlab(NULL)+ylab("prevalence (%)")+coord_cartesian(ylim=c(0,15))+
    theme(axis.text.x.bottom = element_text(angle=45,vjust = 1,hjust=1))+
    guides(colour= guide_legend(title=element_blank()))+
    scale_color_grey(start = 0, end = 0.6)
  
  p5 = p3+p4+no_y()+patchwork::plot_layout(nrow=1,widths = c(5,1),guides="collect")&theme(legend.position = "bottom")
  
  return(list(p3,p4,p5))
}

p2 = do_bar_plot(combined)

save_as(p2[[3]], here::here("vignettes/latex/s2/fig/simulation_result_prev_10"),size = std_size$third)

# Other scenarios ----

# intro in supplementary
p5 = do_bar_plot(combined, show=NULL)
save_as(p5[[3]], here::here("vignettes/latex/s2/fig/simulation_setup_prev_10"),size = std_size$third)

# main paper
p5 = do_bar_plot(combined, show="bayes")
save_as(p5[[3]]&guides(colour=guide_none()), here::here("vignettes/latex/main/fig/simulation_result_bayes"),size = std_size$third)


## Sens = 60 ----

serotype_tests_60 = do_scenario(
  spec = 0.9975,n_controls = 800,
  sens = 0.6,n_diseased = 260
)

model_60 = serotype_tests_60 %>% 
  unnest(test_pcv_group) %>%
  filter(panel_id == "PCV20") %>%
  group_by(group, panel_id) %>%
  group_modify(do_bayesian_model, sens=0.6, spec = 0.9975, n_controls=800, n_diseased=26)

p_60 = model_60 %>% 
  do_combine_results() %>%
  do_plots()

# p_60
# save_as(p_60[[3]], here::here("vignettes/latex/s2/fig/simulation_result_sens_80_80"),size = std_size$third)

## Sens = 75 ----

serotype_tests_75 = do_scenario(
  spec = 0.9975,n_controls = 800,
  sens = 0.75,n_diseased = 260
)

model_75 = serotype_tests_75 %>% 
  unnest(test_pcv_group) %>%
  filter(panel_id == "PCV20") %>%
  group_by(group, panel_id) %>%
  group_modify(do_bayesian_model, sens=0.75, spec = 0.9975, n_controls=800, n_diseased=26)

p_75 = model_75 %>% 
  do_combine_results() %>%
  do_plots()

# p_75

## Sens = 90 ----

serotype_tests_90 = do_scenario(
  spec = 0.9975,n_controls = 800,
  sens = 0.90,n_diseased = 260
)

model_90 = serotype_tests_90 %>% 
  unnest(test_pcv_group) %>%
  filter(panel_id == "PCV20") %>%
  group_by(group, panel_id) %>%
  group_modify(do_bayesian_model, sens=0.90, spec = 0.9975, n_controls=800, n_diseased=26)

p_90 = model_90 %>% 
  do_combine_results() %>%
  do_plots()

# p_90

p_90_3 = model_90 %>% 
  do_combine_results() %>%
  do_bar_plot()

# p_90_3

## Combination of alternative scenario plots ----

p = p_60[[1]]+ no_x()+ #facet_grid("sens: 0.60"~"components")+
  p_60[[2]]+ no_x()+ #facet_grid("sens: 0.60"~"panel")+
  p_75[[1]]+ no_x()+ #facet_grid("sens: 0.75"~"components")+
  p_75[[2]]+ no_x()+ #facet_grid("sens: 0.75"~"panel")+
  p_90[[1]]+ #facet_grid("sens: 0.90"~"components")+
  p_90[[2]]+ #facet_grid("sens: 0.90"~"panel")+
  plot_layout(ncol=2, guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")

save_as(p, here::here("vignettes/latex/s2/fig/simulation-result-same-sens"),size = std_size$full)


# p2 = 
#   p_60[[2]]+ #facet_grid("sens: 0.60"~"panel")+
#   p_75[[2]]+ no_y()+ #facet_grid("sens: 0.75"~"panel")+
#   p_90[[2]]+ no_y()+ #facet_grid("sens: 0.90"~"panel")+
#   plot_layout(ncol=3)+plot_annotation(tag_levels = "A")
# 
# p6 = (p2)+p[[3]]+plot_layout(design = "
# AAAAAAABBBBBBBCCCCCCC
# AAAAAAABBBBBBBCCCCCCC
# AAAAAAABBBBBBBCCCCCCC
# AAAAAAABBBBBBBCCCCCCC
# AAAAAAABBBBBBBCCCCCCC
# AAAAAAABBBBBBBCCCCCCC
# AAAAAAABBBBBBBCCCCCCC
# DDDDDDDDDDDDDDDDDDDEE
# DDDDDDDDDDDDDDDDDDDEE
# DDDDDDDDDDDDDDDDDDDEE
# DDDDDDDDDDDDDDDDDDDEE
# DDDDDDDDDDDDDDDDDDDEE
# DDDDDDDDDDDDDDDDDDDEE
# ")
#                            
# p6                           
# 
# save_as(p2, here::here("vignettes/latex/main/figs/bayesian-sim"),size = std_size$third)
# 
# save_as(p_90[[3]], here::here("vignettes/latex/main/figs/sim-90-example"),size = std_size$third)

# Sensitivity mismatch ----

serotype_tests_60_90 = do_scenario(
  spec = 0.9975,n_controls = 800,
  sens = 0.60,n_diseased = 260
)

model_60_90 = serotype_tests_60_90 %>% 
  unnest(test_pcv_group) %>%
  filter(panel_id == "PCV20") %>%
  group_by(group, panel_id) %>%
  group_modify(do_bayesian_model, sens=0.90, spec = 0.9975, n_controls=800, n_diseased=26)

p_60_90 = model_60_90 %>%
  do_combine_results() %>%
  do_plots()

# p_60_90

serotype_tests_90_60 = do_scenario(
  spec = 0.9975,n_controls = 800,
  sens = 0.90,n_diseased = 260
)

model_90_60 = serotype_tests_90_60 %>% 
  unnest(test_pcv_group) %>%
  filter(panel_id == "PCV20") %>%
  group_by(group, panel_id) %>%
  group_modify(do_bayesian_model, sens=0.60, spec = 0.9975, n_controls=800, n_diseased=26)

p_90_60 = model_90_60 %>%
  do_combine_results() %>%
  do_plots() #show="bayes")

# p_90_60

# p = p_60_90[[2]]+ #facet_grid("sens: 0.60"~"panel")+
#   p_90_60[[2]]+no_y()+ #facet_grid("sens: 0.90"~"panel")+
#   plot_layout(ncol=2, guides="collect")+plot_annotation(tag_levels = "A")

p = p_60_90[[1]]+p_60_90[[2]]+ #facet_grid("sens: 0.60"~"panel")+
  p_90_60[[1]]+p_90_60[[2]]+ #facet_grid("sens: 0.90"~"panel")+
  plot_layout(ncol=2,guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")


save_as(p, here::here("vignettes/latex/s2/fig/bayesian-sim-mismatch-sens"),size = std_size$two_third)

# Specificity mismatch ----

serotype_tests_9975_99 = do_scenario(
  spec = 0.9975,n_controls = 800,
  sens = 0.80,n_diseased = 260
)

model_9975_99 = serotype_tests_9975_99 %>% 
  unnest(test_pcv_group) %>%
  filter(panel_id == "PCV20") %>%
  group_by(group, panel_id) %>%
  group_modify(do_bayesian_model, sens=0.80, spec = 0.99, n_controls=100, n_diseased=26)

p_9975_99 = model_9975_99 %>%
  do_combine_results() %>%
  do_plots()



serotype_tests_99_9975 = do_scenario(
  spec = 0.99,n_controls = 800,
  sens = 0.80,n_diseased = 260
)

model_99_9975 = serotype_tests_99_9975 %>% 
  unnest(test_pcv_group) %>%
  filter(panel_id == "PCV20") %>%
  group_by(group, panel_id) %>%
  group_modify(do_bayesian_model, sens=0.80, spec = 0.9975, n_controls=100, n_diseased=26)

p_99_9975 = model_99_9975 %>%
  do_combine_results() %>%
  do_plots()

p = p_9975_99[[1]]+p_9975_99[[2]]+ #facet_grid("sens: 0.60"~"panel")+
  p_99_9975[[1]]+p_99_9975[[2]]+ #facet_grid("sens: 0.90"~"panel")+
  plot_layout(ncol=2,guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")

save_as(p, here::here("vignettes/latex/s2/fig/bayesian-sim-mismatch-spec"),size = std_size$two_third)


p = p_60_90[[2]]+
  p_90_60[[2]]+p_9975_99[[2]]+p_99_9975[[2]]+ #facet_grid("sens: 0.90"~"panel")+
  plot_layout(ncol=2,guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")&coord_cartesian(xlim=c(0,0.25),ylim=c(0,0.4))

save_as(p, here::here("vignettes/latex/s2/fig/bayesian-sim-mismatch"),size = std_size$two_third)

