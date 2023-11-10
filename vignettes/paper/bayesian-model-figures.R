# title: "Uncertainty propagation in VE studies"
# date: "2023-03-19"
# saving files to vignettes/latex/s2/fig

library(tidyverse)
library(patchwork)
library(rstan)
devtools::load_all("~/Git/ggrrr")
devtools::load_all("~/Git/avoncap")
devtools::load_all()
here::i_am("vignettes/paper/bayesian-model-figures.R")
source(here::here("vignettes/vignette-utils.R"))
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
ggrrr::gg_pedantic()

## Supp 2 figure 1 ----

app=apparent_prevalence_plot(p = c(0.05,0.3),top_left = "")

# not vectorised
.pos_test_p = function(test, prev, sens, spec, size, ...) {
  sum(
    # TRUE POSITIVES:
    dbinom(test:0, prob = sens*prev, size=size)*
      dbinom(0:test, prob = (1-spec)*(1-prev), size=size))
}

tmp = crossing(
  prev=seq(0,1,length.out=1001),
  test=0:100
) %>% mutate(
  sens = 0.8,
  spec = 0.95,
  size = 100,
  ap = prev*sens+(1-prev)*(1-spec),
  #p_test = dbinom(test, prob = ap, size=max(test)),
  test_n = test/max(test)
) %>% mutate(
  p_test = purrr::pmap_dbl(., .pos_test_p)
)

marks = tibble(
  p_mark = c(0.05, 0.3)
) %>% mutate(
  ap_mark = apparent_prevalence(p_mark, 0.8, 0.95),
  test_mark = round(ap_mark*100),
  binom::binom.confint(test_mark,100,method="wilson")
)

thres = testerror::underestimation_threshold(0.8,0.95)

#p_cols = scales::brewer_pal(palette = "Dark2")(2)
p_cols = c("blue","red")
names(p_cols) = as.character(marks$p_mark)
test_cols = p_cols
names(test_cols) = as.character(marks$test_mark)

marginal_x = tmp %>% inner_join(marks, by = c("test"="test_mark"))
marginal_y = tmp %>% inner_join(marks, by = c("prev"="p_mark"))

app2 = ggplot(tmp, aes(x=prev,y=ap*100))+
  geom_abline(colour="grey40",slope = 100)+
  geom_line()+
  geom_tile(aes(x=prev, y=test, fill = p_test),data = tmp, inherit.aes = FALSE)+
  scale_fill_gradient(trans="sqrt",high = "#0A0AFF" , low = "#FFFFFF00", guide="none",name="probability", oob=scales::squish)+
  scale_y_continuous(sec.axis = sec_axis( trans=~./100, name="apparent prevalence"), expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  theme(axis.text.y.left = element_blank(), axis.title.y.left = element_blank())+
  ggrrr::gg_hide_X_axis()+
  geom_hline(data = marks, mapping=aes(yintercept=ap_mark*100, color = as.character(p_mark)), linetype="dotted")+
  geom_segment(data = marks, mapping=aes(x=p_mark, y=0, xend = p_mark, yend=ap_mark*100, color = as.character(p_mark)), linetype="dashed")+
  geom_point(data = marks, mapping=aes(x=p_mark, y=ap_mark*100, color = as.character(p_mark)), size=1)+
  geom_point(x=thres,y=thres*100,colour="black", size=1)+
  # geom_errorbarh(data= marks, mapping=aes(y=ap_mark*100, x=mean, xmin=lower, xmax=upper, colour = as.character(p_mark)))+
  scale_color_manual(values = p_cols, guide="none")
  

# ggrrr::gg_save_as(app2,  here::here("vignettes/latex/s2/fig/rogan-gladen-v2"), size = std_size$third)

mx = ggplot(marginal_x, aes(x=prev, y=p_test, colour=as.character(test),fill=as.character(test)))+geom_area(position="dodge",alpha=0.1)+
  scale_color_manual(values = test_cols, name = "count", aesthetics = c("colour","fill"))+
  geom_vline(data=marks, mapping=aes(xintercept=p_mark, colour=as.character(test_mark)), linetype = "dashed")+
  geom_rect(data= marks, mapping=aes(ymin=-0.015, ymax=-0.005, xmin=lower, xmax=upper, fill = as.character(test_mark), colour = as.character(test_mark)), alpha=0.2, inherit.aes = FALSE)+
  scale_x_continuous(expand = c(0, 0), name="true prevalence")+scale_y_reverse(expand = c(0.005, 0.005),breaks = NULL)+ggrrr::gg_hide_Y_axis()

  
my = ggplot(marginal_y, aes(x=test, y=p_test, colour=as.character(prev),fill=as.character(prev)))+
  scale_color_manual(values = p_cols, name="prevalence",aesthetics = c("colour","fill"))+
  geom_bar(stat="identity",position=position_identity(),colour=NA,alpha=0.1)+
  geom_step(direction="mid")+
  geom_vline(data=marks, mapping=aes(xintercept=ap_mark*100,colour=as.character(p_mark)), linetype = "dotted")+
  ggplot2::scale_x_continuous(expand = c(0, 0), name="positive count (N=100)")+
  coord_flip()+scale_y_reverse(expand = c(0, 0.005),breaks = NULL)+ggrrr::gg_hide_X_axis()

sf1 = my + app2 + patchwork::guide_area() + mx + plot_layout(design = 
"ABBBB
ABBBB
ABBBB
ABBBB
CDDDD
",guides="collect")

ggrrr::gg_save_as(sf1,  here::here("vignettes/latex/s2/fig/rogan-gladen-v3"), maxWidth = 4, maxHeight = 4)

## Set up scenario / IPD test data ----

# ipd PCV panels:
panels = ipd_distribution() %>% unnest(pcv_group) %>% 
  transmute(panel_name = forcats::as_factor(panel_id), comp_test = pneumo.phe_serotype)

scenario = panel_example(
  n_comp = 20,
  n_samples = 4000,
  n_controls = 800,
  n_diseased = 200,
  comp_spec = 0.9975,
  comp_sens = 0.8,
  comp_prev = ipd_distribution() %>% pull(prev),
  comp_test = ipd_distribution() %>% pull(pneumo.phe_serotype) %>% forcats::as_factor(),
  panel_spec = panels,
  exact_controls = TRUE,
  exact_samples = TRUE
)

caption = c(
scenario$design %>% 
  filter(n_components == 1) %>% 
  transmute(label = sprintf("simulated component test performance: sensitivity %1.2f%%, specificity %1.2f%%",design_sens*100, design_spec*100)) %>% 
  pull(label) %>% unique(),

scenario$design %>% 
  filter(n_components > 1) %>% 
  transmute(label = sprintf("%s \u2013 %1.2f%%",test,design_prev*100)) %>% 
  pull(label) %>% paste0(collapse = "; ") %>% sprintf(fmt = "simulated panel prevalence: %s")#,
) %>% paste0(collapse = "\n")

## Supp 2 figure 2 ----

p1 = demo_bar_plot_base(scenario$summary %>% filter(n_components == 1))+facet_wrap(~"serotypes")+coord_cartesian(ylim = c(0,15))
p2 = demo_bar_plot_base(scenario$summary %>% filter(n_components > 1))+facet_wrap(~"PCV")+coord_cartesian(ylim = c(0,15))

p = p1+theme(axis.title.y.right = element_blank(),axis.text.y.right = element_blank())+xlab(NULL)+
  p2+theme(axis.title.y.left = element_blank(),axis.text.y.left = element_blank())+xlab(NULL)+
  patchwork::plot_layout(nrow=1, widths = c(20,4))+
  patchwork::plot_annotation(caption=caption)

ggrrr::gg_save_as(p,  here::here("vignettes/latex/s2/fig/simulation_setup_prev_10_v2"), size = std_size$third)

## Modelled output for all PCV groups ----

sim_out_1 = bind_rows(
  purrr::map(c("bayes","lang-reiczigel","rogan-gladen"), function(m, ...) { 
    bind_rows(purrr::map(c("PCV7","PCV13","PCV15","PCV20"), function(p, ...) {
      obs = scenario$samples %>% inner_join(panels %>% filter(panel_name == p), by=c("test"="comp_test")) %>% rename(result=observed)
      control = scenario$performance %>% inner_join(panels %>% filter(panel_name == p), by=c("test"="comp_test"))
      true_panel_prevalence(
        test_results = obs,
        false_pos_controls = control$false_pos_controls,
        n_controls = control$n_controls,
        false_neg_diseased = control$false_neg_diseased,
        n_diseased = control$n_diseased,
        panel_name = p,
        method = m
      ) %>% 
        mutate(panel_group = p) %>%  
        inner_join(scenario$summary, by="test")
    },.progress = TRUE))
  }, .progress = TRUE)
)

## Supp 2 final figure ----


p1 = demo_bar_plot(sim_out_1 %>% filter(n_components==1 & panel_group=="PCV20"))
p2 = demo_bar_plot(sim_out_1 %>% filter(n_components>1))

p = p1+theme(axis.title.y.right = element_blank(),axis.text.y.right = element_blank())+xlab(NULL)+
  p2+theme(axis.title.y.left = element_blank(),axis.text.y.left = element_blank())+xlab(NULL)+
    patchwork::plot_layout(nrow=1, widths = c(20,4),guides = "collect")+
  patchwork::plot_annotation(caption=caption) & theme(legend.position = 'bottom',legend.justification = "center") &
  coord_cartesian(ylim=c(0,15))

ggrrr::gg_save_as(p,  here::here("vignettes/latex/s2/fig/simulation_result_prev_10_v2"), size = std_size$half)

## Main paper figure 2 ----

p1 = demo_bar_plot(sim_out_1 %>% filter(n_components==1 & panel_group=="PCV20" & stringr::str_starts(prevalence.method,"bayes")))
p2 = demo_bar_plot(sim_out_1 %>% filter(n_components>1 & stringr::str_starts(prevalence.method,"bayes")))

p = p1+theme(axis.title.y.right = element_blank(),axis.text.y.right = element_blank())+xlab(NULL)+
  p2+theme(axis.title.y.left = element_blank(),axis.text.y.left = element_blank())+xlab(NULL)+
  patchwork::plot_layout(nrow=1, widths = c(20,4),guides = "collect")+
  patchwork::plot_annotation(caption=caption) & theme(legend.position = 'bottom',legend.justification = "center") &
  coord_cartesian(ylim=c(0,15))
p
ggrrr::gg_save_as(p,  here::here("vignettes/latex/main/fig/simulation_result_bayes_v2"), size = std_size$half)

## Scenario with multiple panels with different prevalences ----

scenario2 = multi_panel_example(
  n_comp = 20,
  n_controls = 800,
  n_diseased = 26,
  n_samples = 4000,
  comp_spec = 0.9975,
  comp_sens = 0.8,
  comp_dist = ipd_distribution("PCV20") %>% pull(x),
  comp_test = ipd_distribution("PCV20") %>% pull(pneumo.phe_serotype) %>% forcats::as_factor(),
  exact_controls = TRUE
)

# scenario2 %>% unnest(summary) %>% filter(test == "Panel") %>% select(scenario_prev,true_prev) %>% view()

analyse_scenario = function(scenario, methods=c("bayes","lang-reiczigel","rogan-gladen"), comp_sens, comp_spec, cap = 5,comp_lim = c(0,0.10),
                            panel_lim = c(0,0.25)) {
  
  sim_out = bind_rows(
      purrr::map(methods, function(m, ...) { 
        scenario %>% mutate(
          method = m,
          modelled = purrr::map2(samples, performance, 
              ~ true_panel_prevalence(
                   test_results = .x %>% rename(result=observed),
                   # false_pos_controls = .y$false_pos_controls,
                   # n_controls = .y$n_controls,
                   # false_neg_diseased = .y$false_neg_diseased,
                   # n_diseased = .y$n_diseased,
                   sens = comp_sens,
                   spec = comp_spec,
                   panel_name = "Panel",
                   method = m
              ), .progress = TRUE)
          )
      }, .progress = TRUE)
  )
  
  caption = c(
    scenario %>% unnest(design) %>% 
      filter(n_components == 1) %>% 
      transmute(label = sprintf("simulated component sensitivity %1.2f%%, specificity %1.2f%%",design_sens*100, design_spec*100)) %>% 
      pull(label) %>% unique(),
    
    sprintf("prior component sensitivity %s, specificity %s", format(comp_sens), format(comp_spec)),
    
    sim_out %>% unnest(modelled) %>% 
      filter(test == "Panel") %>% 
      group_by(prevalence.method) %>%
      summarise(
        sens.lower = mean(sens.lower),
        sens.median = mean(sens.median),
        sens.upper = mean(sens.upper),
        spec.lower = mean(spec.lower),
        spec.median = mean(spec.median),
        spec.upper = mean(spec.upper)
      ) %>%
      transmute(label = sprintf("%s: estimated panel sensitivity %1.2f%% [%1.2f%% \u2013 %1.2f%%], specificity %1.2f%% [%1.2f%% \u2013 %1.2f%%]",
                                prevalence.method,sens.median*100,sens.lower*100,sens.upper*100,spec.median*100,spec.lower*100,spec.upper*100)) %>% 
      pull(label)
  
  ) %>% `[`(1:cap) %>% paste0(collapse="\n")
  
  sim_plot = sim_out %>% transmute(
    scenario_name = scenario_name,
    data = purrr::map2(summary, modelled, ~ .x %>% select(-boot) %>% inner_join(.y, by="test"))
  ) %>% unnest(data)
  
  p1 = testerror:::demo_qq_plot(sim_plot %>% filter(n_components == 1))
  p2 = testerror:::demo_qq_plot(sim_plot %>% filter(n_components != 1))+facet_wrap(~"panel")
  p = p1+coord_fixed(xlim=comp_lim, ylim=comp_lim)+
    p2+coord_fixed(xlim=panel_lim, ylim=panel_lim)+
    patchwork::plot_layout(guides="collect")+
    patchwork::plot_annotation(caption = caption) & 
    theme(legend.position = 'bottom',legend.justification = "center")
  
  return(list(p, p1, p2, caption))

}


tmp = analyse_scenario(scenario2,
                       comp_sens = beta_dist(p=0.8, n=26), comp_spec = beta_dist(p=0.9975, n=800))
ggrrr::gg_save_as(tmp[[1]], here::here("vignettes/latex/s2/fig/simulation_result_sens_80_80_v2"),size = std_size$half)

# Sens = 60 ----

scenario3 = multi_panel_example(
  n_comp = 20,
  n_controls = 800,
  n_diseased = 200,
  n_samples = 4000,
  comp_spec = 0.9975,
  comp_sens = 0.6,
  comp_dist = ipd_distribution("PCV20") %>% pull(x),
  comp_test = ipd_distribution("PCV20") %>% pull(pneumo.phe_serotype) %>% forcats::as_factor(),
  exact_controls = TRUE
)

p_60 = analyse_scenario(scenario3, comp_sens = beta_dist(p=0.6, n=26), comp_spec = beta_dist(p=0.9975, n=800), cap=2)

# Sens = 75 ----

scenario4 = multi_panel_example(
  n_comp = 20,
  n_controls = 800,
  n_diseased = 200,
  comp_spec = 0.9975,
  comp_sens = 0.75,
  comp_dist = ipd_distribution("PCV20") %>% pull(x),
  comp_test = ipd_distribution("PCV20") %>% pull(pneumo.phe_serotype) %>% forcats::as_factor(),
  exact_controls = TRUE
)

p_75 = analyse_scenario(scenario4, comp_sens = beta_dist(p=0.75, n=26), comp_spec = beta_dist(p=0.9975, n=800), cap=2)

# Sens = 90 ----

scenario5 = multi_panel_example(
  n_comp = 20,
  n_controls = 800,
  n_diseased = 200,
  comp_spec = 0.9975,
  comp_sens = 0.9,
  comp_dist = ipd_distribution("PCV20") %>% pull(x),
  comp_test = ipd_distribution("PCV20") %>% pull(pneumo.phe_serotype) %>% forcats::as_factor(),
  exact_controls = TRUE
)

p_90 = analyse_scenario(scenario5, comp_sens = beta_dist(p=0.9, n=26), comp_spec = beta_dist(p=0.9975, n=800), cap=2)


# ggrrr::gg_save_as(tmp[[1]], here::here("vignettes/latex/s2/fig/simulation_result_sens_60_60_v2"),size = std_size$half)

comp_lim = c(0,0.10)
panel_lim = c(0,0.25)

p = p_60[[2]] + coord_fixed(xlim=comp_lim, ylim=comp_lim) + no_x() + #facet_grid("sens: 0.60"~"components")+
  p_60[[3]] + coord_fixed(xlim=panel_lim, ylim=panel_lim) + no_x() + facet_grid("component sens: 60%"~"panel")+
  p_75[[2]] + coord_fixed(xlim=comp_lim, ylim=comp_lim) + no_x() + #facet_grid("sens: 0.75"~"components")+
  p_75[[3]] + coord_fixed(xlim=panel_lim, ylim=panel_lim) + no_x() + facet_grid("component sens: 75%"~"panel")+
  p_90[[2]] + coord_fixed(xlim=comp_lim, ylim=comp_lim) + #facet_grid("sens: 0.90"~"components")+
  p_90[[3]] + coord_fixed(xlim=panel_lim, ylim=panel_lim) + facet_grid("component sens: 90%"~"panel")+
  plot_layout(ncol=2, guides="collect")+
  plot_annotation(tag_levels = "A", caption=sprintf(
"in all scenarios component test specificity is constant at 99.75%%
and the prior for specificity is %s
in the 60%% simulation the prior for sensitivity is %s
in the 75%% simulation the prior for sensitivity is %s
in the 90%% simulation the prior for sensitivity is %s
", format(beta_dist(p=0.9975, n=800)), format(beta_dist(p=0.6, n=26)), format(beta_dist(p=0.75, n=26)), format(beta_dist(p=0.9, n=26))
  ))&
  theme(legend.position = "bottom",legend.justification = "center")

ggrrr::gg_save_as(p, here::here("vignettes/latex/s2/fig/simulation_result_same_sens_v2"),size = std_size$full)


# Sensitivity mismatch ----

# Sens = 75

scenario4 = multi_panel_example(
  n_comp = 20,
  n_controls = 800,
  n_diseased = 200,
  comp_spec = 0.9975,
  comp_sens = 0.75,
  comp_dist = ipd_distribution("PCV20") %>% pull(x),
  comp_test = ipd_distribution("PCV20") %>% pull(pneumo.phe_serotype) %>% forcats::as_factor(),
  exact_controls = TRUE
)

p_60_75 = analyse_scenario(scenario4, comp_sens = beta_dist(p=0.6, n=26), comp_spec = beta_dist(p=0.9975, n=800), cap=2)
p_90_75 = analyse_scenario(scenario4, comp_sens = beta_dist(p=0.9, n=26), comp_spec = beta_dist(p=0.9975, n=800), cap=2)

p = 
  p_60_75[[2]]+ coord_fixed(xlim=comp_lim, ylim=comp_lim)+ no_x()+ 
  p_60_75[[3]]+ coord_fixed(xlim=panel_lim, ylim=panel_lim)+ no_x()+ facet_grid("component sens prior: 60%"~"panel")+
  p_90_75[[2]]+ coord_fixed(xlim=comp_lim, ylim=comp_lim)+ 
  p_90_75[[3]]+ coord_fixed(xlim=panel_lim, ylim=panel_lim)+ facet_grid("component sens prior: 90%"~"panel")+
  plot_layout(ncol=2,guides="collect")+
  plot_annotation(tag_levels = "A", caption=
    "in all scenarios component test specificity is constant at 99.75%, and sensitivity at 75%"
  )&theme(legend.position = "bottom",legend.justification = "center")

ggrrr::gg_save_as(p, here::here("vignettes/latex/s2/fig/bayesian_sim_mismatch_sens_v2"),size = std_size$two_third)

# Specificity mismatch ----

p_999_9975 = analyse_scenario(scenario4, comp_sens = beta_dist(p=0.75, n=26), comp_spec = beta_dist(p=0.999, n=800), cap=2)
p_995_9975 = analyse_scenario(scenario4, comp_sens = beta_dist(p=0.75, n=26), comp_spec = beta_dist(p=0.995, n=800), cap=2)

p = 
  p_995_9975[[2]]+ coord_fixed(xlim=comp_lim, ylim=comp_lim)+ no_x()+ 
  p_995_9975[[3]]+ coord_fixed(xlim=panel_lim, ylim=panel_lim)+ no_x()+ facet_grid("component spec prior: 99.5%"~"panel")+
  p_999_9975[[2]]+ coord_fixed(xlim=comp_lim, ylim=comp_lim)+ 
  p_999_9975[[3]]+ coord_fixed(xlim=panel_lim, ylim=panel_lim)+ facet_grid("component sens prior: 99.9%"~"panel")+
  plot_layout(ncol=2,guides="collect")+
  plot_annotation(tag_levels = "A", caption=
                    "in all scenarios component test specificity is constant at 99.75%, and sensitivity at 75%"
  )&theme(legend.position = "bottom",legend.justification = "center")

ggrrr::gg_save_as(p, here::here("vignettes/latex/s2/fig/bayesian_sim_mismatch_spec_v2"),size = std_size$two_third)

p = p_60_75[[3]]+ facet_grid("sim: sens 75%; spec 99.75%" ~ "priors: sens 60%; spec 99.75%")+no_x()+
  p_90_75[[3]]+ facet_grid("sim: sens 75%; spec 99.75%" ~"priors: sens 90%; spec 99.75%")+no_x()+no_y()+
  p_995_9975[[3]]+facet_grid("sim: sens 75%; spec 99.75%" ~"priors: sens 75%; spec 99.5%")+
  p_999_9975[[3]]+facet_grid("sim: sens 75%; spec 99.75%" ~"priors: sens 75%; spec 99.9%")+no_y()+
  plot_layout(ncol=2,guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom",legend.justification = "center")&coord_fixed(xlim=panel_lim, ylim=panel_lim)

ggrrr::gg_save_as(p, here::here("vignettes/latex/s2/fig/bayesian_sim_mismatch_v2"),size = std_size$two_third)


# scenario = do_scenario(n_controls=800, n_diseased=26)
# 
# # Just use scenario setup to describe sens uncertainty --
# # TODO: need to figure out what te fit looks like.
# 
# scenario_config = scenario %>% select(-id,-boot,-actual,-test) %>% unnest(test_pcv_group) %>% distinct()
# 
# scenario_config_2 = scenario_config %>% 
#   glimpse() %>%
#   group_by(across(c(starts_with("group"),starts_with("panel")))) %>%
#   summarise(
#     panel_size = n(),
#     panel_prevalence = testerror::panel_prevalence(test_prevalence),
#     panel_sens = testerror::panel_sens(test_prevalence, test_sens, test_spec),
#     panel_expected_apparent_prevalence = testerror::panel_prevalence(test_expected_apparent_prevalence),
#     panel_spec = testerror::panel_spec(test_spec),
#     panel_sens_estimator = testerror::panel_sens_estimator(ap = test_expected_apparent_prevalence,sens = test_sens, spec = test_spec),
#     panel_sens_samples = list(testerror::uncertain_panel_sens_estimator(
#       pos_obs = test_expected_apparent_prevalence*1000, n_obs = 1000,
#       false_pos_controls = test_false_pos_controls, n_controls = test_n_controls,
#       false_neg_diseased = test_false_neg_diseased, n_diseased = test_n_diseased
#     )),
#     panel_spec_samples = list(testerror::uncertain_panel_spec(
#       false_pos_controls = test_false_pos_controls, n_controls = test_n_controls
#     ))
#   ) %>%
#   group_by(across(c(starts_with("group"),starts_with("panel")))) %>%
#   mutate(
#     panel_sens_beta = purrr::map(panel_sens_samples, ~ as_tibble(beta_fit(.x))),
#     panel_spec_beta = purrr::map(panel_spec_samples, ~ as_tibble(beta_fit(.x))),
#     testerror::uncertain_rogan_gladen(
#       pos_obs = panel_expected_apparent_prevalence*1000, n_obs = 1000,
#       sens = panel_sens_samples,
#       spec = panel_spec_samples,
#       prefix = "rogan_gladen"
#     )
#   ) %>%
#   unnest(cols = c(panel_sens_beta,panel_spec_beta),names_sep = ".") %>%
#   mutate(
#     testerror::true_prevalence(
#       pos_obs = panel_expected_apparent_prevalence*1000, n_obs = 1000,
#       false_pos_controls = panel_spec_beta.shape2, n_controls = panel_spec_beta.conc,
#       false_neg_diseased = panel_sens_beta.shape2, n_diseased = panel_sens_beta.conc,
#       prefix = "lang_reiczigel"
#     )
#   ) %>%
#   glimpse()

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

# Compile the multiple tests model --

# stan_model_combined = rstan::stan_model(file = here::here("vignettes/multiple-tests.stan"))
# 
# sampling = memoise::memoise(rstan::sampling,cache = cd)
# 
# # The simulated data (which is grouped) needs to be fed to stan.
# # the data has a patient id (id) and a test id (test_id)
# do_bayesian_model = function(d, g = tibble::tibble(), idCol = id, testIdCol = test_id, resultCol = test, sens = 0.8, spec = 0.9975, 
#                              #false_pos_controls = (n_controls-2/3)*(1-spec)+1/3,
#                              false_pos_controls = (n_controls)*(1-spec),
#                              n_controls = 1/(1-spec),
#                              #false_neg_diseased = (n_diseased-2/3)*(1-sens)+1/3,
#                              false_neg_diseased = (n_diseased)*(1-sens),
#                              n_diseased = 1/(1-sens), 
#                              controls = NULL,
#                              model = stan_model_combined,
#                              ...) {
#   idCol = rlang::enexpr(idCol)
#   testIdCol = rlang::enexpr(testIdCol)
#   resultCol = rlang::enexpr(resultCol)
#   
#   tmp = d %>%
#     transmute(
#       id = !!idCol,
#       test_id = !!testIdCol,
#       result = !!resultCol
#     ) 
#   if (is.factor(tmp$test_id)) {
#     tests = levels(forcats::fct_drop(tmp$test_id))
#   } else {
#     tests = unique(tmp$test_id)
#   }
#   n_test = length(tests)
#   
#   tmp = tmp %>%
#     pivot_wider(names_from = test_id, values_from = result,values_fill = NA_integer_) %>%
#     select(-id) %>% 
#     as.matrix()
#   
#   if (!is.null(controls)) {
#     
#     y_control = controls %>%
#       transmute(
#         id = !!idCol,
#         test_id = !!testIdCol,
#         result = !!resultCol
#       ) %>%
#       pivot_wider(names_from = test_id, values_from = result,values_fill = NA_integer_) %>%
#       select(-id) %>% as.matrix()
#     
#   } else {
#     y_control =  matrix(numeric(), ncol = n_test, nrow = 0, byrow = TRUE)
#   }
#   
#   est_data_combined = list(
#     n_test = n_test,
#     k_sample = nrow(tmp),
#     y_sample = tmp,
#     k_control = nrow(y_control),
#     y_control = y_control,
#     # n_disease = nrow(y_disease),
#     # y_disease = y_disease,
#     # Priors on sensitivity / specificity
#     a_spec = n_controls-false_pos_controls,
#     b_spec = false_pos_controls,
#     a_sens = n_diseased-false_neg_diseased,
#     b_sens = false_neg_diseased
#   )
#   
#   sens_prior = do.call(sprintf, c(list(fmt = "%1.4f [%1.4f \u2013 %1.4f]"), 1-false_neg_diseased/n_diseased, 
#                                   qbeta(c(0.025,0.975), shape1 = n_diseased-false_neg_diseased, shape2 = false_neg_diseased)))
#   
#   spec_prior = do.call(sprintf, c(list(fmt = "%1.4f [%1.4f \u2013 %1.4f]"), 1-false_pos_controls/n_controls, 
#                                   qbeta(c(0.025,0.975), shape1 = n_controls-false_pos_controls, shape2 = false_pos_controls)))
#   
#   
#   fit_combined = sampling(
#       model,
#       data = est_data_combined,
#       chains = 4,
#       warmup = 1000,          # number of warmup iterations per chain
#       iter = 2000,            # total number of iterations per chain
#       show_messages = TRUE
#     )
#   
#   summ = summary(fit_combined, pars = c("p","sens","spec"))$summary
#   
#   summ2 = as_tibble(summ,rownames = "param") %>% 
#     mutate(
#       test_id = tests[as.integer(stringr::str_extract(param,"[0-9]+"))],
#       param = stringr::str_extract(param,"[a-zA-Z]+")
#     )
#   
#   comb_summ = summary(fit_combined, pars = c("p_combined","sens_combined","spec_combined"))$summary
#   comb_summ2 = as_tibble(comb_summ,rownames = "param") %>% 
#     mutate(
#       param = stringr::str_extract(param,"[a-zA-Z]+")
#     )
#   
#   return(
#     tibble(
#       data = list(d),
#       standata = list(est_data_combined),
#       stanfit = list(fit_combined),
#       component = list(summ2),
#       panel = list(comb_summ2),
#       tests = list(tests),
#       sens_prior = sens_prior,
#       spec_prior = spec_prior
#     )
#   )
# }


# Test the multiple tests model --
# 
# bayesian_result = scenario %>% 
#   unnest(test_pcv_group) %>%
#   #glimpse() %>%
#   #filter(panel_id == "PCV20") %>%
#   group_by(group, panel_id) %>%
#   group_modify(do_bayesian_model, sens=0.8, spec = 0.9975, n_controls=800, n_diseased=26)
# 
# 
# do_combine_results = function(bayesian_result) {
#   
#   grp = bayesian_result %>% groups()
#   serotype_tests = bayesian_result %>% 
#     mutate(
#       prior_false_pos_controls = purrr::map_dbl(standata, ~ .x$b_spec),
#       prior_n_controls = purrr::map_dbl(standata, ~ .x$a_spec+.x$b_spec),
#       prior_false_neg_diseased = purrr::map_dbl(standata, ~ .x$b_sens),
#       prior_n_diseased = purrr::map_dbl(standata, ~ .x$a_sens+.x$b_sens),
#       prior_spec = purrr::map_dbl(standata, ~ .x$a_spec/(.x$a_spec+.x$b_spec)),
#       prior_sens = purrr::map_dbl(standata, ~ .x$a_sens/(.x$a_sens+.x$b_sens)),
#     ) %>% 
#     select(!!!grp, data, starts_with("prior")) %>% 
#     unnest(data)
#   
#   disagg_prev = serotype_tests %>%
#     group_by(!!!grp, across(c(starts_with("prior"), starts_with("group"),starts_with("test"),-test))) %>%
#     summarise(
#       actual = sum(actual),
#       test = sum(test),
#       total = n()
#     ) %>%
#     group_by(!!!grp, across(c(starts_with("prior"), starts_with("group"),starts_with("test"),-test))) %>%
#     mutate(
#       actual_prevalence = actual/total,
#       prev.est = testerror::rogan_gladen(test/total, test_sens, test_spec),
#       testerror::true_prevalence(pos_obs = test, n_obs = total, false_pos_controls = prior_false_pos_controls, n_controls = prior_n_controls, false_neg_diseased = prior_false_neg_diseased, n_diseased = prior_n_diseased,),
#       testerror::true_prevalence(pos_obs = test, n_obs = total, false_pos_controls = prior_false_pos_controls, n_controls = prior_n_controls, false_neg_diseased = prior_false_neg_diseased, n_diseased = prior_n_diseased,method = "rogan-gladen", prefix="rogan_gladen")
#     )
#   
#   combined_prev = serotype_tests %>%
#     mutate(prevalence = test_prevalence) %>%
#     group_by(!!!grp, across(starts_with("group"))) %>%
#     estimate_panel_performance_uncertain(test_id) %>%
#     mutate(
#       prev_est.0.5 = testerror::rogan_gladen(panel_apparent_prevalence, panel_sens_est, panel_spec),
#       panel_actual_prevalence = actual/total
#     )
#   
#   comp_tmp = bayesian_result %>% ungroup() %>% select(!!!grp,group,component) %>% unnest(component) %>% filter(param=="p") %>%
#     left_join(disagg_prev %>% mutate(test_id = as.character(test_id)), by=join_by(!!!grp, test_id)) %>%
#     mutate(binom_ci_2(test,total,"apparent"))
#   
#   panel_tmp = bayesian_result %>% ungroup() %>% select(!!!grp,group,panel) %>% unnest(panel) %>%  filter(param=="p") %>%
#     inner_join(combined_prev,by=join_by(!!!grp)) %>%
#     mutate(binom_ci_2(test,total,"apparent"))
#   
#   return(list(
#     components = comp_tmp,
#     panels = panel_tmp,
#     priors = sprintf("prior:\nsens: %s\nspec: %s", unique(bayesian_result$sens_prior), unique(bayesian_result$spec_prior)),
#     params = sprintf("simulation:\nsens: %1.4f\nspec: %1.4f", unique(serotype_tests$test_sens),unique(serotype_tests$test_spec)),
#     raw_result = bayesian_result
#   ))
# }
# 
# combined = do_combine_results(bayesian_result)

# do_plots = function(combined, show=c("bayes","lang-reiczigel","rogan-gladen")) {
# 
#   comp_tmp = combined$components %>% filter(panel_id == "PCV20")
#   panel_tmp = combined$panels %>% filter(panel_id == "PCV20")
#   
#   sep = 0.002
#   
#   # b = "bayes" %in% show
#   # l = "lang-reiczigel" %in% show
#   
#   offset = seq(-(length(show)-1)*sep/2,(length(show)-1)*sep/2,length.out=length(show))
#   
#   b = c(offset[which(show=="bayes")],0)[[1]]
#   l = c(offset[which(show=="lang-reiczigel")],0)[[1]]
#   r = c(offset[which(show=="rogan-gladen")],0)[[1]]
#   
#   points = bind_rows(
#     comp_tmp %>% transmute(x=actual_prevalence+b/2.5, test_id = test_id, type="bayes", lower = `2.5%`, mid = `50%`, upper=`97.5%`),
#     comp_tmp %>% transmute(x=actual_prevalence+r/2.5, test_id = test_id, type="rogan-gladen", lower = rogan_gladen.lower, mid = rogan_gladen.median, upper=rogan_gladen.upper),
#     comp_tmp %>% transmute(x=actual_prevalence+l/2.5, test_id = test_id, type="lang-reiczigel", lower = prevalence.lower, mid = prevalence.median, upper=prevalence.upper)
#   )
#   
#   points2 = bind_rows(
#     panel_tmp %>% transmute(x=panel_actual_prevalence+b, panel_id = panel_id, type="bayes", lower = `2.5%`, mid = `50%`, upper=`97.5%`),
#     panel_tmp %>% transmute(x=panel_actual_prevalence+r, panel_id = panel_id, type="rogan-gladen", lower = rogan_gladen.lower, mid = rogan_gladen.median, upper=rogan_gladen.upper),
#     panel_tmp %>% transmute(x=panel_actual_prevalence+l, panel_id = panel_id, type="lang-reiczigel", lower = prevalence.lower, mid = prevalence.median, upper=prevalence.upper)
#   )
#   
#   points = points %>% filter(type %in% show)
#   points2 = points2 %>% filter(type %in% show)
#   
#   cilim = binom::binom.confint(seq(0,0.25,length.out=51)*1000,rep(1000,51),methods="wilson")
#   
#   p1 = ggplot(points, aes(x=x, y= mid, ymin=lower, ymax=upper,colour=type))+
#     geom_point(size=0.5)+
#     geom_errorbar(width = sep/2.5, size=0.25)+
#     coord_fixed(xlim=c(0,0.1),ylim=c(0,0.1))+
#     #geom_point(aes(y=prev.est), colour="magenta", size=1)+
#     #geom_errorbar(aes(ymin=apparent.0.025, ymax=apparent.0.975), colour="red",position=position_nudge(x=-0.001/2.5), width = 0, alpha=0.4)+
#     # geom_point(aes(y=prevalence.median),position=position_nudge(x=-sep/2.5),  size=0.5)+
#     # geom_errorbar(aes(ymin=prevalence.lower, ymax=prevalence.upper,linetype="lang-reiczigel"),position=position_nudge(x=-sep/2.5), width = 0, alpha=0.4)+
#     annotate("text", x=0.005, y=0.095, label = combined$params, vjust="inward", hjust="inward")+
#     annotate("text", x=0.095, y=0.005, label = combined$priors, vjust="inward", hjust="inward")+
#     geom_point(aes(x=actual_prevalence, y=apparent.0.5), data=comp_tmp, colour="red",shape=4, size=2, inherit.aes = FALSE)+
#     geom_abline(colour="#8080FF")+
#     # geom_line(aes(x=x/n,y=upper),data = cilim, colour="#8080FF",linetype="dotted", inherit.aes = FALSE)+
#     # geom_line(aes(x=x/n,y=lower),data = cilim, colour="#8080FF",linetype="dotted", inherit.aes = FALSE)+
#     xlab("simulation prevalence")+
#     ylab("estimated prevalence")+
#     facet_wrap(~"components")+
#     guides(colour= guide_legend(title=element_blank()))+
#     scale_color_grey(start = 0, end = 0.6)+
#     theme(legend.position = "bottom")
#   
#   p2 = ggplot(points2, aes(x=x, y= mid, ymin=lower, ymax=upper,colour=type))+
#     geom_point(size=0.5)+
#     geom_errorbar(width = sep, size=0.25)+
#     coord_fixed(xlim=c(0,0.25),ylim=c(0,0.25))+
#     #geom_point(aes(y=prev_est.0.5), colour="magenta", size=1)+
#     #geom_errorbar(aes(ymin=apparent.0.025, ymax=apparent.0.975), colour="red",position=position_nudge(x=-0.001), width = 0, alpha=0.4)+
#     # geom_point(aes(y=prevalence.median),position=position_nudge(x=-sep),  size=0.5)+
#     # geom_errorbar(aes(ymin=prevalence.lower, ymax=prevalence.upper, linetype="lang-reiczigel"),position=position_nudge(x=-sep), width = 0, alpha=0.4)+
#     annotate("text", x=0.01, y=0.24, label = combined$params, vjust="inward", hjust="inward")+
#     annotate("text", x=0.24, y=0.01, label = combined$priors, vjust="inward", hjust="inward")+
#     geom_point(aes(x=panel_actual_prevalence, y=apparent.0.5), data=panel_tmp, colour="red",shape=4,  size=2, inherit.aes = FALSE)+
#     geom_abline(colour="#8080FF")+
#     # geom_line(aes(x=x/n,y=upper),data = cilim, colour="#8080FF",linetype="dotted", inherit.aes = FALSE)+
#     # geom_line(aes(x=x/n,y=lower),data = cilim, colour="#8080FF",linetype="dotted", inherit.aes = FALSE)+
#     xlab("simulation prevalence")+
#     ylab("estimated prevalence")+
#     facet_wrap(~"panel")+
#     guides(colour= guide_legend(title=element_blank()))+
#     scale_color_grey(start = 0, end = 0.6)+
#     theme(legend.position = "bottom")
#   
#   p3 = p1+p2+patchwork::plot_layout(guides="collect") & theme(legend.position = 'bottom',legend.justification = "center")
#   # browser()
#   return(list(p1,p2,p3))
# }
# 
# p = do_plots(combined) #, show="bayes")
# p
# 
# save_as(p[[3]], here::here("vignettes/latex/s2/fig/simulation_result_sens_80_80"),size = std_size$third)
# 
# do_bar_plot = function(combined, prev_filter = 0.1, show=c("bayes","lang-reiczigel","rogan-gladen")) {
#   
#   tests = # combined$raw_result %>% filter(panel_id == "PCV20") %>% pull(tests) %>% `[[`(1)
#     avoncap::serotype_data$xr %>% filter(order==1) %>% pull(serotype)
#   
#   comp_tmp = combined$components %>% 
#     filter(group_prevalence == prev_filter) %>%
#     mutate(test_id = factor(test_id, tests)) %>% filter(panel_id == "PCV20" & !is.na(test_id))
#   panel_tmp = combined$panels %>% filter(group_prevalence == prev_filter) %>%
#     mutate(panel_id = factor(panel_id, c("PCV7","PCV13","PCV15","PCV20")))
#   
#   points = bind_rows(
#     comp_tmp %>% transmute(test_id = test_id, type="bayes", lower = `2.5%`*100, mid = `50%`*100, upper=`97.5%`*100),
#     comp_tmp %>% transmute(test_id = test_id, type="lang-reiczigel", lower = prevalence.lower*100, mid = prevalence.median*100, upper=prevalence.upper*100),
#     comp_tmp %>% transmute(test_id = test_id, type="rogan-gladen", lower = rogan_gladen.lower*100, mid = rogan_gladen.median*100, upper=rogan_gladen.upper*100)
#   )
#   
#   points2 = bind_rows(
#     panel_tmp %>% transmute(panel_id = panel_id, type="bayes", lower = `2.5%`*100, mid = `50%`*100, upper=`97.5%`*100),
#     panel_tmp %>% transmute(panel_id = panel_id, type="lang-reiczigel", lower = prevalence.lower*100, mid = prevalence.median*100, upper=prevalence.upper*100),
#     panel_tmp %>% transmute(panel_id = panel_id, type="rogan-gladen", lower = rogan_gladen.lower*100, mid = rogan_gladen.median*100, upper=rogan_gladen.upper*100)
#   )
#   
#   points = points %>% filter(type %in% show)
#   points2 = points2 %>% filter(type %in% show)
#   
#   p3 = ggplot(comp_tmp, aes(x=test_id, y=test/total*100))+
#     geom_bar(stat="identity",fill="grey90", colour="black",size=0.05)+
#     geom_errorbar(aes(ymin=test/total*100,ymax=test/total*100),colour="red")+
#     geom_errorbar(aes(ymin=actual/total*100,ymax=actual/total*100),colour="#8080FF")+
#     geom_point(aes(x=test_id, y=mid, colour=type),data=points,inherit.aes = FALSE,size=0.5,position = position_dodge(width=0.4))+
#     geom_errorbar(aes(x=test_id, ymin=lower, ymax=upper, colour=type),data=points,inherit.aes = FALSE,width=0,position = position_dodge(width=0.4))+
#     
#     ggpp::annotate("text_npc", npcx = 0.05, npcy = 0.95, label = sprintf("prevalence: %1.2g\n\n%s \n\n%s", prev_filter, combined$params, combined$priors), size=6/ggplot2::.pt )+
#     xlab(NULL)+ylab("prevalence (%)")+coord_cartesian(ylim=c(0,15))+
#     guides(colour= guide_legend(title=element_blank()))+
#     scale_color_grey(start = 0, end = 0.6)
#   
#   p4 = ggplot(panel_tmp, 
#               aes(x=panel_id, y=test/total*100))+
#     geom_bar(stat="identity",fill="grey90", colour="black",size=0.05)+
#     geom_errorbar(aes(ymin=test/total*100,ymax=test/total*100),colour="red")+
#     geom_errorbar(aes(ymin=actual/total*100,ymax=actual/total*100),colour="#8080FF")+
#     geom_point(aes(x=panel_id, y=mid, colour=type),data=points2,inherit.aes = FALSE,size=0.5,position = position_dodge(width=0.4))+
#     geom_errorbar(aes(x=panel_id, ymin=lower, ymax=upper, colour=type),data=points2,inherit.aes = FALSE,width=0,position = position_dodge(width=0.4))+
#     
#     xlab(NULL)+ylab("prevalence (%)")+coord_cartesian(ylim=c(0,15))+
#     theme(axis.text.x.bottom = element_text(angle=45,vjust = 1,hjust=1))+
#     guides(colour= guide_legend(title=element_blank()))+
#     scale_color_grey(start = 0, end = 0.6)
#   
#   p5 = p3+p4+no_y()+patchwork::plot_layout(nrow=1,widths = c(5,1),guides="collect")&theme(legend.position = "bottom")
#   
#   return(list(p3,p4,p5))
# }
# 
# p2 = do_bar_plot(combined)
# 
# save_as(p2[[3]], here::here("vignettes/latex/s2/fig/simulation_result_prev_10"),size = std_size$third)
# 
# # Other scenarios --
# 
# # intro in supplementary
# p5 = do_bar_plot(combined, show=NULL)
# save_as(p5[[3]], here::here("vignettes/latex/s2/fig/simulation_setup_prev_10"),size = std_size$third)
# 
# # main paper
# p5 = do_bar_plot(combined, show="bayes")
# save_as(p5[[3]]&guides(colour=guide_none()), here::here("vignettes/latex/main/fig/simulation_result_bayes"),size = std_size$third)

## Sens = 60 --

# serotype_tests_60 = do_scenario(
#   spec = 0.9975,n_controls = 800,
#   sens = 0.6,n_diseased = 260
# )
# 
# model_60 = serotype_tests_60 %>% 
#   unnest(test_pcv_group) %>%
#   filter(panel_id == "PCV20") %>%
#   group_by(group, panel_id) %>%
#   group_modify(do_bayesian_model, sens=0.6, spec = 0.9975, n_controls=800, n_diseased=26)
# 
# p_60 = model_60 %>% 
#   do_combine_results() %>%
#   do_plots()
# 
# # p_60
# # save_as(p_60[[3]], here::here("vignettes/latex/s2/fig/simulation_result_sens_80_80"),size = std_size$third)
# 
# ## Sens = 75 --
# 
# serotype_tests_75 = do_scenario(
#   spec = 0.9975,n_controls = 800,
#   sens = 0.75,n_diseased = 260
# )
# 
# model_75 = serotype_tests_75 %>% 
#   unnest(test_pcv_group) %>%
#   filter(panel_id == "PCV20") %>%
#   group_by(group, panel_id) %>%
#   group_modify(do_bayesian_model, sens=0.75, spec = 0.9975, n_controls=800, n_diseased=26)
# 
# p_75 = model_75 %>% 
#   do_combine_results() %>%
#   do_plots()
# 
# # p_75
# 
# ## Sens = 90 --
# 
# serotype_tests_90 = do_scenario(
#   spec = 0.9975,n_controls = 800,
#   sens = 0.90,n_diseased = 260
# )
# 
# model_90 = serotype_tests_90 %>% 
#   unnest(test_pcv_group) %>%
#   filter(panel_id == "PCV20") %>%
#   group_by(group, panel_id) %>%
#   group_modify(do_bayesian_model, sens=0.90, spec = 0.9975, n_controls=800, n_diseased=26)
# 
# p_90 = model_90 %>% 
#   do_combine_results() %>%
#   do_plots()
# 
# # p_90
# 
# p_90_3 = model_90 %>% 
#   do_combine_results() %>%
#   do_bar_plot()
# 
# # p_90_3
# 
# ## Combination of alternative scenario plots --
# 
# p = p_60[[2]]+ no_x()+ facet_grid("sens: 0.60"~"components")+
#   p_60[[3]]+ no_x()+ facet_grid("sens: 0.60"~"panel")+
#   p_75[[1]]+ no_x()+ #facet_grid("sens: 0.75"~"components")+
#   p_75[[2]]+ no_x()+ #facet_grid("sens: 0.75"~"panel")+
#   p_90[[1]]+ #facet_grid("sens: 0.90"~"components")+
#   p_90[[2]]+ #facet_grid("sens: 0.90"~"panel")+
#   plot_layout(ncol=2, guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")
# 
# save_as(p, here::here("vignettes/latex/s2/fig/simulation-result-same-sens"),size = std_size$full)


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

# Sensitivity mismatch --
# 
# serotype_tests_60_90 = do_scenario(
#   spec = 0.9975,n_controls = 800,
#   sens = 0.60,n_diseased = 260
# )
# 
# model_60_90 = serotype_tests_60_90 %>% 
#   unnest(test_pcv_group) %>%
#   filter(panel_id == "PCV20") %>%
#   group_by(group, panel_id) %>%
#   group_modify(do_bayesian_model, sens=0.90, spec = 0.9975, n_controls=800, n_diseased=26)
# 
# p_60_90 = model_60_90 %>%
#   do_combine_results() %>%
#   do_plots()
# 
# # p_60_90
# 
# serotype_tests_90_60 = do_scenario(
#   spec = 0.9975,n_controls = 800,
#   sens = 0.90,n_diseased = 260
# )
# 
# model_90_60 = serotype_tests_90_60 %>% 
#   unnest(test_pcv_group) %>%
#   filter(panel_id == "PCV20") %>%
#   group_by(group, panel_id) %>%
#   group_modify(do_bayesian_model, sens=0.60, spec = 0.9975, n_controls=800, n_diseased=26)
# 
# p_90_60 = model_90_60 %>%
#   do_combine_results() %>%
#   do_plots() #show="bayes")
# 
# p_90_60

# p = p_60_90[[2]]+ #facet_grid("sens: 0.60"~"panel")+
#   p_90_60[[2]]+no_y()+ #facet_grid("sens: 0.90"~"panel")+
#   plot_layout(ncol=2, guides="collect")+plot_annotation(tag_levels = "A")

# p = p_60_90[[1]]+p_60_90[[2]]+ #facet_grid("sens: 0.60"~"panel")+
#   p_90_60[[1]]+p_90_60[[2]]+ #facet_grid("sens: 0.90"~"panel")+
#   plot_layout(ncol=2,guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")
# 
# 
# save_as(p, here::here("vignettes/latex/s2/fig/bayesian-sim-mismatch-sens"),size = std_size$two_third)

# Specificity mismatch --

# serotype_tests_9975_99 = do_scenario(
#   spec = 0.9975,n_controls = 800,
#   sens = 0.80,n_diseased = 260
# )
# 
# model_9975_99 = serotype_tests_9975_99 %>% 
#   unnest(test_pcv_group) %>%
#   filter(panel_id == "PCV20") %>%
#   group_by(group, panel_id) %>%
#   group_modify(do_bayesian_model, sens=0.80, spec = 0.99, n_controls=100, n_diseased=26)
# 
# p_9975_99 = model_9975_99 %>%
#   do_combine_results() %>%
#   do_plots()
# 
# 
# 
# serotype_tests_99_9975 = do_scenario(
#   spec = 0.99,n_controls = 800,
#   sens = 0.80,n_diseased = 260
# )
# 
# model_99_9975 = serotype_tests_99_9975 %>% 
#   unnest(test_pcv_group) %>%
#   filter(panel_id == "PCV20") %>%
#   group_by(group, panel_id) %>%
#   group_modify(do_bayesian_model, sens=0.80, spec = 0.9975, n_controls=100, n_diseased=26)
# 
# p_99_9975 = model_99_9975 %>%
#   do_combine_results() %>%
#   do_plots()
# 
# p = p_9975_99[[1]]+p_9975_99[[2]]+ #facet_grid("sens: 0.60"~"panel")+
#   p_99_9975[[1]]+p_99_9975[[2]]+ #facet_grid("sens: 0.90"~"panel")+
#   plot_layout(ncol=2,guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")
# 
# save_as(p, here::here("vignettes/latex/s2/fig/bayesian-sim-mismatch-spec"),size = std_size$two_third)
# 
# 
# p = p_60_90[[2]]+
#   p_90_60[[2]]+p_9975_99[[2]]+p_99_9975[[2]]+ #facet_grid("sens: 0.90"~"panel")+
#   plot_layout(ncol=2,guides="collect")+plot_annotation(tag_levels = "A")&theme(legend.position = "bottom")&coord_cartesian(xlim=c(0,0.25),ylim=c(0,0.4))
# 
# save_as(p, here::here("vignettes/latex/s2/fig/bayesian-sim-mismatch"),size = std_size$two_third)
# 
