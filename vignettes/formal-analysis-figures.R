
# Setup ----

library(tidyverse)
library(patchwork)
here::i_am("vignettes/formal-analysis-figures.R")
source(here::here("vignettes/vignette-utils.R"))
ggrrr::gg_pedantic()
devtools::load_all()

# Tabular combination of tests ----

# .md = function(x) {
#   tmp = sapply(x, function(y) markdown::markdownToHTML(text=y,fragment.only = TRUE))
#   return(tmp)
# }
# 
# .tex = function(x) {
#   tmp = sapply(x, function(y) {
#     t = tth::ttm(y, L = TRUE)
#     t = paste0(t,collapse="")
#     sprintf("<p>%s</p>",t)
#   })
#   return(tmp)
# }
# 
# .eq = function(x) {
#   .tex(sapply(x, sprintf, fmt="$%s$"))
# }
# 
# map = tibble::tribble(
#   ~test1, ~test2, ~combination,
#   "TP", "TP", "TP",
#   "TP", "FP", "TP",
#   "TP", "FN", "TP",
#   "TP", "TN", "TP",
#   "FP", "TP", "TP",
#   "FP", "FP", "FP",
#   "FP", "FN", "TP+",
#   "FP", "TN", "FP",
#   "FN", "TP", "TP",
#   "FN", "FP", "TP+",
#   "FN", "FN", "FN",
#   "FN", "TN", "FN",
#   "TN", "TP", "TP",
#   "TN", "FP", "FP",
#   "TN", "FN", "FN",
#   "TN", "TN", "TN",
# ) %>% mutate(
#   p1 = case_when(
#     test1 == "TP" ~ "P(O_i|A_i) P(A_i)",
#     test1 == "FP" ~ "P(O_i|\\neg A_i) P(\\neg A_i)",
#     test1 == "FN" ~ "P(\\neg O_i|A_i) P(A_i)",
#     test1 == "TN" ~ "P(\\neg O_i|\\neg A_i) P(\\neg A_i)",
#   ),
#   p2 = case_when(
#     test2 == "TP" ~ "P(O_j|A_j) P(A_j)",
#     test2 == "FP" ~ "P(O_j|\\neg A_j) P(\\neg A_j)",
#     test2 == "FN" ~ "P(\\neg O_j|A_j) P(A_j)",
#     test2 == "TN" ~ "P(\\neg O_j|\\neg A_j) P(\\neg A_j)",
#   ),
#   test1group = .eq(ifelse(test1 %in% c("TP","FP"), "O", "\\neg O")),
#   test2group = .eq(ifelse(test2 %in% c("TP","FP"), "O", "\\neg O")),
#   p_combined = sprintf("%s %s",p1,p2)
#   #p_combined = sprintf("\\biggl(%s\\biggr) \\times \\biggl(%s\\biggr)",p1,p2)
# ) 
# 
# t = map %>% 
#   select(` ` = test1group, `  ` = test1, test2group, test2, combination) %>%
#   ggrrr::hux_tidy(vars(` `, `  `),vars(test2group, test2)) %>% 
#   huxtable::map_escape_contents(huxtable::everywhere, huxtable::everywhere, huxtable::by_function(function(s) !stringr::str_starts(s,"<p"))) 
# %>% 
#   huxtable::to_html() %>%
#   htmltools::HTML()
# t %>% ggrrr::hux_save_as(here::here("vignettes/latex/probability-table.pdf"))


# Characterisation ----

## Aggregate vs component sensitivity as with varying prevalence -----
 
# * Assume 10 identically prevalent diseases
# * Fix sensitivity to 75% of each component test; and fix prevalence of all subtypes 
# to 2%
# * Fix sensitivity of components to either 60% or 90%.
# * Varying number of tests with low verus high sensitivity; 


plot_data = tidyr::crossing(
    test = 1:10,
    p = seq(0,0.25,length.out=101),
    sens = seq(0.5,1,length.out=101),
    spec = 0.975
  ) %>% 
  # group by the thing we are varying
  group_by(p,sens) %>%
  mutate(
   `combined sens` = panel_sens(p, sens, spec),
  )

p1 = ggplot(plot_data, aes(x=sens, y=`combined sens`, colour=p, z=p)) + 
  geom_point(size=0.25) + 
  geom_abline(colour="grey50") + coord_fixed(ylim = c(0.5,1))+
  scale_color_viridis_c(option = "magma",name="prevalence")+
  geom_line(data = plot_data %>% filter(p %in% seq(0,0.25,length.out=5)), mapping = aes(x=sens, y=`combined sens`, group=factor(p)), colour="cyan", size=0.25)+
  xlab("component sensitivity")+ylab("combined sensitivity") + 
  facet_wrap(~"prevalence")+theme(legend.box = "horizontal")

#p1

## Aggregate vs component sensitivity as with varying lo/high components -----
 
# * Assume 10 identically prevalent diseases
# * Fix specificity to 99.75% of each component test; 
# * Fix sensitivity of components to either 60% or 90%.
# * Varying number of tests with low verus high sensitivity; 
# * allow prevalence of components to vary from 0 to 25%

plot_data = tidyr::crossing(
    test = 1:10,
    n_low_versus_high = 0:10,
    p = 0.05,
    sens_lo = seq(0.5,1,length.out = 1001),
    spec = 0.975
  ) %>% 
  mutate(
    sens_hi = (1-(1-sens_lo)*0.2),
    sens = ifelse(test <= n_low_versus_high, sens_lo, sens_hi),
    ratio = factor(n_low_versus_high+1, labels = sapply(0:10, function(x) sprintf("%d:%d",x,10-x)))
  ) %>%
  # group by the thing we are varying
  group_by(n_low_versus_high, sens_lo, sens_hi) %>%
  mutate(
   `combined sens` = panel_sens(p, sens, spec),
  )

p2 = ggplot(plot_data, aes(x=sens_lo, y=`combined sens`, colour=factor(ratio))) + 
  geom_line() + 
  geom_abline(colour="grey50") + 
  scale_color_viridis_d(option=  "turbo", name="low:high ratio")+
  geom_abline(colour="black",intercept = 0.8,slope = 0.2) +
  coord_fixed(ylim = c(0.5,1)) +
  xlab("low component sensitivity")+ylab("combined sensitivity")+
  facet_wrap(~"composition")+theme(legend.box = "horizontal")

#p2


## Aggregate vs component sensitivity as with varying specificity -----

# * Assume 10 identically prevalent diseases
# * Fix sensitivity to 75% of each component test; and fix prevalence of all subtypes 
# to 2%
# * Fix sensitivity of components to either 60% or 90%.
# * Varying number of tests with low verus high sensitivity; 


plot_data = tidyr::crossing(
    test = 1:10,
    p = 0.05,
    sens = seq(0.5,1,length.out=101),
    spec = seq(0.75,1,length.out=101)
  ) %>% 
  # group by the thing we are varying
  group_by(spec,sens) %>%
  mutate(
   `combined sens` = panel_sens(p, sens, spec)
  )

p3 = ggplot(plot_data, aes(x=sens, y=`combined sens`, colour=spec, z=spec)) + 
  geom_point(size=0.25) + 
  geom_abline(colour="grey50") + coord_fixed(ylim = c(0.5,1))+
  scale_color_viridis_c(name="specificity")+
  geom_line(data = plot_data %>% filter(spec %in% seq(0.75,1,length.out=5)), mapping = aes(x=sens, y=`combined sens`, group=factor(spec)), colour="grey20", size=0.25)+
  xlab("component sensitivity")+ylab("combined sensitivity")+
  facet_wrap(~"specificity")+theme(legend.box = "horizontal")


## characterisation plot output ----
p = p1+p2+p3+patchwork::plot_layout(ncol=2,guides = "collect")+patchwork::plot_annotation(tag_levels = "A")
save_as(p, here::here("vignettes/latex/s1/fig/component-vs-combined-specificity"),size = std_size$half)

# Simulation validation ----

set.seed(101)

.norm = function(x) {
  if (all(x==0)) x[1] = 1
  x/sum(x)
}

## Setup random set of experiments ----

serotype_prevalence_2 = tibble::tibble(
    # Number of comonetns in group. One of 2:9 in groups of 100
    group = c(1:800,unlist(sapply(seq(800,100,-100), function(x) 1:x))),
    
    # 4400 is 100*2+100*3+100+4+...
    test_id = 1:(sum(100*2:9)),
    
    # FNR ~ 0.2 - sens ~ 80%
    false_neg_rate = rbeta(4400,0.2*100,0.8*100),
    
    
    # FPR ~ 0.025 - spec ~ 0.975
    false_pos_rate = rbeta(4400,0.025*100,0.975*100), 
    
    
  ) %>% inner_join(
    tibble::tibble(
      
      group = 1:800,
      # prevalence ~ 2% per test.
      group_prevalence = rbeta(800,1+1,9+1), # this gives mode of 10%
      n_diseased = 100+sample.int(200,800,replace = TRUE), 
      n_controls = 600+sample.int(400,800,replace = TRUE)
    )
  ) %>%
  mutate( 
    false_neg_diseased = false_neg_rate * n_diseased, 
    false_pos_controls = false_pos_rate * n_controls, 
    # group = factor(sprintf("group %d",group))
  ) %>%
  group_by(group) %>%
  group_modify(function(d,g,...) {
    d %>% mutate(
      group_size = n(),
      distribution = .norm(rpois(group_size,5)*rpois(group_size,1)),
      prevalence = distribution*group_prevalence
    )
  })
  

# ggplot(serotype_prevalence_2, aes(x=false_neg_rate))+geom_density()
# ggplot(serotype_prevalence_2, aes(x=false_pos_rate))+geom_density()


serotype_tests_2 = serotype_prevalence_2 %>%
  group_by(across(everything())) %>%
  group_modify(sim_pop, n=10000, boots=1) %>%
  group_modify(sim_test) 

combined_prev = serotype_tests_2 %>% 
  group_by(across(starts_with("group"))) %>% 
  estimate_panel_performance(test_id)


tmp = combined_prev %>% 
  filter(actual > 0) # & panel_sens_est > 0  & panel_sens_est < 1 )
tmp = tmp %>% mutate(
  panel_actual_prevalence = actual/total,
  prev_est.0.5 = testerror::rogan_gladen(panel_apparent_prevalence, panel_sens_est, panel_spec),
  prev_est_diff.0.5 = panel_actual_prevalence - prev_est.0.5,
  
  sens_est_diff.0.5 = (panel_sens_est - sens_est.0.5),
  sens_est_diff.0.975= (panel_sens_est - sens_est.0.975),
  sens_est_diff.0.025= (panel_sens_est - sens_est.0.025),
  sens_est_calib = panel_sens_est > sens_est.0.025 & panel_sens_est < sens_est.0.975,
  
  sens_diff.0.5 = (panel_sens - sens_est.0.5),
  sens_diff.0.975 = (panel_sens - sens_est.0.975),
  sens_diff.0.025 = (panel_sens - sens_est.0.025),
  sens_calib = panel_sens > sens_est.0.025 & panel_sens < sens_est.0.975,
  
  spec_diff.0.5 = (panel_spec - spec_est.0.5),
  spec_diff.0.975 = (panel_spec - spec_est.0.975),
  spec_diff.0.025 = (panel_spec - spec_est.0.025),
  spec_calib = panel_spec > spec_est.0.025 & panel_spec < spec_est.0.975
) %>% inner_join(serotype_prevalence_2 %>% group_by(group) %>% count(), by="group")


# p1 = ggplot(combined_prev, aes(x=group_prevalence, y=panel_apparent_prevalence, colour=as.factor(group_size)))+geom_point(alpha=0.4)+
#   geom_abline(colour="grey")+coord_fixed(y=c(0,0.6),x=c(0,0.6))+scale_color_viridis_d(option=  "turbo",name="components")+
#   xlab("actual panel prevalence")+ylab("apparent panel prevalence")+facet_wrap(~"components")
# 
# p2 = ggplot(combined_prev, aes(x=group_prevalence, y=panel_apparent_prevalence, colour=sens_est.0.5))+geom_point(alpha=0.4)+
#   geom_abline(colour="grey")+coord_fixed(y=c(0,0.6),x=c(0,0.6))+scale_color_viridis_c(option="viridis", name="sensitivity")+
#   xlab("actual panel prevalence")+ylab("apparent panel prevalence")+facet_wrap(~"panel sensitivity")
# 
# p3 = ggplot(combined_prev, aes(x=group_prevalence, y=panel_apparent_prevalence, colour=spec_est.0.5))+geom_point(alpha=0.4)+
#   geom_abline(colour="grey")+coord_fixed(y=c(0,0.6),x=c(0,0.6))+scale_color_viridis_c(option="magma", name="specificity")+
#   xlab("actual panel prevalence")+ylab("apparent panel prevalence")+facet_wrap(~"panel specificity")
# 
# p = p1+p2+p3+patchwork::plot_layout(ncol=3,guides="collect")+patchwork::plot_annotation(tag_levels = "A")
# save_as(p, here::here("vignettes/latex/s1/fig/simulation-apparent-prevalence"),size = std_size$third)

disagg_prev = serotype_tests_2 %>%
  group_by(group, test_id, sens, spec, prevalence) %>%
  summarise(
    actual = sum(actual),
    test = sum(test),
    total = n()
  ) %>% 
  mutate(
    prev.est = testerror::rogan_gladen(test/total, sens, spec)
  )


p1 = ggplot(disagg_prev, aes(x=prevalence, y=test/total))+geom_abline(colour="red")+
  geom_point(size=0.25,alpha=0.5)+
  coord_fixed(y=c(0,0.6),x=c(0,0.6))+
  xlab("actual prevalence")+ylab("estimated prevalence")+facet_grid("component"~"apparent prevalence")


p2 = ggplot(disagg_prev, aes(x=prevalence, y=prev.est))+geom_abline(colour="red")+
  geom_point(size=0.25,alpha=0.5)+
  coord_fixed(y=c(0,0.6),x=c(0,0.6))+
  xlab("actual prevalence")+ylab("estimated prevalence")+facet_grid("component"~"Rogan-Gladen")




p3 = ggplot(tmp, aes(x=panel_actual_prevalence, y=panel_apparent_prevalence))+geom_abline(colour="red")+
  geom_point(size=0.25,alpha=0.5)+
  coord_fixed(y=c(0,0.6),x=c(0,0.6))+
  xlab("actual prevalence")+ylab("estimated prevalence")+facet_grid("panel"~"apparent prevalence")

p4 = ggplot(tmp,aes(x=panel_actual_prevalence, y=prev_est.0.5))+geom_abline(colour="red")+
  geom_point(size=0.25,alpha=0.5)+
  coord_fixed(xlim=c(0,0.6),ylim=c(0,0.6))+
  xlab("actual prevalence")+ylab("estimated prevalence")+facet_grid("panel"~"Rogan-Gladen (estimated)")

p= p1+no_x()+p2+no_y()+no_x()+p3+p4+no_y()+patchwork::plot_layout(ncol=2)+patchwork::plot_annotation(tag_levels = "A")
save_as(p, here::here("vignettes/latex/s1/fig/qq-prevalence-prediction-v-simulation"),size = std_size$two_third)

## QQ plot ----


p1 = ggplot(tmp,aes(x=panel_sens_est, y=sens_est.0.5, ymin=sens_est.0.025, ymax=sens_est.0.975, colour=sens_est_calib))+geom_abline(colour="red")+
  geom_point(size=0.25,alpha=0.5)+geom_errorbar(size=0.25,alpha=0.05)+
  coord_fixed(xlim=c(0.5,1),ylim=c(0.5,1))+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"sensitivity (theory)")+
  xlab("predicted")+ylab("simulation")

p2 = ggplot(tmp,aes(x=panel_sens, y=sens_est.0.5, ymin=sens_est.0.025, ymax=sens_est.0.975, colour=sens_calib))+geom_abline(colour="red")+
  geom_point(size=0.25,alpha=0.5)+geom_errorbar(size=0.25,alpha=0.05)+
  coord_fixed(xlim=c(0.5,1),ylim=c(0.5,1))+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"sensitivity (estimate)")+
  xlab("predicted")+ylab("simulation")

p3 = ggplot(tmp,aes(x=panel_spec, y=spec_est.0.5, ymin=spec_est.0.025, ymax=spec_est.0.975, colour=spec_calib))+geom_abline(colour="red")+
  geom_point(size=0.25,alpha=0.5)+geom_errorbar(size=0.25,alpha=0.05)+
  coord_fixed(xlim=c(0.5,1),ylim=c(0.5,1))+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"specificity")+
  xlab("predicted")+ylab("simulation")

p = p3+p1+p2+patchwork::plot_layout(ncol=3)+patchwork::plot_annotation(tag_levels = "A")

save_as(p, here::here("vignettes/latex/s1/fig/qq-prediction-v-simulation"),size = std_size$third)


## Absolute errors plot ----


p1 = ggplot(tmp, aes(x=panel_actual_prevalence, y=sens_diff.0.5, ymin=sens_diff.0.025, ymax=sens_diff.0.975, colour=sens_calib))+geom_point(size=0.25)+geom_errorbar(size=0.25,alpha=0.2)+geom_hline(yintercept=0, colour="red")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"sensitivity (theory)")+
  xlab("prevalence")+ylab("error")+ylim(c(-1,1))
p2 = ggplot(tmp, aes(x=panel_actual_prevalence, y=sens_est_diff.0.5, ymin=sens_est_diff.0.025, ymax=sens_est_diff.0.975,colour=sens_est_calib))+geom_point(size=0.25)+geom_errorbar(size=0.25,alpha=0.2)+geom_hline(yintercept=0, colour="red")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"sensitivity (estimate)")+
  xlab("prevalence")+ylab("error")+ylim(c(-1,1))
p3 = ggplot(tmp, aes(x=panel_actual_prevalence, y=spec_diff.0.5, ymin=spec_diff.0.025, ymax=spec_diff.0.975,colour=spec_calib))+geom_point(size=0.25)+geom_errorbar(size=0.25,alpha=0.2)+geom_hline(yintercept=0, colour="red")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"specificity")+
  xlab("prevalence")+ylab("error")+ylim(c(-1,1)*0.02)

p4 = ggplot(tmp, aes(x=panel_spec, y=sens_diff.0.5, ymin=sens_diff.0.025, ymax=sens_diff.0.975, colour=sens_calib))+geom_point(size=0.25)+geom_errorbar(size=0.25,alpha=0.2)+geom_hline(yintercept=0, colour="red")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"sensitivity (theory)")+
  xlab("specificity")+ylab("error")+ylim(c(-1,1))
p5 = ggplot(tmp, aes(x=panel_spec, y=sens_est_diff.0.5, ymin=sens_est_diff.0.025, ymax=sens_est_diff.0.975,colour=sens_est_calib))+geom_point(size=0.25)+geom_errorbar(size=0.25,alpha=0.2)+geom_hline(yintercept=0, colour="red")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"sensitivity (estimate)")+
  xlab("specificity")+ylab("error")+ylim(c(-1,1))
p6 = ggplot(tmp, aes(x=panel_spec, y=spec_diff.0.5, ymin=spec_diff.0.025, ymax=spec_diff.0.975,colour=spec_calib))+geom_point(size=0.25)+geom_errorbar(size=0.25,alpha=0.2)+geom_hline(yintercept=0, colour="red")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"), guide="none")+facet_wrap(~"specificity")+
  xlab("specificity")+ylab("error")+ylim(c(-1,1)*0.02)


p7 = ggplot(tmp, aes(x=factor(n), y=sens_diff.0.5,colour=sens_calib))+ggbeeswarm::geom_beeswarm(size=0.25,corral = "random")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"),guide="none")+geom_hline(yintercept=0, colour="red")+facet_wrap(~"sensitivity (theory)")+
  xlab("components")+ylab("error")+ylim(c(-1,1))
p8 = ggplot(tmp, aes(x=factor(n), y=sens_est_diff.0.5,colour=sens_est_calib))+ggbeeswarm::geom_beeswarm(size=0.25,corral = "random")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"),guide="none")+geom_hline(yintercept=0, colour="red")+facet_wrap(~"sensitivity (estimate)")+ xlab("components")+ylab("error")+ylim(c(-1,1))
p9 = ggplot(tmp, aes(x=factor(n), y=spec_diff.0.5,colour=spec_calib))+ggbeeswarm::geom_beeswarm(size=0.25,corral = "random")+scale_color_manual(values = list("TRUE" = "#000000", "FALSE" = "#FF00FF"),guide="none")+geom_hline(yintercept=0, colour="red")+facet_wrap(~"specificity")+
  xlab("components")+ylab("error")+ylim(c(-1,1)*0.02)

p = p1+no_x()+p4+no_x()+no_y()+p7+no_x()+no_y()+
  p2+p5+no_y()+p8+no_y()+
# p2+no_x()+p5+no_x()+no_y()+p8+no_x()+no_y()+
# p3+p6+no_y()+p9+no_y()+
patchwork::plot_layout(ncol=3)&scale_y_continuous(trans=scales::modulus_trans(-10),limits = c(-1,1), breaks=c(-0.25,-0.1,-0.05,-0.02,0,0.02,0.05,0.1,0.25))


save_as(p, here::here("vignettes/latex/s1/fig/error-prediction-v-simulation"),size = std_size$half)

# p1+p4+no_y()+p7+no_y()+patchwork::plot_layout(ncol=3)&scale_y_continuous(trans=scales::modulus_trans(-10),limits = c(-1,1), breaks=c(-0.5,-0.1,-0.05,-0.02,0,0.02,0.05,0.1,0.5))
# p2+p5+no_y()+p8+no_y()+patchwork::plot_layout(ncol=3)&scale_y_continuous(trans=scales::modulus_trans(-10),limits = c(-1,1), breaks=c(-0.5,-0.1,-0.05,-0.02,0,0.02,0.05,0.1,0.5))
# p3+p6+no_y()+p9+no_y()+patchwork::plot_layout(ncol=3)&scale_y_continuous(trans=scales::modulus_trans(-10),limits = c(-1,1), breaks=c(-0.5,-0.1,-0.05,-0.02,0,0.02,0.05,0.1,0.5))

# serotype_prevalence_2 %>% group_by(group) %>% count() %>% group_by(n) %>% count()

## Calibration ----

p1 = ggplot(tmp %>% filter(panel_actual_prevalence < 0.3), aes(x=panel_actual_prevalence, y=1-sens_calib))+geom_smooth(method = "glm", formula = y ~ splines::ns(x,df=4), method.args=list(family="binomial"), colour="black", size=0.5)+facet_wrap(~"sensitivity (theory)")+
  xlab("prevalence")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))
p2 = ggplot(tmp %>% filter(panel_actual_prevalence < 0.3), aes(x=panel_actual_prevalence, y=1-sens_est_calib))+geom_smooth(method = "glm", formula = y ~ splines::ns(x,df=4), method.args=list(family="binomial"), colour="black", size=0.5)+facet_wrap(~"sensitivity (estimate)")+
  xlab("prevalence")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))
p3 = ggplot(tmp %>% filter(panel_actual_prevalence < 0.3), aes(x=panel_actual_prevalence,  y=1-spec_calib))+geom_smooth(method = "glm", formula = y ~ splines::ns(x,df=4), method.args=list(family="binomial"), colour="black", size=0.5)+facet_wrap(~"specificity")+
  xlab("prevalence")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))


p4 = ggplot(tmp %>% filter(panel_spec > 0.8), aes(x=panel_spec, y=1-sens_calib))+geom_smooth(method = "glm", formula = y ~ splines::ns(x,df=4), method.args=list(family="binomial"), colour="black", size=0.5)+facet_wrap(~"sensitivity (theory)")+
  xlab("specificity")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))
p5 = ggplot(tmp %>% filter(panel_spec > 0.8), aes(x=panel_spec, y=1-sens_est_calib))+geom_smooth(method = "glm", formula = y ~ splines::ns(x,df=4), method.args=list(family="binomial"), colour="black", size=0.5)+facet_wrap(~"sensitivity (estimate)")+
  xlab("specificity")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))
p6 = ggplot(tmp %>% filter(panel_spec > 0.8), aes(x=panel_spec, y=1-spec_calib))+geom_smooth(method = "glm", formula = y ~ splines::ns(x,df=4), method.args=list(family="binomial"), colour="black", size=0.5)+facet_wrap(~"specificity")+
  xlab("specificity")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))

tmp2 = tmp %>% group_by(n) %>% summarise(
  binom_ci_2(sum(sens_calib),n = n(),name = "sens_calib"),
  binom_ci_2(sum(sens_est_calib),n = n(),name = "sens_est_calib"),
  binom_ci_2(sum(spec_calib),n = n(),name = "spec_calib")
) %>% glimpse()

p7 = ggplot(tmp2, aes(x=factor(n), y=1-sens_calib.0.5, ymin=1-sens_calib.0.975, ymax=1-sens_calib.0.025))+geom_errorbar(width=0.6,colour="grey50")+geom_point()+facet_wrap(~"sensitivity (theory)")+
  xlab("components")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))
p8 = ggplot(tmp2, aes(x=factor(n), y=1-sens_est_calib.0.5, ymin=1-sens_est_calib.0.975, ymax=1-sens_est_calib.0.025)) + geom_errorbar(width=0.6,colour="grey50")+geom_point()+facet_wrap(~"sensitivity (estimate)")+
  xlab("components")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))
p9 = ggplot(tmp2, aes(x=factor(n), y=1-spec_calib.0.5, ymin=1-spec_calib.0.975, ymax=1-spec_calib.0.025)) + geom_errorbar(width=0.6,colour="grey50")+geom_point()+facet_wrap(~"specificity")+
  xlab("components")+ylab("disagreement")+geom_hline(yintercept=0.05, colour="red")+coord_cartesian(ylim=c(0,0.2))

#p1+p2+no_y()+p3+no_y()+p4+p5+no_y()+p6+no_y()+p7+p8+no_y()+p9+no_y()+patchwork::plot_layout(ncol=3)

# p1+p4+no_y()+p7+no_y()+patchwork::plot_layout(ncol=3)
# p2+p5+no_y()+p8+no_y()+patchwork::plot_layout(ncol=3)
# p3+p6+no_y()+p9+no_y()+patchwork::plot_layout(ncol=3)

# p = p1+no_x()+p4+no_y()+no_x()+p7+no_y()+no_x()+
# p2+no_x()+p5+no_x()+no_y()+p8+no_x()+no_y()+
# p3+p6+no_y()+p9+no_y()+

## errors & calibration pplot output ----

p = p1+no_x()+p4+no_y()+no_x()+p7+no_y()+no_x()+
p2+p5+no_y()+p8+no_y()+
patchwork::plot_layout(ncol=3)
save_as(p, here::here("vignettes/latex/s1/fig/calibration-prediction-v-simulation"),size = std_size$half)

# serotype_prevalence_2 %>% group_by(group) %>% count() %>% group_by(n) %>% count()

