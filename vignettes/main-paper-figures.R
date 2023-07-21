# "Test error analytical rates"

library(tidyverse)
devtools::load_all("~/Git/ggrrr")
devtools::load_all("~/Git/avoncap")
devtools::load_all("~/Git/testerror")
here::i_am("vignettes/analytical-error.Rmd")
source(here::here("vignettes/vignette-utils.R"))
ggrrr::gg_pedantic()
library(patchwork)

# controls = 800
# samples = 1000
# prev = 0.1
# calculate the false positive distribution
fp_dist = function(prev, samples, control_pos = (1-spec)*controls, control_neg= spec*controls, controls=800, spec=798/800, ...) {
  negatives = round((1-prev)*samples)
  # fp = negatives*(1-sens)
  x= negatives #min(c(fp*10,negatives))
  # adj is done on positives (aka prev)
  #TODO: positives as attribute (i.e. prevalence)
  tibble(
      count = 0:x,
      adj = prev,
      max = x,
      p = extraDistr::dbbinom(0:x, negatives, control_pos, control_neg)
  ) %>%
  mutate(
    p = p/sum(p)
  )
}

empirical_bernoulli = function(dist, p) {
  
  # dist is a vector of probabilities associated with the count, starting at zero
  # we are randomly excluding a proportion of these cases defined by p.
  # when you exclude 1 from the 9 column you get more probabuility in the 8 column
  # when you exclude 2 from the 9 column you get more probabuility in the 7 column
  
  #dist = dbinom(0:5,5,0.1)
  #p=0.001
  
  # e.g. p0,p1,p2,p3,p4,p5
  # think about this in terms of keeping
  # to get 0
  # (c(dbinom(0,0,0.9), dbinom(0,1,0.9), dbinom(0,2,0.9), dbinom(0,3,0.9), dbinom(0,4,0.9), dbinom(0,5,0.9))*dist) %>% sum()
  # to get 1
  # (c(dbinom(1,0,0.9), dbinom(1,1,0.9), dbinom(1,2,0.9), dbinom(1,3,0.9), dbinom(1,4,0.9), dbinom(1,5,0.9))*dist) %>% sum()
  # to get 2
  # (c(dbinom(2,0,0.9), dbinom(2,1,0.9), dbinom(2,2,0.9), dbinom(2,3,0.9), dbinom(3,4,0.9), dbinom(4,5,0.9))*dist) %>% sum()
  n = length(dist)-1
  if (length(dist)==0) return(numeric())
  if (length(dist)==1) return(c(1))
  apply(sapply(0:n, dbinom, 0:n, (1-p))*dist,MARGIN = 2,FUN=sum)
  
  
}

# fp_dist(0.001,1000)

fn_dist = function(prev, samples, spec, disease_pos = sens*diseased, disease_neg= (1-sens)*diseased, diseased=260, sens=0.75,  ...) {
  positives = round(prev*samples)
  # we need the test positive rate
  test_positive_rate = (1-prev)*(1-spec) + prev*sens
  #TODO: positives? as attribute
  # fn = positives*(1-sens)
  x= positives # min(c(fn*10,positives))
  tibble(
      count = 0:x,
      max = x,
      adj = test_positive_rate,
      p = extraDistr::dbbinom(0:x, positives, disease_neg, disease_pos)
  ) %>%
  mutate(
    p = p/sum(p)
  )
}

# TODO: think this through:
# When combining FP or FN distributions we have to adjust for the fact that 
# some FP on the LHS will be FN, TP (=condition positives) on the RHS and vice versa.
# also some FP on the LHS will be FP on the RHS.
# This essentially reduces the effective count on each side that we are combining.
# sharpening the distribution and shifting the probability towards zero counts.
# 
# basically I'm trying to adjust for the situation where a false positive 
# in one test is also an actual positive in another test. The true positive
# trumping the false positive (or really the false positive becomes reinterpreted as a true positive)


# e.g. p(count=1) = p(count=1) + p(count=2 and 1 collision) + p(count=3 and 2 collisions) etc.
# p collision is just prevalence in other part of combination. If high lots of collisions if low not very many.

# FP only are counted as a FP in the combination when combined with a TN or FP.
# FP+FP=FP; FP+FN=TP (for the wrong reason); FP+TP=TP; FP+TN=FP
# FP colliding with FN or TP need to be excluded completely (actual positives = prev)
# FP colliding with FP need to be counted only once

# This is almost the same for FN. They only remain a FN and are counted in combination when combined with a TN.
# FN+FN = FN; FN+TP = TP; FN+TN = FN; FN+FP = TP (although true for the wrong reason)
# FN colliding with FP or TP need to be excluded completely (test positives = prev*sens)
# FN colliding with FN need to be counted only once

# Both these scenarios are the same as exclude all collisions of LHS with not(TN) 
# I.e. ((1-prev)*spec)
# and add back in FP or FN although this is not particularly helpful I think.

# therefore we need to check for collisions with TP and FP (i.e. test positives) and 
# I.e. P(not FN or a TN) which is a P(test positive) = FP+TP = (1-prev)*(1-spec)+prev*sens

dists_sum = function(a_data, b_data, samples = 1000) {
  
  a_frac = a_data %>% summarise(f = sum(count*p)) %>% pull(f)/samples
  b_frac = b_data %>% summarise(f = sum(count*p)) %>% pull(f)/samples
  # fraction from a & b that overlap (i.e. the FP rate on LHS and RHS)
  a_b_adj = a_frac * b_frac
  
  # b_adj is the probability that a single observation in B is also in A as a 
  # positive or a false positive (for false positives) or as 
  # a test positive or false negative in B
  b_adj = unique(a_data$adj) + a_frac
  
  # browser()
  
  b_data_2 = b_data %>% arrange(count) %>% mutate(
    
    # adjust the probability distribution based on the possibility of collision
    # of the new FP and existing TP/FP by adjusting by b_adj
    p = empirical_bernoulli(p, b_adj)
      
  )
  
  # having performed adjustment on b we can combine the new distributions as a sum assuming
  # no collision.
  a_data %>% 
    #filter(p>0.0000001) %>%
    # cross join
    inner_join(
      b_data_2, #%>% filter(p>0.0000001), 
      by=character(), suffix = c(".lhs",".rhs")
    ) %>%
    mutate(
      count = count.lhs + count.rhs,
      p = p.lhs * p.rhs
    ) %>%
    group_by(count) %>%
    summarise(p = sum(p),.groups = "drop") %>%
    filter(count <= samples) %>%
    mutate(
      p = p/sum(p),
      # the overlap / collision adjustment term. This is all positives for the FP
      # distribution and all true positives for the FN distribution.
      # The positives (FP) can be combined; the test positives (FN) can also be combined.
      adj = 1-(1-unique(a_data$adj))*(1-unique(b_data$adj))) %>%
    arrange(count)
}

dists_diff = function(a_data, b_data) {
  a_data %>% 
    # filter(p>0.0000001) %>%
    # cross join
    cross_join(b_data, # %>% filter(p>0.0000001), 
               suffix = c(".lhs",".rhs")) %>%
    mutate(
      count = count.lhs - count.rhs,
      p = p.lhs * p.rhs
    ) %>%
    group_by(count) %>%
    summarise(p = sum(p)) %>%
    arrange(count)
}

# combined_fp = lapply(1:20, function(x) fp_dist(0.005,1000)) %>% reduce(dists_sum)
# combined_fn = lapply(1:20, function(x) fn_dist(0.005,1000,0.9975)) %>% reduce(dists_sum)
# 
# lapply(1:2, function(x) fp_dist(0.2,1000,control_pos = 10,control_neg = 90)) %>% reduce(dists_sum)
# 
# sum(combined_fp$p)
# sum(combined_fn$p)
# 
# example_plot(prevalence = combined_p, samples = 1000, fp=combined_fp, fn=combined_fn)

discrete_quantile = function(dist, q) {
  x = dist %>% arrange(count) %>% mutate(c = cumsum(p))
  tmp = sapply(q, function(q1) min(which(x$c > q1)))
  tmp = x$count[tmp]
  names(tmp) = q
  return(tmp)
}

specificity_est = function(fp_dist, prevalence, samples, ci=NULL) {
  n1 = (1-prevalence)*samples
  if (is.null(ci)) {
    fp1 = fp_dist %>% summarise(x = sum(count*p)) %>% pull(x)
  } else {
    fp1 = discrete_quantile(fp_dist, ci)
  }
  spec = rev(1-fp1/n1)
  names(spec) = ci
  return(spec)
}

sensitivity_est = function(fn_dist, prevalence, samples, ci=NULL) {
  p1 = prevalence*samples
  if (is.null(ci)) {
    fn1 = fn_dist %>% summarise(x = sum(count*p)) %>% pull(x)
  } else {
    fn1 = discrete_quantile(fn_dist, ci)
  }
  sens = rev(1-fn1/p1)
  names(sens) = ci
  return(sens)
}


summary_test = function(
    prevalence, samples, spec = 0.9975,
    ...,
    fp = fp_dist(prevalence, samples, ...),
    fn = fn_dist(prevalence, samples, spec, ...)) 
  {
    p1 = prevalence*samples
    spec = specificity_est(fp,prevalence,samples)
    sens = sensitivity_est(fn,prevalence,samples)
    fp1 = fp %>% summarise(f = sum(count*p)) %>% pull(f)
    fn1 = fn %>% summarise(f = sum(count*p)) %>% pull(f)
    
    comb = dists_diff(fp, fn) 
    comb_f = comb %>% summarise(f = sum(count*p)) %>% pull(f)
    
  tribble(
    
    
    ~characteristic, ~value,
    "N", sprintf("%d",samples),
    # "Prevalence", sprintf("%1.4g%%", prevalence*100),
    # "Sensitivity", sprintf("%1.4f",sens),
    # "Specificity", sprintf("%1.4f",spec),
    "E(FP)",  sprintf("%1.2f [%d \u2014 %d]",
                fp1,
                discrete_quantile(fp, 0.025),
                discrete_quantile(fp, 0.975)
    ),
    "E(FN)",  sprintf("%1.2f [%d \u2014 %d]",
                fn1,
                discrete_quantile(fn, 0.025),
                discrete_quantile(fn, 0.975)
    ),
    "Condition pos", sprintf("%d", floor(p1)),
    "E(Test pos)", sprintf("%1.2f [%d \u2014 %d]",
                p1+comb_f,
                round(p1+discrete_quantile(comb, 0.025)),
                round(p1+discrete_quantile(comb, 0.975))
    )
  )
}


example_plot = function(
    prevalence = 0.1, 
    samples = 1000,
    spec = 0.9975,
    ...,
    fp = fp_dist(prevalence, samples, ...),
    fn = fn_dist(prevalence, samples, spec, ...),
    xlim = c(NA,NA)
) {
  
  p1 = round(prevalence*samples)
  fp1 = fp %>% summarise(f = sum(count*p)) %>% pull(f)
  fn1 = fn %>% summarise(f = sum(count*p)) %>% pull(f)
  
  title = sprintf("p=%1.2g%%",prevalence*100)
  
  tmp3 = bind_rows(
     fp %>% mutate(style = "fp", x=prevalence*samples+count),
     fn %>% mutate(style = "fn", x=prevalence*samples-count, p = -p)
    ) %>% 
    mutate(prevalence = title) %>%
    filter(abs(p) > 0.0001)
  
  ylims = c(-1,1)*max(abs(tmp3$p))
  
  comb = dists_diff(fp, fn) %>% 
    mutate(prevalence = title) %>%
    filter(abs(p) > 0.0001)
  
  comb_f = comb %>% summarise(f = sum(count*p)) %>% pull(f)
  # comb_ylims = c(0,max(conb$p))
  
  list(
    
    ggplot(tmp3, aes(x=x, y=p, fill = style, colour = style))+geom_bar(stat="identity")+
      geom_vline(xintercept = p1, colour="black")+
      # geom_text(x=p1, y=Inf, label="real positives", vjust=1, hjust=1, angle=90 )+
      geom_vline(xintercept = p1+fp1, colour="blue")+
      geom_vline(xintercept = p1-fn1, colour="red")+
      geom_hline(yintercept = 0,colour= "black")+
      guides(fill=guide_none(), colour = guide_none())+
      ylab("probability")+
      xlab("count")+
      facet_wrap(~prevalence) +
      coord_cartesian(ylim=ylims, xlim=xlim),
    
    ggplot(comb, aes(x=count+p1, y=p))+geom_bar(stat="identity",fill="grey50",colour="grey50")+
      geom_vline(xintercept = p1, colour="black")+
      geom_vline(xintercept = p1+comb_f, colour="magenta")+
      facet_wrap(~prevalence) +
      ylab("probability")+
      xlab("count")+
      coord_cartesian(xlim=xlim),
    
    ggrrr::gg_simple_table(summary_test(prevalence = prevalence,samples = samples,fp=fp,fn=fn),pts = 5)
  )
}



# patchwork::wrap_plots(example_plot(0.16), nrow=1) %>% ggrrr::gg_save_as(tempfile(fileext = ".png"), size = std_size$quarter)

# patchwork::wrap_plots(example_plot(0.16), nrow=1) #%>% plot_to_google("Combined panel tests",1, size = std_size$third)
p = patchwork::wrap_plots(example_plot(0.02), nrow=1) /
patchwork::wrap_plots(example_plot(0.005), nrow=1) /
# patchwork::wrap_plots(example_plot(0.002), nrow=1) /
patchwork::wrap_plots(example_plot(0), nrow=1) +
  patchwork::plot_annotation(tag_levels = "A")

# p %>% ggrrr::gg_save_as(tempfile(fileext = ".png"), size = std_size$two_third)
p %>% save_as(here::here("vignettes/latex/main/fig/low-prevalence-sensitivity-specificity"), size = std_size$half)


# tablular views of combination scenario ----
              
# combined_p = 1-(1-0.005)^20
# # this is an overestimate as it fails to account for false positves in one test
# # that are true positives in another test.
# 
# combined_fp = lapply(1:20, function(x) fp_dist(0.005,1000)) %>% reduce(dists_sum)
# combined_fn = lapply(1:20, function(x) fn_dist(0.005,1000,0.9975)) %>% reduce(dists_sum)
# 
# combined_fp %>% specificity_est(combined_p, 1000, ci = c(0.025,0.975))
# combined_fn %>% sensitivity_est(combined_p, 1000, ci = c(0.025,0.975))
# 
# tmp = example_plot(prevalence = combined_p, samples = 1000, fp=combined_fp, fn=combined_fn)
# 
# p = 
#   patchwork::wrap_plots(example_plot(0.005), nrow=1) /
#   patchwork::wrap_plots(tmp, nrow=1)+
#   patchwork::plot_annotation(tag_levels = "A")
# 
# p
# # p %>% ggrrr::gg_save_as(tempfile(fileext = ".png"), size = std_size$half)
# # p %>% plot_to_google("Combined panel tests", 3, size = std_size$half)
# 
# t = bind_rows(
#   summary_test(prevalence = 0.16, samples = 1000) %>% mutate(configuration = "1 x 16%"),
#   # summary_test(prevalence = 0.005, samples = 1000) %>% mutate(configuration = "1 x 0.5%"),
#   summary_test(prevalence = combined_p, samples = 1000, fp = combined_fp, fn = combined_fn) %>% mutate(configuration = "20 x 0.5%")
# ) %>% ggrrr::hux_tidy(rowGroupVars = vars(characteristic), colGroupVars = vars(configuration))
# 
# t
# # t %>% table_to_google("Combined panel tests", 1)
# 
# 
# # plot for combination of 4 groups as per the
# 
# pv = c(0,0.001,0.0025,0.025) %>% lapply(rep,5) %>% unlist()
# 
# combined_summary = function(pv, name, spec=0.9975) {
#   combined_p = 1-prod((1-pv))
#   combined_fp = lapply(pv, function(x) fp_dist(x,1000)) %>% reduce(dists_sum)
#   combined_fn = lapply(pv, function(x) fn_dist(x,1000,spec)) %>% reduce(dists_sum)
#   summary_test(prevalence = combined_p, samples = 1000, fp = combined_fp, fn = combined_fn) %>% mutate(configuration = name)
# }
# 
# 
# 
# t2 = bind_rows(
#   combined_summary(rep(0,15), "15x0%"),
#   combined_summary(rep(0.03,5), "5x3%"),
#   combined_summary(c(rep(0,15),rep(0.03,5)), "15x0% + 5x3%")
# ) %>% ggrrr::hux_tidy(rowGroupVars = vars(characteristic), colGroupVars = vars(configuration))
# t2 # %>% table_to_google("Combined panel tests", 2)
# 
# t3 = bind_rows(
#   combined_summary(rep(0,13), "13x0%"),
#   combined_summary(rep(0.05,1), "1x5%"),
#   combined_summary(c(rep(0,13),rep(0.05,1)), "13x0% + 1x5%")
# ) %>% ggrrr::hux_tidy(rowGroupVars = vars(characteristic), colGroupVars = vars(configuration))
# t3 # %>% table_to_google("Combined panel tests", 3)

# The main paper

dist = ipd_distribution("PCV20") %>% 
  mutate(prev =  distribute(distribution, 0.1)) %>%
  pull(prev) %>%
  unique()


p1 = apparent_prevalence_plot(
  p=dist, 
  cols = scales::brewer_pal(palette="Dark2")(length(p)), 
  lim=c(0,0.05),
  sens=0.80,
  spec = 0.9975,
  bottom_right = NULL,
  top_left = "IPD distribution\n20 components\ncomponent sens: 0.8\ncomponent spec: 0.9975",
  annotate_cols = list(top_left="black")
)+facet_wrap(~"component")

plot_layer = ipd_distribution("PCV20") %>% 
  cross_join(tibble(panel_prev = seq(0,1,length.out=101))) %>%
  group_by(panel_prev) %>%
  mutate(
    prev =  distribute(distribution, panel_prev),
    spec = 0.9975,
    sens = 0.8
  ) %>%
  summarise(
    panel_sens = panel_sens(prev, sens = sens, spec = spec),
    panel_spec = panel_spec(spec)
  ) %>% mutate(
    panel_apparent_prev = apparent_prevalence(panel_prev, panel_sens, panel_spec)
  ) %>%
  glimpse()

pl = plot_layer %>% 
  filter(panel_prev == 0.1) %>%
  mutate(xaxis = 0, yaxis = 0, xlabel="0.1", ylabel=sprintf("%1.3f",panel_apparent_prev))

pl2 = plot_layer %>% 
  filter(panel_prev > panel_apparent_prev) %>%
  filter(panel_prev == min(panel_prev)) %>%
  mutate(label=sprintf("(%1.3f, %1.3f)",panel_prev,panel_prev))

p2 = ggplot(plot_layer, aes(x=panel_prev, y=panel_apparent_prev))+
  coord_fixed(xlim=c(0,1),ylim=c(0,1), clip="off")+
  geom_abline(colour="grey")+
  geom_line()+
  geom_segment(aes(x=panel_prev, y=yaxis, xend=panel_prev, yend=panel_apparent_prev), data=pl, linetype="dashed")+
  geom_text(aes(x=panel_prev, y=yaxis, label=xlabel),data=pl, hjust=1.1, vjust=0.5, angle=90)+
  geom_segment(aes(x=xaxis, y=panel_apparent_prev, xend=panel_prev, yend=panel_apparent_prev), data=pl, linetype="dashed")+
  geom_text(aes(x=xaxis, y=panel_apparent_prev, label=ylabel),data=pl, hjust=1.1, align="right")+
  geom_point(aes(x=panel_prev, y=panel_apparent_prev),data=pl, size=0.5)+
  scale_x_continuous(breaks = c(0,1),expand = c(0, 0))+
  scale_y_continuous(breaks = c(0,1),expand = c(0, 0))+
  xlab("True panel prevalence")+
  ylab("E(Apparent panel prevalence)")+
  geom_point(aes(x=panel_prev, y=panel_apparent_prev), data=pl2, size=0.5)+
  geom_text(aes(x=panel_prev, y=panel_apparent_prev, label=label), vjust = 1.1, hjust = -0.1, data=pl2)+
  #annotate("text",x=0.1, y=0.9, vjust="inward", hjust="inward", label="IPD distribution\n20 components\ncomponent spec: 0.9975\ncomponent sens: 0.8",size = 6/ggplot2:::.pt)+
  theme(
    plot.margin = unit(c(1,1,1,1), "lines"),
    #axis.title.x = element_text(margin=margin(t=2, unit="lines")),
    #axis.title.y = element_text(margin=margin(r=2, unit="lines"))
    axis.title.x = element_text(hjust = 1),
    axis.title.y = element_text(hjust = 1)
  ) + facet_wrap(~"panel")

p = p1+p2+patchwork::plot_layout(ncol=2)+patchwork::plot_annotation(tag_levels = "A")

save_as(p, here::here("vignettes/latex/main/fig/true-apparent-prevalence-component-panels"),size = std_size$third)


# Nottingham scenario ----

specs = c(0.000001, 1-1.01^-(1:1000), 0.9999999)

p_spec = do_scenario(
    spec = specs, 
    n_controls = 800,
    group_prevalence = 0.1, samples=FALSE,
    group = sprintf("specificity: %1.3f", specs )
  ) %>% 
  unnest(pcv_group) %>%
  group_by(group, group_prevalence, panel_id, false_pos_rate) %>%
  summarise(
    panel_prevalence = testerror::panel_prevalence(prevalence),
    panel_sensitivity = testerror::panel_sens(prevalence, 1-false_neg_rate, 1-false_pos_rate),
    panel_specificity = testerror::panel_spec(1-false_pos_rate),
  ) %>%
  mutate(
    panel_apparent_prevalence = testerror::apparent_prevalence(panel_prevalence, panel_sensitivity, panel_specificity),
    panel_id = factor(panel_id, c("PCV7","PCV13","PCV15","PCV20"))
  ) %>%
  glimpse()

ggplot(p_spec, aes(x=1-false_pos_rate, y= panel_prevalence, colour=panel_id))+geom_line()+scale_x_continuous(trans="logit", breaks=c(0.5, 0.75,0.9, 0.95, 0.975, 0.99, 0.999, 0.9999))+
  coord_cartesian(xlim=c(0.4, 0.9999))

ggplot(p_spec, aes(x=1-false_pos_rate, y= panel_apparent_prevalence, colour=panel_id))+geom_line()+scale_x_continuous(trans="logit", breaks=c(0.5, 0.75,0.9, 0.95, 0.975, 0.99, 0.999, 0.9999))+
  coord_cartesian(xlim=c(0.4, 0.9999))+
  geom_line(aes(x=1-false_pos_rate, y= panel_prevalence, colour=panel_id), linetype="dashed")+
  xlab("component specificity")+
  ylab("panel apparent prevalence")
  

# sensitivity

senses = seq(0.000001, 0.9999999, length.out= 1001)

p_spec = do_scenario(
  spec = 0.9975, 
  sens = senses,
  n_controls = 800,
  n_diseased = 26,
  group_prevalence = 0.1, samples=FALSE,
  group = sprintf("sensitivity: %1.3f", senses )
) %>% 
  unnest(pcv_group) %>%
  group_by(group, group_prevalence, panel_id, false_neg_rate) %>%
  summarise(
    panel_prevalence = testerror::panel_prevalence(prevalence),
    panel_sensitivity = testerror::panel_sens(prevalence, 1-false_neg_rate, 1-false_pos_rate),
    panel_specificity = testerror::panel_spec(1-false_pos_rate),
  ) %>%
  mutate(
    panel_apparent_prevalence = testerror::apparent_prevalence(panel_prevalence, panel_sensitivity, panel_specificity),
    panel_id = factor(panel_id, c("PCV7","PCV13","PCV15","PCV20"))
  ) %>%
  glimpse()


ggplot(p_spec, aes(x=1-false_neg_rate, y= panel_apparent_prevalence, colour=panel_id))+geom_line()+
  # scale_x_continuous(trans="logit", breaks=c(0.5, 0.75,0.9, 0.95, 0.975, 0.99, 0.999, 0.9999))+
  coord_cartesian(xlim=c(0, 0.9999))+
  geom_line(aes(x=1-false_neg_rate, y= panel_prevalence, colour=panel_id), linetype="dashed")+
  xlab("component sensitivity")+
  ylab("panel apparent prevalence")

## Surface plot

senses2 = senses[10*0:100+1]
specs2 = specs[10*0:100+1]

matching_sens = function(toMatch) sapply(toMatch, function(x) which(abs(senses2-x) <= .Machine$double.eps))
matching_spec = function(toMatch) sapply(toMatch, function(x) which(abs(specs2-x) <= .Machine$double.eps))

## assuming prevalence is 10%

uncache = memoise::drop_cache(do_scenario)

uncache(    sens = senses2, 
            spec = specs2, 
            n_controls = 800,
            n_diseased = 26,
            group_prevalence = 0.1, 
            samples=FALSE,
            group = sprintf("grp: %d", 1:10201))

p_surface = do_scenario(
    sens = senses2, 
    spec = specs2, 
    n_controls = 800,
    n_diseased = 26,
    group_prevalence = 0.1, 
    samples=FALSE,
    group = sprintf("grp: %d", 1:10201)
) %>% 
  unnest(pcv_group) %>%
  group_by(group, group_prevalence, panel_id, false_neg_rate, false_pos_rate) %>%
  summarise(
    panel_prevalence = testerror::panel_prevalence(prevalence),
    panel_sensitivity = testerror::panel_sens(prevalence, 1-false_neg_rate, 1-false_pos_rate),
    panel_specificity = testerror::panel_spec(1-false_pos_rate),
  ) %>%
  mutate(
    panel_apparent_prevalence = testerror::apparent_prevalence(panel_prevalence, panel_sensitivity, panel_specificity),
    panel_id = factor(panel_id, c("PCV7","PCV13","PCV15","PCV20")),
    sens_tmp = 1-false_neg_rate, 
    spec_tmp = 1-false_pos_rate
  )  

p_surface = p_surface %>% mutate(
    sens_low = lag(senses2,default = 0)[matching_sens(sens_tmp)],
    sens_high = lead(senses2,default = 1)[matching_sens(sens_tmp)],
    spec_low = lag(specs2,default = 0)[matching_spec(spec_tmp)],
    spec_high = lead(specs2,default = 1)[matching_spec(spec_tmp)],
    sens_low = (sens_low+sens_tmp)/2,
    sens_high = (sens_high+sens_tmp)/2,
    spec_low = (spec_low+spec_tmp)/2,
    spec_high = (spec_high+spec_tmp)/2,
    prev_diff = panel_apparent_prevalence-panel_prevalence,
    rel_diff = panel_apparent_prevalence/panel_prevalence
  ) %>%
  glimpse()

p1 = ggplot(p_surface, aes(xmin=spec_low, xmax=spec_high, ymin=sens_low, ymax=sens_high, fill=rel_diff))+
  geom_rect()+
  coord_cartesian(xlim=c(0.95, 0.9999),ylim=c(0, 0.9999))+
  xlab("component specificity")+
  ylab("component sensitivity")+
  scale_fill_gradientn(limits=c(0,4), colours=c("blue","white","red"), values = scales::rescale(c(0,1,4)), name="relative error", oob=scales::squish, guide = "none") +
  geom_vline(xintercept = c(0.98,0.9975),linetype="dotted",colour="grey20")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(0), inherit.aes = FALSE, colour="black")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(-0.01), inherit.aes = FALSE, colour="#8080FF")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(0.01), inherit.aes = FALSE, colour="#FF8080")+
  metR::geom_contour2(aes(x=spec_tmp, y = sens_tmp, z=rel_diff, label=sprintf("\u00d7%1.3g",after_stat(level))), 
                      breaks=c(0.5, 0.8,1,1.25,2,4,8,16), 
                      inherit.aes = FALSE, 
                      label.placer = metR::label_placer_fraction(frac = 0.3),
                      skip=0,
                      label_size = 6/ggplot2:::.pt,
                      margin = grid::unit(c(2, 2, 2, 2), "pt")
  )+
  geom_point(aes(x=0.9975, y=0.8),colour="blue",size=2, shape=4)+
  # scale_fill_gradient2(low = "blue",mid = "white",high = "red") +
  scale_x_continuous(trans="logit", breaks=c(0.5, 0.75,0.9, 0.95, 0.98, 0.99, 0.9975, 0.9999)) +
  
  facet_wrap(~panel_id)

p1 %>% save_as(file=here::here("vignettes/latex/s3/fig/relative-error-by-sensitivity-specificity"),size = std_size$half)

## Apparent VE:


p_surface_2 = do_scenario(
  sens = senses2, 
  spec = specs2, 
  n_controls = 800,
  n_diseased = 26,
  group_prevalence = 0.05, 
  samples=FALSE,
  group = sprintf("grp: %d", 1:10201)
) %>% 
  unnest(pcv_group) %>%
  group_by(group, group_prevalence, panel_id, false_neg_rate, false_pos_rate) %>%
  summarise(
    panel_prevalence = testerror::panel_prevalence(prevalence),
    panel_sensitivity = testerror::panel_sens(prevalence, 1-false_neg_rate, 1-false_pos_rate),
    panel_specificity = testerror::panel_spec(1-false_pos_rate),
  ) %>%
  mutate(
    panel_apparent_prevalence = testerror::apparent_prevalence(panel_prevalence, panel_sensitivity, panel_specificity),
    panel_id = factor(panel_id, c("PCV7","PCV13","PCV15","PCV20")),
    sens_tmp = 1-false_neg_rate, 
    spec_tmp = 1-false_pos_rate
  )  

# Lower prevalence ----

p_surface_2 = p_surface_2 %>% mutate(
  sens_low = lag(senses2,default = 0)[matching_sens(sens_tmp)],
  sens_high = lead(senses2,default = 1)[matching_sens(sens_tmp)],
  spec_low = lag(specs2,default = 0)[matching_spec(spec_tmp)],
  spec_high = lead(specs2,default = 1)[matching_spec(spec_tmp)],
  sens_low = (sens_low+sens_tmp)/2,
  sens_high = (sens_high+sens_tmp)/2,
  spec_low = (spec_low+spec_tmp)/2,
  spec_high = (spec_high+spec_tmp)/2,
  prev_diff = panel_apparent_prevalence-panel_prevalence,
  rel_diff = panel_apparent_prevalence/panel_prevalence
) %>%
  glimpse()

p1b = ggplot(p_surface_2, aes(xmin=spec_low, xmax=spec_high, ymin=sens_low, ymax=sens_high, fill=rel_diff))+
  geom_rect()+
  coord_cartesian(xlim=c(0.95, 0.9999),ylim=c(0, 0.9999))+
  xlab("component specificity")+
  ylab("component sensitivity")+
  scale_fill_gradientn(limits=c(0,4), colours=c("blue","white","red"), values = scales::rescale(c(0,1,4)), name="relative error", oob=scales::squish, guide = "none") +
  geom_vline(xintercept = c(0.98,0.9975),linetype="dotted",colour="grey20")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(0), inherit.aes = FALSE, colour="black")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(-0.01), inherit.aes = FALSE, colour="#8080FF")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(0.01), inherit.aes = FALSE, colour="#FF8080")+
  metR::geom_contour2(aes(x=spec_tmp, y = sens_tmp, z=rel_diff, label=sprintf("\u00d7%1.3g",after_stat(level))), 
                      breaks=c(0.5, 0.8,1,1.25,2,4,8,16), 
                      inherit.aes = FALSE, 
                      label.placer = metR::label_placer_fraction(frac = 0.3),
                      skip=0,
                      label_size = 6/ggplot2:::.pt,
                      margin = grid::unit(c(2, 2, 2, 2), "pt")
  )+
  geom_point(aes(x=0.9975, y=0.8),colour="blue",size=2, shape=4)+
  # scale_fill_gradient2(low = "blue",mid = "white",high = "red") +
  scale_x_continuous(trans="logit", breaks=c(0.5, 0.75,0.9, 0.95, 0.98, 0.99, 0.9975, 0.9999)) +
  
  facet_wrap(~panel_id)

# VE ----

p1b %>% save_as(file=here::here("vignettes/latex/s3/fig/relative-error-by-sensitivity-specificity_5_percent"),size = std_size$half)

ve_surface = p_surface %>% 
inner_join(
  p_surface_2 %>% select(group, panel_id, panel_prevalence,panel_apparent_prevalence),
  by = join_by(group, panel_id),
  suffix = c(".hi",".lo")
) %>% mutate(
  true_ve = 1-panel_prevalence.lo/panel_prevalence.hi,
  apparent_ve = 1-panel_apparent_prevalence.lo/panel_apparent_prevalence.hi,
  abs_diff = apparent_ve - true_ve,
  rel_diff = apparent_ve/true_ve,
  true_preventable = panel_prevalence.hi-panel_prevalence.lo,
  apparent_preventable = panel_apparent_prevalence.hi-panel_apparent_prevalence.lo,
  abs_diff_prevent = apparent_preventable - true_preventable,
  rel_diff_prevent = apparent_preventable/true_preventable,
) %>% glimpse()
  
# VE ----

p2 = ggplot(ve_surface, aes(xmin=spec_low, xmax=spec_high, ymin=sens_low, ymax=sens_high, fill=rel_diff))+
  geom_rect()+
  coord_cartesian(xlim=c(0.95, 0.9999),ylim=c(0, 0.9999))+
  xlab("component specificity")+
  ylab("component sensitivity")+
  scale_fill_gradient(limits=c(0,1), low = "violet" , high = "white", name="relative error", oob=scales::squish, guide = "none") +
  geom_vline(xintercept = c(0.98,0.9975),linetype="dotted",colour="grey20")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(0), inherit.aes = FALSE, colour="black")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(-0.01), inherit.aes = FALSE, colour="#8080FF")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(0.01), inherit.aes = FALSE, colour="#FF8080")+
  metR::geom_contour2(aes(x=spec_tmp, y = sens_tmp, z=rel_diff, label=sprintf("\u00d7%1.3g",after_stat(level))), 
                      breaks=c(0.95, 0.9, 0.75, 0.5, 0.25, 0.1), 
                      inherit.aes = FALSE, 
                      label.placer = metR::label_placer_fraction(frac = 0.70),
                      skip=0,
                      label_size = 6/ggplot2:::.pt,
                      margin = grid::unit(c(2, 2, 2, 2), "pt")
  )+
  geom_point(aes(x=0.9975, y=0.8),colour="blue",size=2, shape=4)+
  # scale_fill_gradient2(low = "blue",mid = "white",high = "red") +
  scale_x_continuous(trans="logit", breaks=c(0.5, 0.75,0.9, 0.95, 0.98, 0.99, 0.9975, 0.9999)) +
  facet_wrap(~panel_id)

p2 %>% save_as(file=here::here("vignettes/latex/s3/fig/ve-error-by-sensitivity-specificity"),size = std_size$half)


# Preventable disease

p3 = ggplot(ve_surface, aes(xmin=spec_low, xmax=spec_high, ymin=sens_low, ymax=sens_high, fill=abs_diff_prevent))+
  geom_rect()+
  coord_cartesian(xlim=c(0.95, 0.9999),ylim=c(0, 0.9999))+
  xlab("component specificity")+
  ylab("component sensitivity")+
  scale_fill_gradient(limits=c(-0.05,0), low = "orange" , high = "white", name="relative error", oob=scales::squish, guide = "none") +
  geom_vline(xintercept = c(0.98,0.9975),linetype="dotted",colour="grey20")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(0), inherit.aes = FALSE, colour="black")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(-0.01), inherit.aes = FALSE, colour="#8080FF")+
  # geom_contour(aes(x=spec_tmp, y = sens_tmp, z=prev_diff),breaks = c(0.01), inherit.aes = FALSE, colour="#FF8080")+
  metR::geom_contour2(aes(x=spec_tmp, y = sens_tmp, z=abs_diff_prevent, label=sprintf("%+1.1f%%",after_stat(level)*100)), 
                      breaks=c(-0.05,-0.04,-0.03,-0.02,-0.01,-0.001,0), 
                      inherit.aes = FALSE, 
                      label.placer = metR::label_placer_fraction(frac = 0.3),
                      skip=0,
                      label_size = 6/ggplot2:::.pt,
                      margin = grid::unit(c(2, 2, 2, 2), "pt")
  )+
  # scale_fill_gradient2(low = "blue",mid = "white",high = "red") +
  geom_point(aes(x=0.9975, y=0.8),colour="blue",size=2, shape=4)+
  scale_x_continuous(trans="logit", breaks=c(0.5, 0.75,0.9, 0.95, 0.98, 0.99, 0.9975, 0.9999)) +
  facet_wrap(~panel_id)

p3 %>% save_as(file=here::here("vignettes/latex/s3/fig/preventable-error-by-sensitivity-specificity"),size = std_size$half)

# restict to PCV20 for main paper ----

p1a = p1                                                                
p1a$data <- p1a$data %>% filter(panel_id=="PCV20")
p1ba = p1b                                                                
p1ba$data <- p1ba$data %>% filter(panel_id=="PCV20")
p2a = p2                                                                
p2a$data <- p2a$data %>% filter(panel_id=="PCV20")
p3a = p3                                                                
p3a$data <- p3a$data %>% filter(panel_id=="PCV20")

p4 = p1a + facet_wrap(~"10% prevalence relative error") + no_x() +
  p1ba + facet_wrap(~"5% prevalence relative error") + no_y() + no_x() +
  p2a + facet_wrap(~"50% VE relative error") + 
  p3a + facet_wrap(~"5% preventable disease absolute error") + no_y() +
  plot_layout(nrow=2)
                                                                  
p4 %>% save_as(file=here::here("vignettes/latex/s3/fig/impact-error-by-sensitivity-specificity-4-panel"),size = std_size$half)


(p1a+facet_wrap(~"10% prevalence relative error")) %>% save_as(file=here::here("vignettes/latex/main/fig/impact-error-by-sensitivity-specificity"),size = std_size$half)

