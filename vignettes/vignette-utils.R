## Simulation ----



ipd_distribution = function(pcv_groups = c("PCV7","PCV13","PCV15","PCV20")) {
  ipd = readRDS(here::here("vignettes/ipd-serotype-distribution.rds"))
  tmp = avoncap::serotype_data$map %>% pivot_longer(cols = -serotype, names_to = "pneumo.group") %>%
    filter(pneumo.group %in% pcv_groups & value==TRUE) %>% 
    select(-value) %>%
    rename(panel_id = pneumo.group) %>%
    nest(pcv_group = panel_id) %>%
    rename(pneumo.phe_serotype = serotype) %>%
    left_join(ipd$by_serotype %>% select(pneumo.phe_serotype, x) %>% dtrackr::untrack(), by = join_by(pneumo.phe_serotype)) %>% 
    mutate(
      x= ifelse(is.na(x),0,x),
      n= sum(x),
      distribution = x/n
    )
  return(tmp)
}


  
  # tmp_ipd = avoncap::serotype_data$xr %>% filter(order==1) %>% 
  #   select() %>% 
  #   left_join(ipd$by_serotype %>% dtrackr::untrack()) %>% 
  #   transmute(
  #     pneumo.phe_serotype=factor(pneumo.phe_serotype),
  #     pneumo.group = ifelse(is.na(pneumo.group),"PCV13-7",as.character(pneumo.group)),
  #     x= ifelse(is.na(x),0,x)
  #   )
  # 
  # serotype_prevalence = tmp_ipd %>% 
  #   # pick out only the relevant test results for the serotypes in the scenario from the count data.
  #   filter(pneumo.group %in% pcv_group$pneumo.group) %>%
  #   nest_join(pcv_group, by = "pneumo.group", name="pcv_group", unmatched = "drop") %>% 
  #   # normalise the distribution
  #   mutate(n= sum(x)) %>%
  #   mutate(binom::binom.confint(x,n,methods="wilson")) %>%
  #   transmute(
  #     test_id = pneumo.phe_serotype,
  #     pcv_id = pneumo.group,
  #     pcv_group = pcv_group,
  #     distribution = mean,
  #     false_neg_rate = 1-sens, # + rlnorm(n(), log(0.01), 0.5),
  #     n_diseased = n_diseased,
  #     false_pos_rate = 1-spec, # + rlnorm(n(), log(0.005), 0.5),
  #     n_controls = n_controls
  #   )
    


# Create a sample with exactly n*prev positives.
rfixed = function(boots, n, prev) {
  pos = round(n*prev)
  neg = n-pos
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

# group modify helper function
# g must contain grouping info only - group identifier and group_prevalence
# d is a data frame of with a test_id and a distribution column and a group
# result is a 
sim_dist = function(d, g, ..., n= 1000, boots = 1, exact=TRUE) {
  if (!"group_prevalence" %in% colnames(g)) stop("data must be grouped by at least the group_prevalence column")
  if (!"distribution" %in% colnames(d)) stop("the test data must define a distribution column")
  if (!"test_id" %in% colnames(d)) d = d %>% mutate(test_id = row_number())
  if (anyDuplicated(d$test_id)) stop("test_id must be unique (in each group)")
  n_tests = length(d$distribution)
  d = d %>% mutate(prevalence = g$group_prevalence*distribution/sum(distribution))
  
  # browser()
  # uniq = d %>% select(-distribution) %>% summarise(across(everything(), ~ n_distinct(.x,test_id)),.groups = "keep") %>% pivot_longer(cols = -test_id, names_to = "col", values_to = "uniqueness") %>% group_by(col) %>% filter(test_id == uniqueness) %>% pull(col) %>% unique()
  # pres = d %>% select(test_id, all_of(uniq))
  
  if (exact) {
    tmp = rfixed_mnom(boots, n, dist = d$distribution, g$group_prevalence)
  } else {
    td = c(1-g$group_prevalence, d$prevalence)
    tmp = lapply(1:boots, function(...) {
      (rmultinom(n, 1, td) %>% apply(MARGIN = 2, function(x) which(x==1)))-1
    }) %>% unlist()
  }
  tmp2 = tibble(
      boot = unlist(lapply(1:boots, rep, n*n_tests)),
      id = rep(rep(1:n, boots), n_tests),
      test_id = unlist(lapply(d$test_id, rep, n*boots)),
      actual = lapply(1:n_tests, function(x) tmp==x) %>% unlist() %>% as.integer()
    ) %>% 
    left_join(d %>% rename_with(~paste0("test_",.x), .cols=-test_id), by="test_id",suffix = c("",".old")) %>% 
    select(-ends_with(".old")) 
  return(tmp2 %>% group_by(across(starts_with("test"))))
}
# e.g.
# tmp = tibble(group_prevalence = 0.5, distribution = c(1,2,3,4), name = c("a","b","c","d"), note = c("x","x","y","y")) %>% group_by(group_prevalence) %>% group_modify(sim_dist, n=500, boots=1)
# tmp %>% group_by(boot,test_id) %>% summarise(pos = sum(actual))

# group_modify helper function
# Create a simulated population from a grouped dataframe which has a 
# is grouped by metadata including a prevalence column defining the 
# prevalence in the group.
# The result will have 1:boots replicates of 1:n samples containing
# id (of sample), boot (id of replicate) and actual (test positivity)
sim_pop = function(d, g, ..., n=1000, boots=100, exact=FALSE) {
  if (!"prevalence" %in% colnames(g)) stop("data must be grouped by at least the prevalence column")
  return(tibble(
    boot = unlist(lapply(1:boots, rep, n)),
    id = rep(1:n, boots),
    actual = if(!exact) {
      stats::rbinom(n=n*boots, size=1, p=g$prevalence)
    } else {
      rfixed(boots, n, g$prevalence)
    }
  ))
}



# group_modify helper function
# Modify a grouped dataframe to include a generate a test positivity column (and
# a sens & spec column). The grouped dataframe should be grouped by metadata
# (e.g. boot, id) including a test_id
#
# * false_pos_controls: the number of positives that appeared in the specificity
#   disease-free control group. These are by definition false positives. This
#   is (1-specificity)*n_controls
# * n_controls the number of controls in the specificity disease-free control group. 
# * false_neg_diseased the number of negatives that appeared in the sensitivity
#   confirmed disease group. These are by definition false negatives. This
#   is (1-sensitivity)*n_controls
# * n_diseased the number of confirmed disease cases in the sensitivity control group.
# 
# The ungrouped columns should contain an "actual" column containing the
# condition positive status (e.g. 1: diseased, 0: disease-free)
sim_test = function(d, g, ...) {
  
  g = g %>% mutate(
    sens_pos = n_diseased-false_neg_diseased,
    sens_neg = false_neg_diseased,
    control_pos = false_pos_controls,
    control_neg = n_controls-false_pos_controls
  )
  
  return(d %>% mutate(
    # unknown specificity defined by control test results. This creates the Beta-binomial
    test = 
      # actual positives * sensitivity gives 
      # actual * stats::rbinom(nrow(.), 1, stats::rbeta(nrow(.),shape1 = g$sens_pos, shape2=g$sens_neg)) + 
      # (1-actual) * stats::rbinom(nrow(.), 1, stats::rbeta(nrow(.),shape1 = g$control_pos, shape2=g$control_neg)),
      # 
      as.integer(
        actual * extraDistr::rbbinom(nrow(.), 1, alpha = g$sens_pos, beta=g$sens_neg) +
        (1-actual) * extraDistr::rbbinom(nrow(.), 1, alpha = g$control_pos, beta=g$control_neg)
      ),
    
    sens = g$sens_pos/(g$sens_neg+g$sens_pos), #TPR
    spec = g$control_neg/(g$control_neg+g$control_pos) #TNR
  ))
}

# group_modify helper function
# Modify a grouped dataframe to include a generate a test positivity column (and
# a sens & spec column). The grouped dataframe should be grouped by metadata
# (e.g. boot, id) including a test_id
#
# The ungrouped columns should contain an "actual" column containing the
# condition positive status (e.g. 1: diseased, 0: disease-free)
sim_test_2 = function(d, g, ..., exact = TRUE) {
  
  g = g %>% mutate(
    sens_pos = test_n_diseased-test_false_neg_diseased,
    sens_neg = test_false_neg_diseased,
    control_pos = test_false_pos_controls,
    control_neg = test_n_controls-test_false_pos_controls,
    test_sens = (test_n_diseased-test_false_neg_diseased) / test_n_diseased,
    test_spec = (test_n_controls-test_false_pos_controls) / test_n_controls
  )
  
  if (exact) {
    pos_count = sum(d$actual == 1)
    neg_count = sum(d$actual == 0)
    # browser()
    fp_tp = rfixed(1, pos_count, g$test_sens)
    fn_tn = rfixed(1, neg_count, 1-g$test_spec)
    pos_idx = cumsum(d$actual == 1)
    neg_idx = cumsum(d$actual == 0)
    test = ifelse(d$actual == 1, fp_tp[pos_idx], fn_tn[neg_idx])
    tmp = d %>% mutate(test = test)
  } else {
    tmp = d %>% mutate(
      test = as.integer(
          actual * extraDistr::rbbinom(nrow(.), 1, alpha = g$sens_pos, beta=g$sens_neg) +
            (1-actual) * extraDistr::rbbinom(nrow(.), 1, alpha = g$control_pos, beta=g$control_neg)
      )
    )
  }
    
  return(tmp %>% mutate(
    test_sens = g$test_sens, #TPR
    test_spec = g$test_spec #TNR
  ))
}

## IPD scenario ----

tmp = fs::path_expand(fs::path(rappdirs::user_cache_dir("testerror"),"vignettes",".cache"))
fs::dir_create(tmp)
cd = memoise::cache_filesystem(tmp)

# Pneumo serotype example ----
# This uses a serotype distribution identified in IPD paper
# for the baseline:
# Assume component sensitivity 0.8
# Assume component specificity 0.9975

# This function generates a realistic serotype distribution 
do_scenario_2 = function(
    pcv_group = c("PCV20","PCV15","PCV13","PCV7"),
    sens = 0.8, 
    spec = 0.9975, 
    # false_pos_controls = (n_controls)*(1-spec),
    n_controls = 2/(1-spec),
    # false_neg_diseased = (n_diseased)*(1-sens),
    n_diseased = 20/(1-sens),
    group_prevalence = seq(0.025, 0.2, 0.025),
    group = sprintf("panel prev: %1.3f", group_prevalence),
    samples = TRUE
) {
  group = enexpr(group)
  serotype_prevalence = ipd_distribution(pcv_group) %>%
    transmute(
      test_id = pneumo.phe_serotype,
      pcv_group = pcv_group,
      distribution = distribution
    ) %>%
    cross_join(
      tidyr::crossing(
        false_neg_rate = 1-sens, # + rlnorm(n(), log(0.01), 0.5),
        n_diseased = n_diseased,
        false_pos_rate = 1-spec, # + rlnorm(n(), log(0.005), 0.5),
        n_controls = n_controls,
        group_prevalence = group_prevalence,
      ) %>% mutate(
        group = !!group
      )
    ) %>%
    mutate( 
      false_neg_diseased = false_neg_rate * n_diseased, 
      false_pos_controls = false_pos_rate * n_controls, 
      # group = factor(sprintf("group %d",group))
    ) %>%
    group_by(group, group_prevalence, false_neg_rate, n_diseased, false_pos_rate, n_controls) %>%
    group_modify(function(d,g,...) {
      d %>% mutate(
        group_size = n(),
        prevalence = distribute(distribution, g$group_prevalence),
        expected_apparent_prevalence = prevalence * (1-g$false_neg_rate) + (1-prevalence) * g$false_pos_rate
      )
    })
  
  if (!samples) return(serotype_prevalence)
  
  # browser()
  set.seed(1001)
  serotype_tests = serotype_prevalence %>%
    group_by(across(c(starts_with("group")))) %>%
    # , false_neg_rate, n_diseased, false_pos_rate, n_controls
    group_modify(sim_dist, n=1000, boots=1, exact=TRUE) %>%
    group_by(across(c(starts_with("group"),starts_with("test")))) %>%
    group_modify(sim_test_2) 
  
  return(serotype_tests %>% group_by(group))
}

do_scenario = memoise::memoise(do_scenario_2,cache = cd)

## Panel functions ----

estimate_panel_performance = function(serotype_tests, test_id_col) {
  test_id_col = ensym(test_id_col)
  grps = serotype_tests %>% groups()
  # calculate the group sens and spec based on observed apparent prevalence
  tmp2 = serotype_tests %>% 
    group_by(boot, !!!grps, !!test_id_col, sens, spec) %>%
    summarise(
      component_ap = sum(test)/n()
    ) %>% 
    group_by(boot, !!!grps) %>% 
    summarise(
      panel_sens_est = testerror::panel_sens_estimator(component_ap, sens = sens, spec = spec)
    )
  # calculate the group sens and spec based on true panel prevalence
  tmp3 = serotype_tests %>%
    group_by(across(c(boot, !!!grps, !!test_id_col, sens, spec, prevalence, starts_with("test"), -test))) %>%#
    summarise(
      # combine single test counts
      actual = sum(actual),
      test = sum(test),
      total = n()
    ) %>%
    group_by(boot, !!!grps) %>%
    summarise(
      panel_prevalence = 1-prod(1-prevalence),
      panel_sens = testerror::panel_sens(prevalence, sens = sens, spec = spec),
      panel_spec = testerror::panel_spec(spec),
      .groups="drop"
    )
  
  tmp = serotype_tests %>%
    group_by(boot, !!!grps, id) %>%
    summarise(
      # combine tests per person
      actual = any(actual == 1),
      test = any(test == 1)
    ) 
  tmp = tmp %>%
    group_by(boot, !!!grps) %>%
    summarise(
      # summarise tests for all people
      # actual is the real ground truth / test is the test result including error
      total = n(),
      TP = sum(test == 1 & actual == 1),
      TN = sum(test == 0 & actual == 0),
      FP = sum(test == 1 & actual == 0),
      FN = sum(test == 0 & actual == 1),
      actual = sum(actual == 1),
      test = sum(test == 1),
      panel_apparent_prevalence = test/total,
      .groups="drop_last"
    ) 
  tmp = tmp %>%
    mutate(
      # TODO: confidence intervals?
      # sens_est = ifelse(TP+FN==0, NA, TP/(TP+FN)),
      # spec_est = ifelse(TN+FP==0, NA, TN/(TN+FP))
      binom_ci_2(TP, TP+FN, "sens_est"),
      binom_ci_2(TN, TN+FP, "spec_est")
    )
  tmp4 = tmp %>% inner_join(tmp2, by=join_by(boot, !!!grps)) %>% inner_join(tmp3, by=join_by(boot, !!!grps))
  tmp4 = tmp4 %>% mutate(
    rogan_gladen = rogan_gladen(panel_apparent_prevalence, sens =  panel_sens,spec = panel_spec)
  )
  return(tmp4)  
}


estimate_panel_performance_uncertain = function(serotype_tests, test_id_col) {
  test_id_col = ensym(test_id_col)
  grps = serotype_tests %>% groups()
  # calculate the group sens and spec based on observed apparent prevalence
  tmp2 = serotype_tests %>% 
    group_by(boot, !!!grps, !!test_id_col, prior_sens, prior_spec) %>%
    summarise(
      component_ap = sum(test)/n()
    ) %>% 
    group_by(boot, !!!grps) %>% 
    summarise(
      panel_sens_est = testerror::panel_sens_estimator(component_ap, sens = prior_sens, spec = prior_spec)
    )
  # calculate the group sens and spec based on true panel prevalence
  tmp3 = serotype_tests %>%
    group_by(across(c(boot, !!!grps, !!test_id_col, starts_with("prior"), starts_with("test"), -test))) %>%#
    summarise(
      # combine single test counts
      actual = sum(actual),
      test = sum(test),
      total = n()
    ) %>%
    group_by(boot, !!!grps) %>%
    summarise(
      panel_prevalence = 1-prod(1-actual/total),
      panel_sens = testerror::panel_sens(actual/total, sens = prior_sens, spec = prior_spec),
      panel_spec = testerror::panel_spec(prior_spec),
      panel_sens_samples = list(testerror::uncertain_panel_sens_estimator(
        pos_obs = test, n_obs = total,
        false_pos_controls = prior_false_pos_controls, n_controls = prior_n_controls,
        false_neg_diseased = prior_false_neg_diseased, n_diseased = prior_n_diseased
      )),
      panel_spec_samples = list(testerror::uncertain_panel_spec(
        false_pos_controls = prior_false_pos_controls, n_controls = prior_n_controls
      ))
    ) %>%
    group_by(across(c(starts_with("group"),starts_with("panel")))) %>%
    mutate(
      panel_sens_beta = purrr::map(panel_sens_samples, ~ as_tibble(beta_fit(.x))),
      panel_spec_beta = purrr::map(panel_spec_samples, ~ as_tibble(beta_fit(.x))),
    ) %>%
    unnest(cols = c(panel_sens_beta, panel_spec_beta),names_sep = ".")
  
  tmp = serotype_tests %>%
    group_by(boot, !!!grps, id) %>%
    summarise(
      # combine tests per person
      actual = any(actual == 1),
      test = any(test == 1)
    ) 
  tmp = tmp %>%
    group_by(boot, !!!grps) %>%
    summarise(
      # summarise tests for all people
      # actual is the real ground truth / test is the test result including error
      total = n(),
      TP = sum(test == 1 & actual == 1),
      TN = sum(test == 0 & actual == 0),
      FP = sum(test == 1 & actual == 0),
      FN = sum(test == 0 & actual == 1),
      actual = sum(actual == 1),
      test = sum(test == 1),
      panel_apparent_prevalence = test/total,
      .groups="drop_last"
    ) 
  tmp = tmp %>%
    mutate(
      # TODO: confidence intervals?
      # sens_est = ifelse(TP+FN==0, NA, TP/(TP+FN)),
      # spec_est = ifelse(TN+FP==0, NA, TN/(TN+FP))
      binom_ci_2(TP, TP+FN, "sens_est"),
      binom_ci_2(TN, TN+FP, "spec_est")
    )
  tmp4 = tmp %>% inner_join(tmp2, by=join_by(boot, !!!grps)) %>% inner_join(tmp3, by=join_by(boot, !!!grps))
  tmp4 = tmp4 %>% 
    group_by(boot,!!!grps) %>%
    mutate(
    rogan_gladen = rogan_gladen(panel_apparent_prevalence, sens =  panel_sens,spec = panel_spec),
    testerror::true_prevalence(
      pos_obs = test, n_obs = total,
      false_pos_controls = panel_spec_beta.shape2,
      n_controls = panel_spec_beta.conc,
      false_neg_diseased = panel_sens_beta.shape2,
      n_diseased = panel_sens_beta.conc
    ),
    testerror::uncertain_rogan_gladen(
      pos_obs = test, n_obs = total,
      sens = panel_sens_samples,
      spec = panel_spec_samples,
      prefix = "rogan_gladen"
    )
  )
  return(tmp4)  
}    

# dataframe function
# takes a dataframe containing boot and id columns and grouped by other metadata
# and combines all the ungrouped columns into a single summary test where a 
# positive result is the result of the logical OR of the tests, calculating
# ground truth error
combine_as_panel = function(serotype_tests) {
  grps = serotype_tests %>% groups()
  serotype_tests %>% 
    group_by(!!!grps, boot, id) %>%
    summarise(
      # actual is the real ground truth / test is the test result including error
      actual = any(actual == 1),
      test = any(test == 1),
      expected = 1-prod(1-prevalence),
      .groups="drop"
    ) %>%
    # this is only grouped by boot at this stage as we are mutating within groups to get
    # p.actual and p.test on a boot by boot basis.
    group_by(boot) %>%
    mutate(
      total = 1,
      TP = test == actual & actual == 1,
      TN = test == actual & actual == 0,
      FP = test != actual & actual == 0,
      FN = test != actual & actual == 1,
      p.actual = actual/sum(actual),
      p.test = test/sum(test)
    ) %>%
    group_by(!!!grps, boot) %>%
    summarise(
      across(-c(id), sum), 
      .groups="drop_last") %>%
    mutate(
      apparent_prevalence = test/total*100,
      prevalence_error = (test-actual)/total*100,
      prevalence_rel_error = (test/actual-1)*100,
      proportional_error = (p.test-p.actual)*100,
      sens_est = ifelse(TP+FN==0, NA, TP/(TP+FN)),
      spec_est = ifelse(TN+FP==0, NA, TN/(TN+FP)),
      ppv = ifelse(TP+FP==0, NA,  TP/(TP+FP)),
      acc = (TP + TN) / (TP+FP+FN+TN)
    ) %>%
    group_by(!!!grps)
} 

# pretty print a summary mean +/- 95% quantiles
mean_sd = function(x, ...) {
  sprintf("%1.3g [%1.3g \u2014 %1.3g]", mean(x, na.rm=TRUE), stats::quantile(x, 0.025, na.rm=TRUE), stats::quantile(x, 0.975, na.rm=TRUE))
}

# binom_ci(c(1,0,3),c(2,0,7))
binom_ci = function(x, n) {
  tmp = binom::binom.confint(x,n,methods = "wilson")
  ifelse(n==0, "\u2014", sprintf("%1.2f [%1.2f\u2013%1.2f]", tmp$mean, tmp$lower, tmp$upper))
}

# binom_ci_2(c(1,0,3),c(2,0,7))
binom_ci_2 = function(x, n, name="proportion") {
  tmp = binom::binom.confint(x,n,methods = "wilson")
  return(tibble::tibble(
    !!(sprintf("%s.0.5",name)) := ifelse(is.finite(tmp$mean), tmp$mean, NA_real_),
    !!(sprintf("%s.0.025",name)) := ifelse(is.finite(tmp$lower), tmp$lower, NA_real_),
    !!(sprintf("%s.0.975",name)) := ifelse(is.finite(tmp$upper), tmp$upper, NA_real_)
  ))
}

# Summarise multiple bootstrapped simulation results
# from 
# TODO: 
summarise_boots = function(panel) {
  grps = panel %>% groups()
  panel %>% 
    group_by(!!!grps, total, expected) %>%
    summarise(across(-boot, ~ mean_sd(.x, na.rm=TRUE)), .groups="drop") %>%
    mutate(true_prevalence = expected/total*100)
}

# For a panel of tests with a relative proportion of prevalence (rel_prop) solve the
# combination $(1-\prod_n(1-\frac{rel_prop_n)}{\sum_n{rel_prop_n}})$ to get a target 
# global prevalence (or just $prev*\frac{rel_prop_n)}{\sum_n{rel_prop_n}}$ the 
# if not global). This expects one row per test but can be used in a group modify
update_prevalence = function(d, g, ..., prevalence, global=TRUE) {
  if (!global) {
    d %>% mutate(
      prevalence = (rel_prop / sum(rel_prop))*prevalence
    )
  } else {
    d %>% mutate(
      prevalence = distibute(rel_prop, prevalence)
    )
  }
}

# abs(1-prod(1-distribute(c(1,2,3),0.25)) - 0.25) < 0.00001


## Plotting ----

# Create an apparent prevalence plot with VE
apparent_prevalence_plot = function(
    p = c(0.05,0.2)*0.5,
    ap = NULL,
    sens = 0.8,
    spec = 0.95,
    cols = c("blue","red"),
    VE = 1-p[1]/p[2],
    aVE = 1-ap(p[1])/ap(p[2]),
    lim = c(0,1),
    bottom_right = "sens: {sprintf('%1.4g',sens)}\nspec: {sprintf('%1.4g',spec)}",
    top_left = "apparent VE: {sprintf('%1.2g%%',aVE*100)}\nspec: {sprintf('%1.2g%%',VE*100)}",
    annotate_cols = list(
      bottom_right = "black",
      top_left = "magenta"
    )
) {
  
  if (!is.null(ap)) p = rogan_gladen(ap, sens, spec)
  ap = function(p) apparent_prevalence(p,sens, spec)
  
  annotate = list(
    bottom_right = paste0(glue::glue(bottom_right), collapse="\n"),
    top_left = paste0(glue::glue(top_left), collapse="\n")
  )
  
  # point ap = p
  # p = p*sens+(1-p)*(1-spec)
  # p = p*(sens+spec-1)+(1-spec)
  px = (1-spec) / (1-sens+1-spec)
  pzero =  rep(lim[[1]], length(p))
  
  segments = tibble(
    x = c(lim[1], p, p ), 
    y = c(ap(lim[1]), pzero, ap(p) ), 
    xend = c(lim[2], p, pzero ), 
    yend = c(ap(lim[2]), ap(p), ap(p)), 
    col = c("black", cols, cols),
    linetype = c("solid", rep("dashed",length(p)*2))
  )
  
  n = 0.95*lim[1]+0.05*lim[2]
  f = 0.05*lim[1]+0.95*lim[2]
  
  ann_pos = list(
    bottom_right = c(f,n,1,0),
    top_left = c(n,f,0,1)
  )
  
  labels = tibble(
    label = c("(1-spec)","sens",sprintf("(%1.2g,%1.2g)",px,px),sprintf("%1.2g",p), sprintf("%1.2g",ap(p))),
    x = c(lim, px, p, pzero),
    y = c(ap(lim), ap(px), pzero, ap(p)),
    hjust=c(1.1, -0.1, -0.1, rep(1.1,length(p)*2)),
    vjust=c(1, 0, 1.1, rep(0.5,length(p)*2)),
    angle=c(0,0,0, rep(90,length(p)), rep(0,length(p))),
    col =c("black","black","black",cols, cols),
    size = c(6,6,6,rep(6,length(p)*2)),
  ) %>% bind_rows(
    tibble(
      label = unlist(annotate),
      x = sapply(names(annotate), function(x) ann_pos[[x]][1]),
      y = sapply(names(annotate), function(x) ann_pos[[x]][2]),
      hjust = sapply(names(annotate), function(x) ann_pos[[x]][3]),
      vjust = sapply(names(annotate), function(x) ann_pos[[x]][4]),
      angle = 0,
      col = unlist(annotate_cols[names(annotate)]),
      size = 6
    )
  )
  
  points = tibble(
    x = c(px, p),
    y = c(ap(px), ap(p)),
    col = c("black", cols)
  )
  
  ggplot()+
    coord_fixed(xlim=lim, ylim=lim, clip="off")+
    xlab("True prevalence")+
    ylab("E(Apparent prevalence)")+
    geom_abline(colour="grey80")+
    geom_segment(aes(x=x,y=y,xend=xend,yend=yend,col=col,linetype=linetype), segments)+
    geom_text(mapping=aes(x=x,y=y,label=label,vjust=vjust,hjust=hjust,angle=angle, colour=col, size = size/ggplot2:::.pt), labels %>% filter(angle==0), direction="y")+
    geom_text(mapping=aes(x=x,y=y,label=label,vjust=vjust,hjust=hjust,angle=angle, colour=col, size = size/ggplot2:::.pt), labels %>% filter(angle==90), direction="x")+
    geom_point(aes(x=x,y=y,colour = col),points, size=0.5)+
    scale_x_continuous(breaks = lim,expand = c(0, 0))+
    scale_y_continuous(breaks = lim,expand = c(0, 0))+
    scale_colour_identity( aesthetics = c("color","linetype"),guide = "none")+
    scale_size_identity()+
    theme(
      plot.margin = unit(c(1,1,1,1), "lines"),
      #axis.title.x = element_text(margin=margin(t=2, unit="lines")),
      #axis.title.y = element_text(margin=margin(r=2, unit="lines"))
      axis.title.x = element_text(hjust = 1),
      axis.title.y = element_text(hjust = 1)
    ) 
  
}

no_y = function() {theme(axis.text.y = element_blank(), axis.title.y = element_blank())}
no_x = function() {theme(axis.text.x.bottom = element_blank(), axis.title.x = element_blank())}
horiz = function() {theme(legend.position = "bottom", legend.direction = "horizontal", 
                                       legend.box = "vertical", legend.justification = "center")}

std_size = list(
  A4 = list(width=8.25,height=11.75,rot=0),
  A5 = list(width=5+7/8,height=8.25,rot=0),
  full =  list(width=5.9,height=8,rot=0),
  landscape =  list(width=9.75,height=5.9,rot=0),
  half =  list(width=5.9,height=4,rot=0),
  third =  list(width=5.9,height=3,rot=0),
  two_third = list(width=5.9,height=6,rot=0),
  quarter = list(width=5.9,height=2,rot=0),
  quarter_portrait = list(width=3,height=4,rot=0),
  sixth = list(width=3,height=3,rot=0),
  slide = list(width=12,height=6,rot=0)
)

save_as = function(p, file, size = std_size$third, width = size$width, height=size$height) {
  tmp = fs::path_ext_remove(file)
  tmp_pdf = fs::path_ext_set(tmp,"pdf")
  fs::dir_create(fs::path_dir(tmp_pdf))
  tmp_png = fs::path_ext_set(tmp,"png")
  ggplot2::ggsave(plot = p, filename = tmp_pdf, width = width, height=height, device = cairo_pdf)
  suppressWarnings(pdftools::pdf_convert(tmp_pdf, dpi = 300,filenames = tmp_png,verbose = FALSE))
  if(.is_running_in_console()) {
    rstudioapi::viewer(tmp_png)
  }
  invisible(knitr::include_graphics(path = tmp_png, auto_pdf = TRUE, dpi=300))
}

.is_running_in_console = function() {
  isTRUE(try(
    rstudioapi::getActiveDocumentContext()$id == "#console" |
      !rstudioapi::getActiveDocumentContext()$path %>% stringr::str_ends("Rmd")
  ))
}
