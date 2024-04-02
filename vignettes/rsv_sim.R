rsv_prev = 0.05
sim_size = 100000

tests = tibble::tribble(
  ~test_id, ~ sens, ~ spec, ~p_performed,
  "lab", 0.5, 0.99, 0.4,
  "poc", 0.8, 0.90, 0.9
  # you can add more tests here
  
  # N.B. p_performed is probability test was done on an individual (assumed random at the moment)
  # possible to split this into a p_performed_given_disease and a p_performed_given_no_disease
  # depending on scenario (e.g. maybe 
  #   lab: p_performed_given_disease == p_performed_given_no_disease,
  #   poc: p_performed_given_disease > p_performed_given_no_disease, 
  # to account for clinical suspicion etc.)
)

# Generation simulated data ----

# data set in long format with one test result per row.
# true_rsv: 0=no rsv; 1=rsv
# test_result: 1=pos; 0=neg; NA=not done
synth_data = tibble::tibble(
  pat_id = 1:sim_size,
  # generate a disease positive flag based on prevalence 
  # this gives you a random sample close to your prevalence
  #   true_rsv = rbinom(sim_size,1,rsv_prev)
  # this gives you exact prevalence in your simulation.
  true_rsv = sample(c(rep(1,sim_size*rsv_prev), rep(0,sim_size* (1-rsv_prev))))
) %>% 
  dplyr::cross_join(tests) %>% 
  dplyr::group_by(test_id) %>%
  dplyr::mutate(
    test_result = ifelse(
      rbinom(sim_size,1,p_performed) == 0, 
        # set test results which were not performed to NA
        NA,
        # otherwise get a result that is consistent with test sens and spec
        true_rsv*rbinom(sim_size,1,sens) + (1-true_rsv)*(1-rbinom(sim_size,1,spec))
        # if we were thinking about test delay we would have to do that here.
    )
  )

# data set with a test per column (1-pos, 0-neg, NA-not done)
wide_synth_data = synth_data %>%
  dplyr::select(pat_id, true_rsv, test_id, test_result) %>%
  tidyr::pivot_wider(names_from = test_id, values_from = test_result)

wide_synth_data %>% dplyr::glimpse()

# Check simulation is doing the right thing. ----

rogan_gladen = function(ap, sens, spec) {
  dplyr::case_when(
    ap <= 1-spec ~ 0,
    sens <= ap ~ 1,
    TRUE ~ (ap + spec -1)/(sens+spec-1)
  )
}

# can reconstruct the input parameters from the data?
reconstruct = synth_data %>% 
  # the input parameters
  dplyr::group_by(test_id, sens, spec, p_performed) %>% 
  # equivalents from the data
  dplyr::summarise(
    rsv_prev = mean(true_rsv),
    test_positivity = mean(test_result,na.rm = TRUE), 
    test_coverage = mean(!is.na(test_result))
  ) %>% dplyr::mutate(
    rsv_prev_est = rogan_gladen(test_positivity, sens, spec)
  )

# not too bad.
reconstruct %>% dplyr::glimpse()
