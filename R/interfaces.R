.input_data = interfacer::iface(
  id = character ~ "the patient identifier",
  test = factor ~ "the test type",
  result = logical ~ "the test result",
  .groups = FALSE
)

.input_panel_data = interfacer::iface(
  id = character ~ "the patient identifier",
  result = logical ~ "the panel result",
  .groups = FALSE
)

.output_data = interfacer::iface(
  test = character ~ "the name of the test or panel",
  prevalence.lower = numeric ~ "the lower estimate",
  prevalence.median = numeric ~ "the median estimate",
  prevalence.upper = numeric ~ "the upper estimate",
  prevalence.method = character ~ "the method of estimation",
  prevalence.label = character ~ "a fomatted label of the true prevalence estimate with CI",
  .groups = FALSE
)