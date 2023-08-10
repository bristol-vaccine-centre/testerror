// Use this when you have a set of samples for a given test and
// you have an estimate of the test sensitivity and specificity as a beta
// prior.
data {
  int<lower = 0> k_sample;
  int<lower = 0, upper=1> y_sample[k_sample];
  real<lower = 0> fp_spec;
  real<lower = 0> tn_spec;
  real<lower = 0> fn_sens;
  real<lower = 0> tp_sens;
}
parameters {
  real<lower=0, upper = 1> p;
  real<lower=0, upper = 1> spec;
  real<lower=0, upper = 1> sens;
}
model {
  real p_sample = p * sens + (1 - p) * (1 - spec);
  y_sample ~ bernoulli(p_sample);
  spec ~ beta(tn_spec, fp_spec);
  sens ~ beta(tp_sens, fn_sens);
}

