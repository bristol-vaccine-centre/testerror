// Original gelman model with renaming
// Use this when you have counts of positives in a test group as well as in 
// a disease negative control group and a disease positive control group
data {
  int<lower = 0> pos_sample;
  int<lower = 0> n_sample;
  int<lower = 0> fp_spec;
  int<lower = 0> n_spec;
  int<lower = 0> fn_sens;
  int<lower = 0> n_sens;
}
parameters {
  real<lower=0, upper = 1> p;
  real<lower=0, upper = 1> spec;
  real<lower=0, upper = 1> sens;
}
model {
  real p_sample = p * sens + (1 - p) * (1 - spec);
  pos_sample ~ binomial(n_sample, p_sample);
  fp_spec ~ binomial(n_spec, 1-spec);
  fn_sens ~ binomial(n_sens, 1-sens);
}

