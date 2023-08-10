// Use this for single test adjustment when we an estimate of sensitivity and specificity as a beta prior and/or data.
// this is vectorised...
data {
  
  int<lower=1> n;
  
  // Counts of samples: 
  // positives per test (pos_sample) and per panel (pos_sample_combined)
  // denominators k_sample and k_sample_combined likely to be the same
  int<lower = 1> k_sample[n];
  int<lower = 0> pos_sample[n];
  
  // component specificity priors based on expert opinion or previous studies
  real<lower = 0> fp_spec_prior[n];
  real<lower = 0> tn_spec_prior[n];
  
  // component specificity priors based on expert opinion or previous studies
  real<lower = 0> fn_sens_prior[n];
  real<lower = 0> tp_sens_prior[n];
  
  // specificity negative control data summary
  int<lower = 0> fp_spec[n];
  int<lower = 0> k_spec[n];
  
  // sensitivity postivie control data summary
  int<lower = 0> fn_sens[n];
  int<lower = 0> k_sens[n];

}
parameters {
  // The real prevalence of each of the components of the test
  vector<lower=0, upper = 1>[n] p;
  // The specificity of the component test
  vector<lower=0, upper = 1>[n] spec;
  // The sensitivity of the component test
  vector<lower=0, upper = 1>[n] sens;
  
}
model {
  
  // the expected value of the apparent prevalence of the panel test:
  vector[n] p_sample = p .* sens + (1 - p) .* (1 - spec);
  
  // sensitivity and specificity priors.
  // set to 1,1 if nothing already known
  target += beta_lpdf(spec | tn_spec_prior, fp_spec_prior);
  target += beta_lpdf(sens | tp_sens_prior, fn_sens_prior);
  
  // group observed positives test results of combined panel.
  // pos_sample_combined[n] ~ binomial(p_sample_combined);
  // weight this part of the model as more important than the components.
  target += binomial_lpmf(pos_sample | k_sample, p_sample);
  
  // The data for sensitivity and specificity
  // this should just work if k_spec is zero
  target += binomial_lpmf(fp_spec | k_spec, 1-spec);
  target += binomial_lpmf(fn_sens | k_sens, 1-sens);
  
}
