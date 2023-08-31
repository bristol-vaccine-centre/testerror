// Use this when we have control group data plus an estimate of panel
// sensitivity and specificity as a beta prior.
data {
  
  int<lower = 1> n;
  
  // Counts of samples: 
  // positives per test (pos_sample) and per panel (pos_sample_combined)
  // denominators k_sample and k_sample_combined likely to be the same
  // these quantities are calculated from raw linelist data.
  int<lower = 0> k_sample[n];
  int<lower = 0> pos_sample[n];
  
  // component priors
  // sensitivity ans specificity.
  // These are going to be determined by expert option but sensible values might
  // be given by `ci_to_logitnorm(median=0.9975, lower=0.9, upper=0.99999)` and
  // `ci_to_logitnorm(median=0.8, lower=0.6, upper=0.95)` in general the 
  // heuristic mu = logit(0.9975) and sigma = abs(mu)/2 is pretty good. If these
  // are unknown then uninformative prior is any mu value with a sigma>50
  real logit_mu_spec_prior[n];
  real logit_mu_sens_prior[n];
  real<lower = 0> logit_sigma_spec_prior[n];
  real<lower = 0> logit_sigma_sens_prior[n];
  
  // component prevalence prior
  // this should be a very weak prior with just enough to 
  // give the sampler something to work with. A good choice
  // would be the logit of lang-reiczigel transform of the observed test positives, but
  // the bounds cannot be equal zero. a sigma less than 5 is likely to be too constrained
  // values from `ci_to_logitnorm(median=??, lower=0.000000001, upper=0.5)` gives a
  // reasonable starting point. another good heuristic is mu = logit(mean) and
  // sigma = abs(mu)/2
  real logit_mu_p_prior[n];
  real<lower = 0> logit_sigma_p_prior[n];
  
  // specificity negative control data summary
  int<lower = 0> fp_spec[n];
  int<lower = 0> k_spec[n];
  
  // sensitivity positive control data summary
  int<lower = 0> fn_sens[n];
  int<lower = 0> k_sens[n];
  
}
parameters {
  // The real prevalence of each of the components of the test
  // The offset and multiplier versions require stan 2.32 for vectors
  // vector<offset = logit_mu_p_prior, multiplier = logit_sigma_p_prior>[n_test] logit_p;
  // vector<offset = logit_mu_spec_prior, multiplier = logit_sigma_spec_prior>[n_test] logit_spec;
  // vector<offset = logit_mu_sens_prior, multiplier = logit_sigma_sens_prior>[n_test] logit_sens;
  vector[n] logit_p;
  vector[n] logit_spec;
  vector[n] logit_sens;

}
transformed parameters {  
  // The specificity of the component test
  vector<lower=0, upper = 1>[n] spec = inv_logit(logit_spec);
  // The sensitivity of the component test
  vector<lower=0, upper = 1>[n] sens = inv_logit(logit_sens);
  // The sensitivity of the component test
  vector<lower=0, upper = 1>[n] p = inv_logit(logit_p);
  
}
model {
  
  // the expected value of the apparent prevalence of the panel test:
  vector[n] p_sample = p .* sens + (1 - p) .* (1 - spec);
    
  // N.b. all vectorised
  
  // sensitivity and specificity priors of each component.
  target += normal_lpdf(logit_spec | logit_mu_spec_prior, logit_sigma_spec_prior);
  target += normal_lpdf(logit_sens | logit_mu_sens_prior, logit_sigma_sens_prior);

  // test-wise sensitivity and specificity are also informed by control group
  // data for each component.
  target += binomial_lpmf(fp_spec | k_spec, 1-spec);
  target += binomial_lpmf(fn_sens | k_sens, 1-sens);
  
  // test-wise prior prevalence assumptions for each component
  target += normal_lpdf(logit_p | logit_mu_p_prior, logit_sigma_p_prior);
  // test-wise observed test results for each component
  target += binomial_lpmf(pos_sample | k_sample, p_sample);

}
