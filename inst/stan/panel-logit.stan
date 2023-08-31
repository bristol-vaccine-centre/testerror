// Use this when we have control group data plus an estimate of panel
// sensitivity and specificity as a beta prior.
data {
  int<lower = 1> n_test;
  
  // Counts of samples: 
  // positives per test (pos_sample) and per panel (pos_sample_combined)
  // denominators k_sample and k_sample_combined likely to be the same
  // these quantities are calculated from raw linelist data.
  int<lower = 0> k_sample[n_test];
  int<lower = 0> pos_sample[n_test];
  int<lower = 0> k_sample_combined;
  int<lower = 0> pos_sample_combined;
  
  // component priors
  // sensitivity ans specificity.
  // These are going to be determined by expert option but sensible values might
  // be given by `ci_to_logitnorm(median=0.9975, lower=0.9, upper=0.99999)` and
  // `ci_to_logitnorm(median=0.8, lower=0.6, upper=0.95)` in general the 
  // heuristic mu = logit(0.9975) and sigma = abs(mu)/2 is pretty good. If these
  // are unknown then uninformative prior is any mu value with a sigma>50
  real logit_mu_spec_prior[n_test];
  real logit_mu_sens_prior[n_test];
  real<lower = 0> logit_sigma_spec_prior[n_test];
  real<lower = 0> logit_sigma_sens_prior[n_test];
  
  // component prevalence prior
  // this should be a very weak prior with just enough to 
  // give the sampler something to work with. A good choice
  // would be the logit of lang-reiczigel transform of the observed test positives, but
  // the bounds cannot be equal zero. a sigma less than 5 is likely to be too constrained
  // values from `ci_to_logitnorm(median=??, lower=0.000000001, upper=0.5)` gives a
  // reasonable starting point. another good heuristic is mu = logit(mean) and
  // sigma = abs(mu)/2
  real logit_mu_p_prior[n_test];
  real<lower = 0> logit_sigma_p_prior[n_test];
  
  // specificity negative control data summary
  int<lower = 0> fp_spec[n_test];
  int<lower = 0> k_spec[n_test];
  
  // sensitivity postivie control data summary
  int<lower = 0> fn_sens[n_test];
  int<lower = 0> k_sens[n_test];
  
  // panel priors
  // the same guidance for component priors holds here. sigmas probably 
  real logit_mu_panel_spec_prior;
  real logit_mu_panel_sens_prior;
  real<lower = 0> logit_sigma_panel_spec_prior;
  real<lower = 0> logit_sigma_panel_sens_prior;
  
}
parameters {
  // The real prevalence of each of the components of the test
  // The offset and multiplier versions require stan 2.32 for vectors
  // vector<offset = logit_mu_p_prior, multiplier = logit_sigma_p_prior>[n_test] logit_p;
  // vector<offset = logit_mu_spec_prior, multiplier = logit_sigma_spec_prior>[n_test] logit_spec;
  // vector<offset = logit_mu_sens_prior, multiplier = logit_sigma_sens_prior>[n_test] logit_sens;
  vector[n_test] logit_p;
  vector[n_test] logit_spec;
  vector[n_test] logit_sens;
  
  
  // The real prevalence of the combination of the panel
  real logit_panel_p;
  real<offset = logit_mu_panel_spec_prior, multiplier = logit_sigma_panel_spec_prior> logit_panel_spec;
  real<offset = logit_mu_panel_sens_prior, multiplier = logit_sigma_panel_sens_prior> logit_panel_sens;

}
transformed parameters {  
  // The specificity of the component test
  vector<lower=0, upper = 1>[n_test] spec = inv_logit(logit_spec);
  // The sensitivity of the component test
  vector<lower=0, upper = 1>[n_test] sens = inv_logit(logit_sens);
  // The sensitivity of the component test
  vector<lower=0, upper = 1>[n_test] p = inv_logit(logit_p);
  
  // The real prevalence of the combination of the panel 
  // this assumes independence (which as also assumed elsewhere)
  // real p_combined = 1-prod(1-p);
  // this is the panel specificity from the formal analysis (supplementary 1)
  real<lower=0, upper = 1> panel_spec = prod(spec);
  // this is the panel sensitivity from the formal analysis (supplementary 1)
  real<lower=0, upper = 1> panel_sens = 1-(
    prod((1-sens) .* p + spec .* (1-p) ) - prod(spec .* (1-p))
  )/(
    1 - prod(1-p)
  );
  // this is the theoretical panel prevalence from the formal analysis
  // we use this as a prior
  real<lower=0, upper = 1> panel_p = 1-prod(1-p);
  
  // The specificity of the panel test
  real<lower=0, upper = 1> spec_combined = inv_logit(logit_panel_spec);
  // The sensitivity of the panel test
  real<lower=0, upper = 1> sens_combined = inv_logit(logit_panel_sens);
  // The prevalence of the combined condition
  real<lower=0, upper = 1> p_combined = inv_logit(logit_panel_p);
  
}
model {
  
  // the expected value of the apparent prevalence of the panel test:
  real p_sample_combined = p_combined * sens_combined + (1 - p_combined) * (1 - spec_combined);
  
  // There are 2 sources of information about the sensitivity and specificity
  // of the panel. It may be directly supplied as a hyperparameter as result of a panel level
  // experiment or may be derived from the component test control data if present.
  // 1) panel_wide sensitivity and specificity priors.
  target += normal_lpdf(logit_panel_spec | logit_mu_panel_spec_prior, logit_sigma_panel_spec_prior);
  target += normal_lpdf(logit_panel_sens | logit_mu_panel_sens_prior, logit_sigma_panel_sens_prior);
  
  // 2) panel wide sensitivity specificity and prevalence from theoretical values
  // calculated from components. 
  // Theoretical sensitivity from components is a strong prior. The sigma value of 
  // 0.05 allowing about a +/-10% variation from theoretical 95% specificity.
  // however we may have very little information about the panel sensitivity
  // in which case the sampling will be fairly variable.
  target += normal_lpdf(logit_panel_spec | logit(panel_spec), 0.05);
  target += normal_lpdf(logit_panel_sens | logit(panel_sens), 0.05);
  
  // Panel prevalence is the true unknown.
  // Again 2 sources of information. The theoretical combined prevalence of 
  // the components (assuming independence), and the observed test positivity from
  // components which may be lower or higher if independence assumption is not
  // true
  
  // 1) Theoretical prevalence from components is a reasonably strong prior. The sigma value of 
  // 0.2 allowing about a +/-50% variation at 10% prevalence. 
  target += normal_lpdf(logit_panel_p | logit(panel_p), 0.2);
  
  // 2) observed test positives of combination panel.
  // weight this part of the model as more important than the components.
  target += n_test * binomial_lpmf(pos_sample_combined | k_sample_combined, p_sample_combined);
  
  
  // Moving onto the components
  for (n in 1:n_test) {
    
    // test-wise expected value of apparent prevalence of each component
    real p_sample = p[n] * sens[n] + (1 - p[n]) * (1 - spec[n]);
    
    // panel_wide sensitivity and specificity priors of each component.
    target += normal_lpdf(logit_spec[n] | logit_mu_spec_prior[n], logit_sigma_spec_prior[n]);
    target += normal_lpdf(logit_sens[n] | logit_mu_sens_prior[n], logit_sigma_sens_prior[n]);
  
    // test-wise sensitivity and specificity are also informed by control group
    // data for each component.
    if (k_spec[n]>0) {
      target += binomial_lpmf(fp_spec[n] | k_spec[n], 1-spec[n]);
    }
    if (k_sens[n]>0) {
      target += binomial_lpmf(fn_sens[n] | k_sens[n], 1-sens[n]);
    }
    
    // test-wise prior prevalence assumptions for each component
    target += normal_lpdf(logit_p[n] | logit_mu_p_prior[n], logit_sigma_p_prior[n]);
    // test-wise observed test results for each component
    target += binomial_lpmf(pos_sample[n] | k_sample[n], p_sample);
    
  }
}
