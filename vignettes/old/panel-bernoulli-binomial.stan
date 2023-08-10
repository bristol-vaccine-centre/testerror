// Use this when we have control group data but no estimates of panel
// sensitivity and specificity (e.g. it is an adhoc panel that hasn;t
// been characterised)
data {
  int<lower = 1> n_test;
  
  // samples are rows; tests are columns although this is transposed in R
  // a positive result is 1 a negative result is 0
  int<lower = 1> k_sample;
  int<lower = 0,upper = 1> y_sample[k_sample,n_test];
  
  // specificity controls data summary
  int<lower = 0> fp_spec[n_test];
  int<lower = 0> k_spec[n_test];
  
  // sensitivity controls data summary
  int<lower = 0> fn_sens[n_test];
  int<lower = 0> k_sens[n_test];
}
transformed data {
  
  // combine the results of n_test tests, for each of n_sample samples 
  // (1=pos, 0=neg) with logical OR.
  // this has to be done by iteration because of the way stan handles integers.
  int<lower = 0,upper = 1> y_sample_combined[k_sample];
  for (k in 1:k_sample) {
    int tmp = 1;
    for (n in 1:n_test) {
      tmp *= 1-y_sample[k,n];
    }
    y_sample_combined[k] = 1-tmp;
  }
}
parameters {
  // The real prevalence of each of the components of the test
  vector<lower=0, upper = 1>[n_test] p;
  // The specificity of the component test
  vector<lower=0, upper = 1>[n_test] spec;
  // The sensitivity of the component test
  vector<lower=0, upper = 1>[n_test] sens;
  // The real prevalence of the combination of the panel
  real<lower=0, upper = 1> p_combined;
}
transformed parameters {
  
  // The real prevalence of the combination of the panel 
  // this assumes independence (which as also assumed elsewhere)
  // real p_combined = 1-prod(1-p);
  // this is the panel specificity from the formal analysis (supplementary 1)
  real spec_combined = prod(spec);
  // this is the panel sensitivity from the formal analysis (supplementary 1)
  real sens_combined = 1-(
    prod((1-sens) .* p + spec .* (1-p) ) - prod(spec .* (1-p))
  )/(
    1 - prod(1-p)
  );
}
model {
  
  // the expected value of the apparent prevalence of the panel test:
  real p_sample_combined = p_combined * sens_combined + (1 - p_combined) * (1 - spec_combined);
  
  // group observed test results of combined panel.
  for (k in 1:k_sample) {
    // y_sample_combined[n] ~ bernoulli(p_sample_combined);
    // weight this part of the model as more important than the components.
    target += n_test * bernoulli_lpmf(y_sample_combined[k] | p_sample_combined);
  }
  
  for (n in 1:n_test) {
    
    // test-wise expected value of apparent prevalence of each component
    real p_sample = p[n] * sens[n] + (1 - p[n]) * (1 - spec[n]);
    
    // test-wise sensitivity and specificity priors of each component.
    target += binomial_lpmf(fp_spec[n] | k_spec[n], 1-spec[n]);
    target += binomial_lpmf(fn_sens[n] | k_sens[n], 1-sens[n]);
    
    // test-wise observed test results for each component
    for (k in 1:k_sample) {
      target += bernoulli_lpmf(y_sample[k,n] | p_sample);
    }
  }
}
