data {
  int<lower = 1> n_test;
  
  // samples are rows; tests are columns although this is transposed in R
  // a positive result is 1 a negative result is 0
  int<lower = 1> k_sample;
  int<lower = 0,upper = 1> y_sample[k_sample,n_test];
  
  // negative disease free controls (for sensitivity)
  // a positive result is 1 a negative result is 0
  int<lower = 0> k_control;
  int<lower = 0,upper = 1> y_control[k_control,n_test];
  
  // positive disease present controls (for sensitivity)
  // a positive result is 1 a negative result is 0
  // int<lower = 0> n_disease;
  // int<lower = 0,upper = 1> y_disease[n_disease,n_test];
  
  // specificity prior hyper-parameters
  // a is true negatives
  // b is false positives
  real<lower = 0> a_spec;
  real<lower = 0> b_spec;
  
  // sensitivity prior hyper-parameters
  // a is true positives
  // b is false negatives
  real<lower = 0> a_sens;
  real<lower = 0> b_sens;
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
  vector<lower=0.75, upper = 1>[n_test] spec;
  // The sensitivity of the component test
  vector<lower=0.5, upper = 1>[n_test] sens;
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
    sens[n] ~ beta(a_sens, b_sens);
    spec[n] ~ beta(a_spec, b_spec);
    
    // test-wise specificity control group results for each component (optional)
    for (k in 1:k_control) {
      target += bernoulli_lpmf(y_control[k,n] | 1-spec[n]);
    }
    
    // test-wise sensitivity control group results for each component (optional)
    // TODO: needs rethink as different subtypes will not all be positive in a 
    // supertype disease positive cohort. To get info about sensitivity we need
    // an expected result per serotype.
    // need to multiply this by an indicator function per seroptype whcih means
    // we need to submit a full matrix of disease positive controls and serotype
    // status and check against that, if we do this at all.
    // for (n in 1:n_disease) {
    //   y_disease[n,t] ~ bernoulli(sens[t]);
    // }
    
    // test-wise observed test results for each component
    for (k in 1:k_sample) {
      target += bernoulli_lpmf(y_sample[k,n] | p_sample);
    }
  }
}
