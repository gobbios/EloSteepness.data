data {
  int<lower=1> N; // total number of observations
  vector[N] the_steepness; // response variable
  
  // data for repeated measurements
  int<lower=1> N_repeated;  // number of grouping levels for repeated measures
  int<lower=1> N_coeff_repeated; // number of coefficients (here: 1)
  int<lower=1> index_rep_measures[N]; // indexing per observation
  vector[N] rep_measures_predictor; // group-level predictor (intercept only)
  
  // data for phylogeny related information
  int<lower=1> N_phyl; // number of grouping levels for phylogney (same as above...)
  int<lower=1> N_coeff_phyl; // number of coefficients (here: 1)
  int<lower=1> spec_index[N]; // indexing per observation
  matrix[N_phyl, N_phyl] chol_mat; // cholesky factor of phylogenetic correlation matrix
  vector[N] phyl_predictor; // group-level predictor values (interecept only)
}

parameters {
  // for overall beta model
  real Intercept;
  real<lower=0> phi;
  
  // specific to repeated measurements
  vector<lower=0>[N_coeff_repeated] sd_repeated_measures;
  vector[N_repeated] blups_repeated_measures[N_coeff_repeated];
  
  // specific to phylogeny
  vector<lower=0>[N_coeff_phyl] sd_phyl;
  vector[N_phyl] blups_phyl[N_coeff_phyl];
}

transformed parameters {
  vector[N_repeated] eff_repeated;  // actual group-level effects
  vector[N_phyl] eff_phyl;  // actual group-level effects
  
  eff_phyl = (sd_phyl[1] * (chol_mat * blups_phyl[1]));
  eff_repeated = (sd_repeated_measures[1] * (blups_repeated_measures[1]));

}

model {
  vector[N] mu = Intercept + rep_vector(0.0, N); // initialize intercept
  for (n in 1:N) {
    mu[n] += eff_repeated[index_rep_measures[n]] * rep_measures_predictor[n] + eff_phyl[spec_index[n]] * phyl_predictor[n];
  }
  for (n in 1:N) {
    mu[n] = inv_logit(mu[n]);
  }
  // the model
  the_steepness ~ beta(mu * phi, (1 - mu) * phi);
  // priors
  Intercept ~ student_t(3, 0, 2.5);
  phi ~ gamma(0.01, 0.01);
  sd_phyl ~ student_t(3, 0, 2.5);
  blups_phyl[1] ~ normal(0, 1);
  sd_repeated_measures ~ student_t(3, 0, 2.5);
  blups_repeated_measures[1] ~ normal(0, 1);
}
generated quantities {
  // just output sample for priors of SDs
  real prior_sd_phyl = student_t_rng(3, 0, 2.5);
  real prior_sd_repeated_measures = student_t_rng(3, 0, 2.5);
  // rejection sampling for SDs
  while (prior_sd_phyl < 0) {
    prior_sd_phyl = student_t_rng(3,0,2.5);
  }
  while (prior_sd_repeated_measures < 0) {
    prior_sd_repeated_measures = student_t_rng(3, 0, 2.5);
  }
}

