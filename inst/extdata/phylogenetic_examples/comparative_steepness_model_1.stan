functions {
  // functions needed
  // ProbFunction: adapted from Goffe et al 2018
  //   returns ratings given the start values, k and interaction outcomes
  // cum_winprob: 
  //   calculate cumulative winning probs given the final ratings, k and interaction outcomes
  // cumwinprob2steep:
  //   calculates steepness from cumulative winning probs
  
  real[] ProbFunction(vector start_vals, real k, int n_interactions, int n_ids, int[] winner_index, int[] loser_index) {
    real result[n_interactions];
    real to_add;
    vector[n_ids] cur_rating = start_vals;
    for (i in 1:n_interactions) {
      // centering:
      cur_rating = cur_rating - mean(cur_rating);
      // likelihood:
      result[i] = 1/(1 + exp(cur_rating[loser_index[i]] - cur_rating[winner_index[i]]));
      // update value:
      to_add = (1 - result[i]) * k;
      // update two ratings:
      cur_rating[winner_index[i]] = cur_rating[winner_index[i]] + to_add;
      cur_rating[loser_index[i]] = cur_rating[loser_index[i]] - to_add;
      }
    return result;
  }

  vector cum_winprob(vector start_vals, real k, int n_interactions, int n_ids, int[] winner_index, int[] loser_index) {
    real single_wp;
    real to_add;
    matrix[n_ids, n_ids] pairwise_winprobs;
    vector[n_ids] cumwinprobs;

    vector[n_ids] cur_rating = start_vals;
    for (i in 1:n_interactions) {
      single_wp = 1/(1 + exp(cur_rating[loser_index[i]] - cur_rating[winner_index[i]]));
      // update value:
      to_add = (1 - single_wp) * k;
      // update two ratings:
      cur_rating[winner_index[i]] = cur_rating[winner_index[i]] + to_add;
      cur_rating[loser_index[i]] = cur_rating[loser_index[i]] - to_add;
    }
    
    // pairwise winprobs
    for (i in 1:(n_ids - 1)) {
      for (j in (i + 1):n_ids) {
        single_wp = 1/(1 + exp(cur_rating[i] - cur_rating[j]));
        pairwise_winprobs[j, i] = single_wp;
        pairwise_winprobs[i, j] = 1.0 - single_wp;
      }
    }
    // individual sums of winprobs
    for (i in 1:n_ids) {
      pairwise_winprobs[i, i] = 0.0;
      cumwinprobs[i] = sum(pairwise_winprobs[i, ]);
    }
    return cumwinprobs;
  }

  vector cumwinprob2steep(vector nds, data int n_ids) {
    // aux objects for slope
    vector[n_ids] A1;
    vector[n_ids] A2;
    vector[n_ids] B1;
    vector[n_ids] AB1;
    // aux objects for intercept
    real sum_y;
    real sum_x2;
    real sum_x;
    real sum_xy;
    // aux for ranks
    vector[n_ids] theranks;
    int r;
    int s;
    // results
    vector[2] xsteep;

    // do the ranks
    for (i in 1:n_ids) {
      r = 1;
      s = 1;
      for (j in 1:i) {
        if (nds[j] < nds[i]) {
          r = r + 1;
        }
        if (nds[j] == nds[i]) {
          s = s + 1;
        }
      }

      for (j in (i + 1) : n_ids) {
        if (nds[j] < nds[i]) {
          r = r + 1;
        }
        if (nds[j] == nds[i]) {
          s = s + 1;
        }
      }
      theranks[i] = r + (s - 1) * 0.5 - 0.5;
    }

    // the intercept
    sum_y = sum(nds);
    sum_x2 = 0.0;
    for (i in 1:n_ids) {
      sum_x2 = sum_x2 + theranks[i] ^ 2;
    }
    sum_x = sum(theranks);
    sum_xy = sum(nds .* theranks);

    xsteep[1] = ((sum_y * sum_x2) - (sum_x * sum_xy)) / (((n_ids * sum_x2) - sum_x ^ 2));

    // the slope
    A1 = theranks - mean(theranks);
    B1 = nds - mean(nds);
    AB1 = A1 .* B1;
    A2 = A1 .* A1;
    xsteep[2] = sum(AB1) / sum(A2);

    return xsteep;
  }
}

data {
  // global information
  int<lower=1> n_datasets; // total number of data sets (interaction matrices)
  int<lower=1> n_total_interactions; // total number of interactions accross all data sets
  int<lower=1> n_total_ids; // total number of individuals accross all data sets

  // per-dataset-info
  int<lower=1> n_ids_per_dataset[n_datasets]; // ids per dataset
  int<lower=1> n_interactions_per_dataset[n_datasets]; // interactions per dataset
  
  // navigating
  int<lower=1> interaction_index[n_total_interactions];
  int<lower=1> individual_index[n_total_ids];
  int<lower=1> index_dataset_interactions_start[n_datasets]; // where do interactions start for a given dataset
  int<lower=1> index_individuals_start[n_datasets]; // start points for individuals for all-individals vector

  // interaction data
  int<lower=1> winner[n_total_interactions, 1]; // winner's index (within data set); columns pertain to n_rand (not yet suported)
  int<lower=1> loser[n_total_interactions, 1]; // losers's index (within data set); columns pertain to n_rand (not yet suported)
  
  // misc
  int<lower=0> y[n_total_interactions]; // outcome, i.e. winner always wins -> all values are 1

  // data for repeated measurements
  int<lower=1> N_repeated; // number of grouping levels for repeated measures
  int<lower=1> N_coeff_repeated; //number of coefficients (here: 1 [SD])
  int<lower=1> index_rep_measures[n_datasets]; // grouping indicator per observation
  vector[n_datasets] rep_measures_predictor; // group-level predictor values (intercept only)
  
  // data for phylogeny related information
  int<lower=1> N_phyl;
  int<lower=1> N_coeff_phyl;
  int<lower=1> spec_index[n_datasets]; 
  matrix[N_phyl, N_phyl] chol_mat; // cholesky factor of phylogenetic correlation matrix
  vector[n_datasets] phyl_predictor; // group-level predictor values (intercept only)
}

parameters {
  matrix[1, n_total_ids] EloStart_raw; // rows pertain to n_rand (not yet suported)
  vector<lower=0>[n_datasets] k_values; // k per data set
  
  // for overall beta model (steepness)
  real Intercept;
  real<lower=0> phi;
  
  vector<lower=0>[N_coeff_repeated] sd_repeated_measures; // group-level standard deviations for repeated measures
  vector[N_repeated] blups_repeated_measures[N_coeff_repeated]; // BLUPS for repeated measurements

  vector<lower=0>[N_coeff_phyl] sd_phyl; // group-level standard deviations for phylogenetic effects
  vector[N_phyl] blups_phyl[N_coeff_phyl];  // BLUPS for phylogenetic effects
}

transformed parameters {
  matrix[1, n_total_ids] EloStart; // rows pertain to n_rand (not yet suported)
  vector[N_repeated] eff_repeated;  // actual group-level effects
  vector[N_phyl] eff_phyl;  // actual group-level effects

  for (d in 1:n_datasets) {
    int idx[n_ids_per_dataset[d]];
    idx = segment(individual_index, index_individuals_start[d], n_ids_per_dataset[d]);
    EloStart[1, idx] = EloStart_raw[1, idx] - mean(EloStart_raw[1, idx]);
  }

  eff_repeated = (sd_repeated_measures[1] * (blups_repeated_measures[1]));
  eff_phyl = (sd_phyl[1] * (chol_mat * blups_phyl[1]));
}

model {
  // initiate vectors for beta model
  vector[n_datasets] steepness_as_response;
  vector[n_datasets] mu = Intercept + rep_vector(0.0, n_datasets);
  
  // estimate steepness for each data set
  for (d in 1:n_datasets) {
    int idx_individual[n_ids_per_dataset[d]];
    int idx_interaction[n_interactions_per_dataset[d]];
    vector[n_ids_per_dataset[d]] cumwinprobs;

    k_values[d] ~ normal(0, 1);
    idx_individual = segment(individual_index, index_individuals_start[d], n_ids_per_dataset[d]);
    idx_interaction = segment(interaction_index, index_dataset_interactions_start[d], n_interactions_per_dataset[d]);

    EloStart_raw[1, idx_individual] ~ normal(0, 1);
    
    y[idx_interaction] ~ bernoulli(ProbFunction(to_vector(EloStart[1, idx_individual]), 
                                                k_values[d],
                                                n_interactions_per_dataset[d],
                                                n_ids_per_dataset[d],
                                                winner[idx_interaction, 1], 
                                                loser[idx_interaction, 1]));
    
    cumwinprobs = cum_winprob(to_vector(EloStart[1, idx_individual]),
                              k_values[d],
                              n_interactions_per_dataset[d],
                              n_ids_per_dataset[d],
                              winner[idx_interaction, 1], 
                              loser[idx_interaction, 1]);
    
    steepness_as_response[d] = cumwinprob2steep(to_vector(cumwinprobs), 
                                                n_ids_per_dataset[d])[2];
  }
  
  // initialize linear predictor term
  for (n in 1:n_datasets) {
    mu[n] += eff_repeated[index_rep_measures[n]] * rep_measures_predictor[n] 
          + eff_phyl[spec_index[n]] * phyl_predictor[n];
  }
  for (n in 1:n_datasets) {
    // inverse link function
    mu[n] = inv_logit(mu[n]);
  }
  
  // likelihood for steepness model
  steepness_as_response ~ beta(mu * phi, (1.0 - mu) * phi);
  
  // priors
  Intercept ~ student_t(3, 0, 2.5);
  phi ~ gamma(0.01, 0.01);
  sd_phyl ~ student_t(3, 0, 2.5);
  blups_phyl[1] ~ normal(0, 1);
  sd_repeated_measures ~ student_t(3, 0, 2.5);
  blups_repeated_measures[1] ~ normal(0, 1);
}

generated quantities {
  // return actual steepness values per data set in the output
  vector[n_datasets] steepness;
  for (d in 1:n_datasets) {
    int idx[n_ids_per_dataset[d]];
    int idx_interaction[n_interactions_per_dataset[d]];
    vector[n_ids_per_dataset[d]] cumwinprobs;
    idx = segment(individual_index, index_individuals_start[d], n_ids_per_dataset[d]);
    idx_interaction = segment(interaction_index, index_dataset_interactions_start[d], n_interactions_per_dataset[d]);

    cumwinprobs = cum_winprob(to_vector(EloStart[1, idx]), 
                                          k_values[d], 
                                          n_interactions_per_dataset[d],
                                          n_ids_per_dataset[d], 
                                          winner[idx_interaction, 1], 
                                          loser[idx_interaction, 1]);
    steepness[d] = cumwinprob2steep(to_vector(cumwinprobs), n_ids_per_dataset[d])[2];
  }

}
