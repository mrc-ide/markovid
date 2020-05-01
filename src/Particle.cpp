
#include "Particle.h"

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(System &s, double beta) {
  
  // pointer to system object
  this->s_ptr = &s;
  
  // thermodynamic power
  this->beta = beta;
  
  // parameters
  d = s_ptr->d;
  
  // spline parameters
  node_y = vector<vector<double>>(s_ptr->n_sitrep, vector<double>(s_ptr->n_node));
  
  // population proportions
  scale_rel_prop = vector<double>(s_ptr->n_age_sitrep);
  
  // transition probabilities
  p_AI = vector<double>(s_ptr->n_age_indlevel);
  p_AD = vector<double>(s_ptr->n_age_indlevel);
  p_ID = vector<double>(s_ptr->n_age_indlevel);
  
  // mean durations
  m_AI = vector<double>(s_ptr->n_age_indlevel);
  m_AD = vector<double>(s_ptr->n_age_indlevel);
  m_AC = vector<double>(s_ptr->n_age_indlevel);
  m_ID = vector<double>(s_ptr->n_age_indlevel);
  m_IS = vector<double>(s_ptr->n_age_indlevel);
  m_SC = vector<double>(s_ptr->n_age_indlevel);
  
  // coefficient of variation of durations
  s_AI = vector<double>(s_ptr->n_age_indlevel);
  s_AD = vector<double>(s_ptr->n_age_indlevel);
  s_AC = vector<double>(s_ptr->n_age_indlevel);
  s_ID = vector<double>(s_ptr->n_age_indlevel);
  s_IS = vector<double>(s_ptr->n_age_indlevel);
  s_SC = vector<double>(s_ptr->n_age_indlevel);
  
  // dynamic lookup tables for interval distributions
  density_AL = vector<double>(s_ptr->lookup_max);
  density_AI = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  density_AD = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  density_AC = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  density_ID = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  density_IS = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  density_SC = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  
  // dynamic lookup tables for complementary cumulative density (ccdf) distributions
  tail_AI = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  tail_AD = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  tail_AC = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  tail_ID = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  tail_IS = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  tail_SC = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  
  // theta is the parameter vector in natural space
  theta = s_ptr->theta_init;
  theta_prop = vector<double>(d);
  
  // phi is a vector of transformed parameters
  phi = vector<double>(d);
  theta_to_phi();
  phi_prop = vector<double>(d);
  
  // run through likelihood calculation once for each parameter to initialise lookup tables
  //for (int i = 0; i < d; ++i) {
  //  double tmp = get_loglike(theta, i, true);
  //}
  //print_matrix(density_AI);
  //Rcpp::stop("header");
  
  // proposal parameters
  bw = vector<double>(d, 1.0);
  bw_index = vector<int>(d, 1);
  bw_stepsize = 1.0;
  
  // likelihoods and priors
  loglike = get_loglike(theta, 0);
  loglike_prop = 0;
  logprior = get_logprior(theta, 0);
  logprior_prop = 0;
  
  // acceptance rates
  accept_count = 0;
}

//------------------------------------------------
// transform phi_prop to theta_prop. See main.R for a key to transformation
// types
void Particle::phi_prop_to_theta_prop(int i) {
  
  switch(s_ptr->trans_type[i]) {
  case 0:
    theta_prop[i] = phi_prop[i];
    break;
  case 1:
    theta_prop[i] = s_ptr->theta_max[i] - exp(phi_prop[i]);
    break;
  case 2:
    theta_prop[i] = exp(phi_prop[i]) + s_ptr->theta_min[i];
    break;
  case 3:
    theta_prop[i] = (s_ptr->theta_max[i]*exp(phi_prop[i]) + s_ptr->theta_min[i]) / (1 + exp(phi_prop[i]));
    break;
  default:
    Rcpp::stop("trans_type invalid");
  }
  
}

//------------------------------------------------
// transform theta to phi. See main.R for a key to transformation types
void Particle::theta_to_phi() {
  
  for (int i = 0; i < d; ++i) {
    switch(s_ptr->trans_type[i]) {
    case 0:
      phi[i] = theta[i];
      break;
    case 1:
      phi[i] = log(s_ptr->theta_max[i] - theta[i]);
      break;
    case 2:
      phi[i] = log(theta[i] - s_ptr->theta_min[i]);
      break;
    case 3:
      phi[i] = log(theta[i] - s_ptr->theta_min[i]) - log(s_ptr->theta_max[i] - theta[i]);
      break;
    default:
      Rcpp::stop("trans_type invalid");
    }
  }
  
}

//------------------------------------------------
// get adjustment factor to account for reparameterisation
double Particle::get_adjustment(int i) {
  
  double ret = 0;
  switch(s_ptr->trans_type[i]) {
  case 0:
    // (no adjustment needed)
    break;
  case 1:
    ret = log(theta_prop[i] - s_ptr->theta_max[i]) - log(theta[i] - s_ptr->theta_max[i]);
    break;
  case 2:
    ret = log(theta_prop[i] - s_ptr->theta_min[i]) - log(theta[i] - s_ptr->theta_min[i]);
    break;
  case 3:
    ret = log(s_ptr->theta_max[i] - theta_prop[i]) + log(theta_prop[i] - s_ptr->theta_min[i]) - log(s_ptr->theta_max[i] - theta[i]) - log(theta[i] - s_ptr->theta_min[i]);
    break;
  default:
    Rcpp::stop("trans_type invalid");
  }
  return ret;
}

//------------------------------------------------
void Particle::update() {
  
  // set theta_prop and phi_prop to current values of theta and phi
  theta_prop = theta;
  phi_prop = phi;
  
  // loop through parameters
  for (int i = 0; i < d; ++i) {
    if (s_ptr->skip_param[i]) {
      continue;
    }
    
    // generate new phi_prop[i]
    phi_prop[i] = rnorm1(phi[i], bw[i]);
    
    // transform phi_prop[i] to theta_prop[i]
    phi_prop_to_theta_prop(i);
    
    // calculate adjustment factor, taking into account forwards and backwards
    // moves
    double adj = get_adjustment(i);
    
    // calculate likelihood and prior of proposed theta
    loglike_prop = get_loglike(theta_prop, i);
    logprior_prop = get_logprior(theta, i);
    
    // calculate Metropolis-Hastings ratio
    double MH = beta*(loglike_prop - loglike) + (logprior_prop - logprior) + adj;
    
    // accept or reject move
    bool MH_accept = (log(runif_0_1()) < MH);
    
    // implement changes
    if (MH_accept) {
      
      // update theta and phi
      theta[i] = theta_prop[i];
      phi[i] = phi_prop[i];
      
      // update likelihoods
      loglike = loglike_prop;
      logprior = logprior_prop;
      
      // Robbins-Monro positive update  (on the log scale)
      bw[i] = exp(log(bw[i]) + bw_stepsize*(1 - 0.234)/sqrt(bw_index[i]));
      bw_index[i]++;
      
      // add to acceptance rate count
      accept_count++;
      
    } else {
      
      // reset theta_prop and phi_prop
      theta_prop[i] = theta[i];
      phi_prop[i] = phi[i];
      
      // Robbins-Monro negative update (on the log scale)
      bw[i] = exp(log(bw[i]) - bw_stepsize*0.234/sqrt(bw_index[i]));
      bw_index[i]++;
      
    } // end MH step
    
  }  // end loop over parameters
    
}  // end update_univar function

//------------------------------------------------
// define cpp loglike function
double Particle::get_loglike(vector<double> &theta, int theta_i, bool init_lookup) {
  
  // ----------------------------------------------------------------
  // unpack parameters and define fixed/derived parameters
  
  // admissions spline y values
  int pi = 0;
  for (int i = 0; i < s_ptr->n_sitrep; ++i) {
    for (int j = 0; j < s_ptr->n_node; ++j) {
      node_y[i][j] = theta[pi++];
    }
  }
  
  // population proportions
  for (int i = 0; i < s_ptr->n_age_sitrep; ++i) {
    scale_rel_prop[i] = theta[pi++];
  }
  
  // transition probabilities
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    p_AI[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    p_AD[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    p_ID[i] = theta[pi++];
  }
  
  // mean durations
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    m_AI[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    m_AD[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    m_AC[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    m_ID[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    m_IS[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    m_SC[i] = theta[pi++];
  }
  
  // coefficient of variation of durations
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    s_AI[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    s_AD[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    s_AC[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    s_ID[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    s_IS[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    s_SC[i] = theta[pi++];
  }
  
  // admission to labtest
  double m_AL = theta[pi++];
  double rho = 1.0 / m_AL;
  
  // misc parameters
  double scale_p_AI = theta[pi++];
  double scale_p_AD = theta[pi++];
  double scale_p_ID = theta[pi++];
  
  // ----------------------------------------------------------------
  // update lookup tables
  
#ifdef LOOK
  // update interval distribution lookup tables
  if (s_ptr->update_param[theta_i] == 0) {
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      density_AL[j] = R::pgamma(j + 1, 1.0, m_AL, true, false) -
                      R::pgamma(j, 1.0, m_AL, true, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 1) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      density_AI[i][j] = R::pgamma(j + 1, 1.0/(s_AI[i]*s_AI[i]), m_AI[i]*s_AI[i]*s_AI[i], true, false) -
                         R::pgamma(j, 1.0/(s_AI[i]*s_AI[i]), m_AI[i]*s_AI[i]*s_AI[i], true, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 2) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      density_AD[i][j] = R::pgamma(j + 1, 1.0/(s_AD[i]*s_AD[i]), m_AD[i]*s_AD[i]*s_AD[i], true, false) -
                         R::pgamma(j, 1.0/(s_AD[i]*s_AD[i]), m_AD[i]*s_AD[i]*s_AD[i], true, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 3) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      density_AC[i][j] = R::pgamma(j + 1, 1.0/(s_AC[i]*s_AC[i]), m_AC[i]*s_AC[i]*s_AC[i], true, false) -
                         R::pgamma(j, 1.0/(s_AC[i]*s_AC[i]), m_AC[i]*s_AC[i]*s_AC[i], true, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 4) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      density_ID[i][j] = R::pgamma(j + 1, 1.0/(s_ID[i]*s_ID[i]), m_ID[i]*s_ID[i]*s_ID[i], true, false) -
                         R::pgamma(j, 1.0/(s_ID[i]*s_ID[i]), m_ID[i]*s_ID[i]*s_ID[i], true, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 5) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      density_IS[i][j] = R::pgamma(j + 1, 1.0/(s_IS[i]*s_IS[i]), m_IS[i]*s_IS[i]*s_IS[i], true, false) -
                         R::pgamma(j, 1.0/(s_IS[i]*s_IS[i]), m_IS[i]*s_IS[i]*s_IS[i], true, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 6) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      density_SC[i][j] = R::pgamma(j + 1, 1.0/(s_SC[i]*s_SC[i]), m_SC[i]*s_SC[i]*s_SC[i], true, false) -
                         R::pgamma(j, 1.0/(s_SC[i]*s_SC[i]), m_SC[i]*s_SC[i]*s_SC[i], true, false);
    }
  }
  
  // update ccdf lookup tables
  if (s_ptr->update_param[theta_i] == 1) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      tail_AI[i][j] = R::pgamma(j + 1, 1.0/(s_AI[i]*s_AI[i]), m_AI[i]*s_AI[i]*s_AI[i], false, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 2) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      tail_AD[i][j] = R::pgamma(j + 1, 1.0/(s_AD[i]*s_AD[i]), m_AD[i]*s_AD[i]*s_AD[i], false, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 3) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      tail_AC[i][j] = R::pgamma(j + 1, 1.0/(s_AC[i]*s_AC[i]), m_AC[i]*s_AC[i]*s_AC[i], false, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 4) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      tail_ID[i][j] = R::pgamma(j + 1, 1.0/(s_ID[i]*s_ID[i]), m_ID[i]*s_ID[i]*s_ID[i], false, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 5) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      tail_IS[i][j] = R::pgamma(j + 1, 1.0/(s_IS[i]*s_IS[i]), m_IS[i]*s_IS[i]*s_IS[i], false, false);
    }
  }
  if (s_ptr->update_param[theta_i] == 6) {
    int i = s_ptr->update_age[theta_i] - 1;
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      tail_SC[i][j] = R::pgamma(j + 1, 1.0/(s_SC[i]*s_SC[i]), m_SC[i]*s_SC[i]*s_SC[i], false, false);
    }
  }
#endif
  
  // force update
  bool force_update = true;
  
  if (force_update) {
    
    for (int j = 0; j < s_ptr->lookup_max; ++j) {
      //density_AL[j] = get_density_gamma(m_AL, 1.0, j);
    }
    for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
      for (int j = 0; j < s_ptr->lookup_max; ++j) {
        density_AI[i][j] = get_density_gamma(m_AI[i], s_AI[i], j);
        density_AD[i][j] = get_density_gamma(m_AD[i], s_AD[i], j);
        density_AC[i][j] = get_density_gamma(m_AC[i], s_AC[i], j);
        density_ID[i][j] = get_density_gamma(m_ID[i], s_ID[i], j);
        density_IS[i][j] = get_density_gamma(m_IS[i], s_IS[i], j);
        density_SC[i][j] = get_density_gamma(m_SC[i], s_SC[i], j);
      }
    }
    
    for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
      for (int j = 0; j < s_ptr->lookup_max; ++j) {
        tail_AI[i][j] = get_tail_gamma(m_AI[i], s_AI[i], j);
        tail_AD[i][j] = get_tail_gamma(m_AD[i], s_AD[i], j);
        tail_AC[i][j] = get_tail_gamma(m_AC[i], s_AC[i], j);
        tail_ID[i][j] = get_tail_gamma(m_ID[i], s_ID[i], j);
        tail_IS[i][j] = get_tail_gamma(m_IS[i], s_IS[i], j);
        tail_SC[i][j] = get_tail_gamma(m_SC[i], s_SC[i], j);
      }
    }
    
  }  // end force update
  
  // option to return after initialising lookups
  if (init_lookup) {
    return 0.0;
  }
  
  // ----------------------------------------------------------------
  // individual level component of likelihood
  
  // sum log-likelihood over all data
  double ret = 0.0;
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    
    // get age-specific parameters
    int age_i = s_ptr->age_group[i] - 1;
    double q_AI = p_AI[age_i];// * scale_p_AI;
    double q_AD = p_AD[age_i];// * scale_p_AD;
    double q_ID = p_ID[age_i];// * scale_p_ID;
    
    // probabilities cannot exceed 1.0
    if (q_AI > 1.0 || q_AD > 1.0 || q_ID > 1.0) {
      ret = -std::numeric_limits<double>::infinity();
    }
    
    // non-icu pathway
    if (s_ptr->icu[i] == 0) {
      
      // admission to death
      if (s_ptr->final_outcome[i] == 1) {
        int delta = s_ptr->date_final_outcome[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * q_AD * density_AD[age_i][delta] );
        //ret += log( (1.0 - q_AI) * q_AD * get_density_gamma(m_AD[age_i], s_AD[age_i], delta) );
        
        //if (!std::isfinite(ret)) {
        //  Rcpp::Rcout << q_AI << " " << q_AD << " " << delta << " " << age_i << " " << density_AD[age_i][delta] << "\n";
        //  print_vector(density_AD[age_i]);
        //  Rcpp::stop("foo1");
        //}
        
      }
      
      // admission to discharge
      if (s_ptr->final_outcome[i] == 2) {
        int delta = s_ptr->date_final_outcome[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * (1.0 - q_AD) * density_AC[age_i][delta] );
        //ret += log( (1.0 - q_AI) * (1.0 - q_AD) * get_density_gamma(m_AC[age_i], s_AC[age_i], delta) );
      }
      
      // open case
      if (s_ptr->final_outcome[i] == -1) {
        int delta = s_ptr->date_censor[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * q_AD * tail_AD[age_i][delta] +
                    (1.0 - q_AI) * (1.0 - q_AD) * tail_AC[age_i][delta] );
        //ret += log( (1.0 - q_AI) * q_AD * get_tail_gamma(m_AD[age_i], s_AD[age_i], delta) +
        //  (1.0 - q_AI) * (1.0 - q_AD) * get_tail_gamma(m_AC[age_i], s_AC[age_i], delta) );
      }
      
    }  // end non-icu pathway
    
    // icu pathway
    if (s_ptr->icu[i] == 1) {
      
      // admission to icu
      {
        int delta = s_ptr->date_icu[i] - s_ptr->date_admission[i];
        ret += log( q_AI * density_AI[age_i][delta] );
        //ret += log( q_AI * get_density_gamma(m_AI[age_i], s_AI[age_i], delta) );
      }
      
      // non-stepdown pathway
      if (s_ptr->stepdown[i] == 0) {
        
        // icu to death
        if (s_ptr->final_outcome[i] == 1) {
          int delta = s_ptr->date_final_outcome[i] - s_ptr->date_icu[i];
          ret += log( q_AI * q_ID * density_ID[age_i][delta] );
          //ret += log( q_AI * q_ID * get_density_gamma(m_ID[age_i], s_ID[age_i], delta) );
        }
        
        // open case in icu
        if (s_ptr->final_outcome[i] == -1) {
          int delta = s_ptr->date_censor[i] - s_ptr->date_icu[i];
          ret += log( q_AI * q_ID * tail_ID[age_i][delta] +
                      q_AI * (1.0 - q_ID) * tail_IS[age_i][delta] );
          //ret += log( q_AI * q_ID * get_tail_gamma(m_ID[age_i], s_ID[age_i], delta) +
          //            q_AI * (1.0 - q_ID) * get_tail_gamma(m_IS[age_i], s_IS[age_i], delta) );
        }
        
      }  // end non-stepdown pathway
      
      // stepdown pathway
      if (s_ptr->stepdown[i] == 1) {
        
        // icu to stepdown
        {
          int delta = s_ptr->date_stepdown[i] - s_ptr->date_icu[i];
          ret += log( q_AI * (1.0 - q_ID) * density_IS[age_i][delta] );
          //ret += log( q_AI * (1.0 - q_ID) * get_density_gamma(m_IS[age_i], s_IS[age_i], delta) );
        }
        
        // stepdown to discharge
        if (s_ptr->final_outcome[i] == 2) {
          int delta = s_ptr->date_final_outcome[i] - s_ptr->date_stepdown[i];
          ret += log( q_AI * (1.0 - q_ID) * density_SC[age_i][delta] );
          //ret += log( q_AI * (1.0 - q_ID) * get_density_gamma(m_SC[age_i], s_SC[age_i], delta) );
        }
        
        // open case in stepdown
        if (s_ptr->final_outcome[i] == -1) {
          int delta = s_ptr->date_censor[i] - s_ptr->date_stepdown[i];
          ret += log( q_AI * (1.0 - q_ID) * tail_SC[age_i][delta] );
          //ret += log( q_AI * (1.0 - q_ID) * get_tail_gamma(m_SC[age_i], s_SC[age_i], delta) );
        }
        
      }  // end stepdown pathway
      
    }  // end icu pathway
    
  }  // end i loop
  
  // ----------------------------------------------------------------
  // return
  
  //print(ret);
  
  // catch underflow
  if (!std::isfinite(ret)) {
    ret = -DBL_MAX/100.0;
  }
  
  return ret;
}

//------------------------------------------------
// define cpp logprior function
double Particle::get_logprior(vector<double> &theta, int theta_i) {
  return 0.0;
}

//------------------------------------------------
double Particle::get_density_gamma(double m, double s, int i) {
  int m_index = ceil(m / 30.0 * s_ptr->lookup_n_m) - 1;
  int s_index = ceil(s * s_ptr->lookup_n_s) - 1;
  if (m_index >= s_ptr->density_gamma.size() ||
      s_index >= s_ptr->density_gamma[m_index].size() ||
      i >= s_ptr->density_gamma[m_index][s_index].size()) {
    //print("\n", m, m_index, s, s_index, i);
    //Rcpp::stop("error in get_density_gamma(): outside lookup range");
    return 0.0;
  }
  double ret = s_ptr->density_gamma[m_index][s_index][i];
  return ret;
}

//------------------------------------------------
double Particle::get_tail_gamma(double m, double s, int i) {
  double ret = s_ptr->tail_gamma[0][0][0];
  return ret;
}
