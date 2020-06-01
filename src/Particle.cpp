
#include "Particle.h"

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(System &s) {
  
  // pointer to system object
  this->s_ptr = &s;
  
  // parameters
  d = s_ptr->d;
  
  // spline parameters
  node_y = vector<vector<double>>(s_ptr->n_region, vector<double>(s_ptr->n_node));
  
  // population proportions
  scale_rel_prop = vector<double>(s_ptr->n_age_sitrep);
  
  // rescaling parameters
  scale_p_AI = vector<double>(s_ptr->n_region);
  scale_p_AD = vector<double>(s_ptr->n_region);
  scale_p_ID = vector<double>(s_ptr->n_region);
  
  // lab test weights
  pos_on_day = vector<double>(s_ptr->lookup_max);
  neg_by_day = vector<double>(s_ptr->lookup_max);
  pos_by_day = vector<double>(s_ptr->lookup_max);
  
  // transition probabilities
  p_AI = vector<double>(s_ptr->n_age_indlevel);
  p_AD = vector<double>(s_ptr->n_age_indlevel);
  p_ID = vector<double>(s_ptr->n_age_indlevel);
  
  // mean durations
  m_AC = vector<double>(s_ptr->n_age_indlevel);
  
  // lookup tables for interval distributions
  density_AL = vector<double>(s_ptr->lookup_max);
  density_AI = vector<double>(s_ptr->lookup_max);
  density_AD = vector<double>(s_ptr->lookup_max);
  density_AC = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  density_ID = vector<double>(s_ptr->lookup_max);
  density_IS = vector<double>(s_ptr->lookup_max);
  density_SC = vector<double>(s_ptr->lookup_max);
  
  // lookup tables for complementary cumulative density (ccdf) distributions
  tail_AI = vector<double>(s_ptr->lookup_max);
  tail_AD = vector<double>(s_ptr->lookup_max);
  tail_AC = vector<vector<double>>(s_ptr->n_age_indlevel, vector<double>(s_ptr->lookup_max));
  tail_ID = vector<double>(s_ptr->lookup_max);
  tail_IS = vector<double>(s_ptr->lookup_max);
  tail_SC = vector<double>(s_ptr->lookup_max);
  
  // objects for storing progression
  admission_incidence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_spline)));
  deaths_incidence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_spline)));
  discharges_incidence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_spline)));
  general_prevalence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_spline)));
  critical_prevalence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_spline)));
  
  delta_stepup = std::vector<double>(s_ptr->n_spline);
  delta_stepdown = std::vector<double>(s_ptr->n_spline);
  delta_deaths_general = std::vector<double>(s_ptr->n_spline);
  delta_discharges_general = std::vector<double>(s_ptr->n_spline);
  delta_open_general = std::vector<double>(s_ptr->n_spline);
  delta_deaths_critical = std::vector<double>(s_ptr->n_spline);
  delta_open_critical = std::vector<double>(s_ptr->n_spline);
  pos_by_day = std::vector<double>(s_ptr->n_spline);
  pos_on_day = std::vector<double>(s_ptr->n_spline);
  
  // theta is the parameter vector in natural space
  theta = s_ptr->theta_init;
  theta_prop = vector<double>(d);
  
  // phi is a vector of transformed parameters
  phi = vector<double>(d);
  theta_to_phi();
  phi_prop = vector<double>(d);
  
  // run through likelihood calculation once for each parameter to initialise lookup tables
  for (int i = 0; i < d; ++i) {
    double tmp = get_loglike(theta, i, true);
  }
  
  // proposal parameters
  bw = vector<double>(d, 1.0);
  bw_index = vector<int>(d, 1);
  bw_stepsize = 1.0;
  
  // likelihoods and priors
  loglike = get_loglike(theta, 0, false);
  loglike_prop = 0;
  logprior = get_logprior(theta, 0);
  logprior_prop = 0;
  
  //print(loglike);
  //Rcpp::stop("header");
  
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
void Particle::update(double beta) {
  
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
    loglike_prop = get_loglike(theta_prop, i, false);
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
      
      // Robbins-Monro positive update (on the log scale)
      bw[i] = exp(log(bw[i]) + bw_stepsize*(1 - 0.234)/sqrt(bw_index[i]));
      bw_index[i]++;
      
      // add to acceptance rate count
      accept_count++;
      
    } else {
      
      // reset theta_prop and phi_prop
      theta_prop[i] = theta[i];
      phi_prop[i] = phi[i];
      
      // reset lookup tables etc.
      double tmp = get_loglike(theta, i, false);
      
      // Robbins-Monro negative update (on the log scale)
      bw[i] = exp(log(bw[i]) - bw_stepsize*0.234/sqrt(bw_index[i]));
      bw_index[i]++;
      
    } // end MH step
    
  }  // end loop over parameters
    
}  // end update_univar function

//------------------------------------------------
// define cpp loglike function
double Particle::get_loglike(vector<double> &theta, int theta_i, bool quick_exit) {
  
  // ----------------------------------------------------------------
  // unpack parameters and define fixed/derived parameters
  
  // admissions spline y values
  int pi = 0;
  for (int i = 0; i < s_ptr->n_region; ++i) {
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
  double m_AI = theta[pi++];
  double m_AD = theta[pi++];
  for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
    m_AC[i] = theta[pi++];
  }
  double m_ID = theta[pi++];
  double m_IS = theta[pi++];
  double m_SC = theta[pi++];
  
  // coefficients of variation of durations
  double s_AI = theta[pi++];
  double s_AD = theta[pi++];
  double s_AC = theta[pi++];
  double s_ID = theta[pi++];
  double s_IS = theta[pi++];
  double s_SC = theta[pi++];
  
  // admission to labtest
  double m_AL = theta[pi++];
  double rho = 1.0 / m_AL;
  
  // rescaling parameters
  for (int i = 0; i < s_ptr->n_region; ++i) {
    scale_p_AI[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_region; ++i) {
    scale_p_AD[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_region; ++i) {
    scale_p_ID[i] = theta[pi++];
  }
  
  // get update rules
  int update_density = s_ptr->update_density[theta_i];
  int update_region = s_ptr->update_region[theta_i];
  int update_indlevel_age = s_ptr->update_indlevel_age[theta_i];
  int update_sitrep_age = s_ptr->update_sitrep_age[theta_i];
  
  // initialise loglikelihood
  double ret = 0.0;
  
  // ----------------------------------------------------------------
  // update density lookup tables
  
  if (update_density == 1) {
    update_gamma_density(density_AI, m_AI, s_AI);
    update_gamma_tail(tail_AI, m_AI, s_AI);
  }
  if (update_density == 2) {
    update_gamma_density(density_AD, m_AD, s_AD);
    update_gamma_tail(tail_AD, m_AD, s_AD);
  }
  if (update_density == 3) {
    if (update_indlevel_age == 0) {
      for (int i = 0; i < s_ptr->n_age_indlevel; ++i) {
        update_gamma_density(density_AC[i], m_AC[i], s_AC);
        update_gamma_tail(tail_AC[i], m_AC[i], s_AC);
      }
    } else {
      int i = update_indlevel_age - 1;
      update_gamma_density(density_AC[i], m_AC[i], s_AC);
      update_gamma_tail(tail_AC[i], m_AC[i], s_AC);
    }
  }
  if (update_density == 4) {
    update_gamma_density(density_ID, m_ID, s_ID);
    update_gamma_tail(tail_ID, m_ID, s_ID);
  }
  if (update_density == 5) {
    update_gamma_density(density_IS, m_IS, s_IS);
    update_gamma_tail(tail_IS, m_IS, s_IS);
  }
  if (update_density == 6) {
    update_gamma_density(density_SC, m_SC, s_SC);
    update_gamma_tail(tail_SC, m_SC, s_SC);
  }
  
  // ----------------------------------------------------------------
  // update progression objects
  
  // calculate lab test weights
  //for (int i = 0; i < s_ptr->n_spline; ++i) {
    //pos_by_day[i] = (i >= floor(m_AL)) ? 1.0 : 0.0;
    //pos_on_day[i] = (i == floor(m_AL)) ? 1.0 : 0.0;
    //pos_on_day[i] = R::pgamma(i + 1, 1.0/(1.0*1.0), m_AL*1.0*1.0, true, false) -
    //                R::pgamma(i, 1.0/(1.0*1.0), m_AL*1.0*1.0, true, false);
    //pos_by_day[i] = R::pgamma(i, 1.0/(1.0*1.0), m_AL*1.0*1.0, true, false);
  //}
  update_gamma_density(pos_on_day, m_AL, 1.0);
  update_gamma_tail(neg_by_day, m_AL, 1.0);
  for (int i = 0; i < s_ptr->lookup_max; ++i) {
    pos_by_day[i] = 1.0 - neg_by_day[i];
  }
  
  // loop through sitrep regions
  for (int region_i = 0; region_i < s_ptr->n_region; ++region_i) {
    //continue;
    
    // skip based on update rules
    if (update_region != (region_i + 1) && update_region != 0) {
      continue;
    }
    
    // get spline gradients between nodes
    vector<double> node_grad(s_ptr->n_node - 1);
    for (int i = 1; i < s_ptr->n_node; ++i) {
      node_grad[i-1] = double(node_y[region_i][i] - node_y[region_i][i-1]) / double(s_ptr->node_x[i] - s_ptr->node_x[i-1]);
    }
    
    // create spline
    vector<double> admissions_spline(s_ptr->n_spline);
    admissions_spline[0] = node_y[region_i][0];
    int node_j = 0;
    for (int i = 1; i < s_ptr->n_spline; ++i) {
      
      // update curve
      admissions_spline[i] = admissions_spline[i-1] + node_grad[node_j];
      
      // update node_j
      if ((s_ptr->node_x[0] + i) >= s_ptr->node_x[node_j+1]) {
        node_j++;
      }
    }
    
    // loop through sitrep age groups
    for (int age_i = 0; age_i < s_ptr->n_age_sitrep; ++age_i) {
      
      // skip based on update rules
      if (update_sitrep_age != (age_i + 1) && update_sitrep_age != 0) {
        continue;
      }
      
      // loop through individual-level age groups within this sitrep age group
      for (unsigned int age_j = 0; age_j < s_ptr->map_age_sitrep[age_i].size(); ++age_j) {
        
        // get individual-level age group
        int age_indlevel = s_ptr->map_age_sitrep[age_i][age_j] - 1;
        
        // get parameters specific to this age group
        double p_AI_i = p_AI[age_indlevel] * scale_p_AI[region_i];
        double p_AD_i = p_AD[age_indlevel] * scale_p_AD[region_i];
        double p_ID_i = p_ID[age_indlevel] * scale_p_ID[region_i];
        
        // probabilities cannot exceed 1.0
        if (p_AI_i > 1.0 || p_AD_i > 1.0 || p_ID_i > 1.0) {
          ret = -std::numeric_limits<double>::infinity();
        }
        
        // reset progression objects
        fill(admission_incidence[region_i][age_i].begin(), admission_incidence[region_i][age_i].end(), 0.0);
        fill(deaths_incidence[region_i][age_i].begin(), deaths_incidence[region_i][age_i].end(), 0.0);
        fill(discharges_incidence[region_i][age_i].begin(), discharges_incidence[region_i][age_i].end(), 0.0);
        fill(general_prevalence[region_i][age_i].begin(), general_prevalence[region_i][age_i].end(), 0.0);
        fill(critical_prevalence[region_i][age_i].begin(), critical_prevalence[region_i][age_i].end(), 0.0);
        
        // store incidence in stepup and stepdown care
        std::vector<double> stepup(s_ptr->n_spline);
        std::vector<double> stepdown(s_ptr->n_spline);
        
        // calculate unscaled progression vectors
        for (int i = 0; i < s_ptr->n_spline; ++i) {
          delta_stepup[i] = p_AI_i * density_AI[i];
          delta_stepdown[i] = (1 - p_ID_i) * density_IS[i];
          delta_deaths_general[i] = (1 - p_AI_i) * p_AD_i * density_AD[i];
          delta_discharges_general[i] =  (1 - p_AI_i) * (1 - p_AD_i) * density_AC[age_indlevel][i];
          delta_open_general[i] =  (1 - p_AI_i) * p_AD_i * tail_AD[i] +
                                   (1 - p_AI_i) * (1 - p_AD_i) * tail_AC[age_indlevel][i] +
                                    p_AI_i * tail_AI[i];
          delta_deaths_critical[i] = p_ID_i * density_ID[i];
          delta_open_critical[i] =  p_ID_i * tail_ID[i] +
                                    (1 - p_ID_i) * tail_IS[i];
        }
        
        // loop through each spline day in turn
        for (int i = 0; i < s_ptr->n_spline; ++i) {
          
          // get true admissions on this day from spline
          double true_admissions = scale_rel_prop[age_i] * s_ptr->rel_prop[age_indlevel] * exp(admissions_spline[i]);
          
          // project true admissions out over future days
          for (int j = i; j < s_ptr->n_spline; ++j) {
            
            // update incidence in stepup and stepdown
            stepup[j] += true_admissions * delta_stepup[j - i];
            stepdown[j] += stepup[i] * delta_stepdown[j - i];
            
            // deaths, discharges and open cases in general ward
            double deaths_general = true_admissions * delta_deaths_general[j - i];
            double discharges_general = true_admissions * delta_discharges_general[j - i];
            double open_general = true_admissions * delta_open_general[j - i];
            
            // deaths and open cases in critical care
            double deaths_critical = stepup[i] * delta_deaths_critical[j - i];
            double open_critical = stepup[i] * delta_open_critical[j - i];
            
            // discharges and open cases in stepdown
            double discharges_stepdown = stepdown[i] * density_SC[j - i];
            double open_stepdown = stepdown[i] * tail_SC[j - i];
            
            // delay from admission to testing
            double pos_by_day_j = pos_by_day[j - i];
            double pos_on_day_j = pos_on_day[j - i];
            
            // incidence of observed admission
            admission_incidence[region_i][age_i][j] += pos_on_day_j * (open_general + open_critical + open_stepdown);
            
            // incidence of death and discharge
            deaths_incidence[region_i][age_i][j] += pos_by_day_j * (deaths_general + deaths_critical);
            discharges_incidence[region_i][age_i][j] += pos_by_day_j * (discharges_general + discharges_stepdown);
            
            // prevalence in general and critical beds
            general_prevalence[region_i][age_i][j] += pos_by_day_j * (open_general + open_stepdown);
            critical_prevalence[region_i][age_i][j] += pos_by_day_j * open_critical;
            
          }  // end j loop over days
          
        }  // end i loop over days
        
      }  // end age_j loop
      
    }  // end age_i loop
    
  }  // end region_i loop
  
  // option to return at this stage
  if (quick_exit) {
    return 0.0;
  }
  
  // ----------------------------------------------------------------
  // individual-level component of likelihood
  
  // sum log-likelihood over individual-level data
  if (true) {
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    //continue;
    
    // get age group of this individual
    int age_i = s_ptr->age_group[i] - 1;
    
    // get age-specific parameters
    double q_AI = p_AI[age_i];
    double q_AD = p_AD[age_i];
    double q_ID = p_ID[age_i];
    
    // probabilities cannot exceed 1.0
    if (q_AI > 1.0 || q_AD > 1.0 || q_ID > 1.0) {
      ret = -std::numeric_limits<double>::infinity();
    }
    
    //if (i == 483) {
    //  print(ret, s_ptr->icu[i], s_ptr->final_outcome[i], s_ptr->stepdown[i]);
    //  Rcpp::stop("i exit");
    //}
    
    // non-icu pathway
    if (s_ptr->icu[i] == 0) {
      
      // admission to death
      if (s_ptr->final_outcome[i] == 1) {
        int delta = s_ptr->date_final_outcome[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * q_AD * density_AD[delta] );
        //ret += log( (1.0 - q_AI) * q_AD );
      }
      
      // admission to discharge
      if (s_ptr->final_outcome[i] == 2) {
        int delta = s_ptr->date_final_outcome[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * (1.0 - q_AD) * density_AC[age_i][delta] );
        //ret += log( (1.0 - q_AI) * (1.0 - q_AD) );
      }
      
      // open case
      if (s_ptr->final_outcome[i] == -1) {
        int delta = s_ptr->date_censor[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * q_AD * tail_AD[delta] +
                    (1.0 - q_AI) * (1.0 - q_AD) * tail_AC[age_i][delta] +
                    q_AI * tail_AI[delta] );
        //ret += log( (1.0 - q_AI) );
      }
      
    }  // end non-icu pathway
    
    // icu pathway
    if (s_ptr->icu[i] == 1) {
      
      // admission to icu
      {
        int delta = s_ptr->date_icu[i] - s_ptr->date_admission[i];
        ret += log( q_AI * density_AI[delta] );
        //ret += log( q_AI );
      }
      
      //if (i == 483) {
      //  print(1, ret);
      //}
      
      // non-stepdown pathway
      if (s_ptr->stepdown[i] == 0) {
        
        // icu to death
        if (s_ptr->final_outcome[i] == 1) {
          int delta = s_ptr->date_final_outcome[i] - s_ptr->date_icu[i];
          ret += log( q_AI * q_ID * density_ID[delta] );
          //ret += log( q_AI * q_ID );
        }
        
        // open case in icu
        if (s_ptr->final_outcome[i] == -1) {
          int delta = s_ptr->date_censor[i] - s_ptr->date_icu[i];
          ret += log( q_AI * q_ID * tail_ID[delta] +
                      q_AI * (1.0 - q_ID) * tail_IS[delta] );
          //ret += log( q_AI );
        }
        
      }  // end non-stepdown pathway
      
      // stepdown pathway
      if (s_ptr->stepdown[i] == 1) {
        
        // icu to stepdown
        {
          int delta = s_ptr->date_stepdown[i] - s_ptr->date_icu[i];
          ret += log( q_AI * (1.0 - q_ID) * density_IS[delta] );
          //ret += log( q_AI * (1.0 - q_ID) );
        }
        
        //if (i == 483) {
        //  print(2, ret);
        //}
        
        // stepdown to discharge
        if (s_ptr->final_outcome[i] == 2) {
          int delta = s_ptr->date_final_outcome[i] - s_ptr->date_stepdown[i];
          ret += log( q_AI * (1.0 - q_ID) * density_SC[delta] );
          //ret += log( q_AI * (1.0 - q_ID) );
          
          //if (i == 483) {
          //  print(3, ret, q_AI, q_ID, delta, density_SC[delta], m_SC, s_SC);
          //  print_vector(density_SC);
          //  Rcpp::stop("i exit");
          //}
          
        }
        
        // open case in stepdown
        if (s_ptr->final_outcome[i] == -1) {
          int delta = s_ptr->date_censor[i] - s_ptr->date_stepdown[i];
          ret += log( q_AI * (1.0 - q_ID) * tail_SC[delta] );
          //ret += log( q_AI * (1.0 - q_ID) );
        }
        
      }  // end stepdown pathway
      
    }  // end icu pathway
    
    //print(i, ret);
    //if(!std::isfinite(ret)) {
    //  Rcpp::stop("non finite");
    //}
    
  }  // end i loop
  }
  
  // ----------------------------------------------------------------
  // SitRep component of likelihood
  
  //print(ret);
  
  // loop through regions
  for (int region_i = 0; region_i < s_ptr->n_region; ++region_i) {
    //continue;
    
    // loop through sitrep age groups
    for (int age_i = 0; age_i < s_ptr->n_age_sitrep; ++age_i) {
      
      // sum log-likelihood over all data
      for (int i = 0; i < s_ptr->n_date_sitrep; ++i) {
        
        // adjust for spline offset
        int j = i - s_ptr->node_x[0];
        
        //double k = 100.0;
        
        // new admissions
        if (s_ptr->daily_influx[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->daily_influx[region_i][age_i][i], admission_incidence[region_i][age_i][j], true);
          //double p = k / (k + admission_incidence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->daily_influx[region_i][age_i][i], k, p, true);
        }
        
        // new deaths
        if (s_ptr->new_deaths[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->new_deaths[region_i][age_i][i], deaths_incidence[region_i][age_i][j], true);
          //double p = k / (k + deaths_incidence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->new_deaths[region_i][age_i][i], k, p, true);
        }
        
        // new discharges
        if (s_ptr->new_discharges[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->new_discharges[region_i][age_i][i], discharges_incidence[region_i][age_i][j], true);
          //double p = k / (k + discharges_incidence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->new_discharges[region_i][age_i][i], k, p, true);
        }
        
        // prevalence in general beds
        if (s_ptr->total_general[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->total_general[region_i][age_i][i], general_prevalence[region_i][age_i][j], true);
          //double p = k / (k + general_prevalence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->total_general[region_i][age_i][i], k, p, true);
        }
        
        // prevalence in critical beds
        if (s_ptr->total_critical[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->total_critical[region_i][age_i][i], critical_prevalence[region_i][age_i][j], true);
          //double p = k / (k + critical_prevalence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->total_critical[region_i][age_i][i], k, p, true);
        }
        
      }
      
    }  // end age_i loop
    
  }  // end region_i loop
  
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
void Particle::update_gamma_density(vector<double> &density_vec, double m, double s) {
  bool use_lookup = false;
  if (use_lookup) {
    
    double m_index = floor(m * 100);
    double s_index = floor(s * 100);
    if (m_index >= s_ptr->gamma_density_lookup.size()) {
      print(m, m_index);
      Rcpp::stop("m outside gamma_density_lookup range");
    }
    if (s_index >= s_ptr->gamma_density_lookup[m_index].size()) {
      print(s, s_index);
      Rcpp::stop("s outside gamma_density_lookup range");
    }
    density_vec = s_ptr->gamma_density_lookup[m_index][s_index];
    
  } else {
    
    for (unsigned int j = 0; j < density_vec.size(); ++j) {
      density_vec[j] = R::pgamma(j + 1, 1.0/(s*s), m*s*s, true, false) -
                       R::pgamma(j, 1.0/(s*s), m*s*s, true, false);
    }
    
  }
}

//------------------------------------------------
void Particle::update_gamma_tail(vector<double> &tail_vec, double m, double s) {
  bool use_lookup = false;
  if (use_lookup) {
    
    double m_index = floor(m * 100);
    double s_index = floor(s * 100);
    if (m_index >= s_ptr->gamma_tail_lookup.size()) {
      print(m, m_index);
      Rcpp::stop("m outside gamma_tail_lookup range");
    }
    if (s_index >= s_ptr->gamma_tail_lookup[m_index].size()) {
      print(s, s_index);
      Rcpp::stop("s outside gamma_tail_lookup range");
    }
    tail_vec = s_ptr->gamma_tail_lookup[m_index][s_index];
    
  } else {
    
    for (unsigned int j = 0; j < tail_vec.size(); ++j) {
      tail_vec[j] = R::pgamma(j + 1, 1.0/(s*s), m*s*s, false, false);
    }
    
  }
}
