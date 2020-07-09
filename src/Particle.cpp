
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
  //node_y = vector<vector<double>>(s_ptr->n_region, vector<double>(s_ptr->n_node));
  
  // rescaling parameters
  scale_p_AI = vector<double>(s_ptr->n_region);
  scale_p_AD = vector<double>(s_ptr->n_region);
  scale_p_ID = vector<double>(s_ptr->n_region);
  
  // lab test weights
  pos_on_day = vector<double>(s_ptr->n_date_sitrep);
  neg_by_day = vector<double>(s_ptr->n_date_sitrep);
  pos_by_day = vector<double>(s_ptr->n_date_sitrep);
  
  // vector over ages for cubic splines
  age_seq = vector<double>(s_ptr->max_indlevel_age + 1);
  for (unsigned int i = 0; i < age_seq.size(); ++i) {
    age_seq[i] = i;
  }
  
  // transition probabilities
  p_AI_node = vector<double>(s_ptr->p_AI_noden);
  p_AI = vector<double>(s_ptr->max_indlevel_age + 1);
  p_AD_node = vector<double>(s_ptr->p_AD_noden);
  p_AD = vector<double>(s_ptr->max_indlevel_age + 1);
  p_ID_node = vector<double>(s_ptr->p_ID_noden);
  p_ID = vector<double>(s_ptr->max_indlevel_age + 1);
  
  // mean durations
  m_AC_node = vector<double>(s_ptr->m_AC_noden);
  m_AC = vector<double>(s_ptr->max_indlevel_age + 1);
  
  // objects for storing progression
  admissions_spline = vector<vector<double>>(s_ptr->n_region, vector<double>(s_ptr->n_date_sitrep));
  admission_incidence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  deaths_incidence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  discharges_incidence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  general_prevalence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  critical_prevalence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  
  delta_stepup = std::vector<double>(s_ptr->n_date_sitrep);
  delta_stepdown = std::vector<double>(s_ptr->n_date_sitrep);
  delta_deaths_general = std::vector<double>(s_ptr->n_date_sitrep);
  delta_discharges_general = std::vector<double>(s_ptr->n_date_sitrep);
  delta_open_general = std::vector<double>(s_ptr->n_date_sitrep);
  delta_deaths_critical = std::vector<double>(s_ptr->n_date_sitrep);
  delta_open_critical = std::vector<double>(s_ptr->n_date_sitrep);
  pos_on_day = std::vector<double>(s_ptr->n_date_sitrep);
  pos_by_day = std::vector<double>(s_ptr->n_date_sitrep);
  neg_by_day = std::vector<double>(s_ptr->n_date_sitrep);
  
  // theta is the parameter vector in natural space
  theta = s_ptr->theta_init;
  theta_prop = vector<double>(d);
  
  // phi is a vector of transformed parameters
  phi = vector<double>(d);
  theta_to_phi();
  phi_prop = vector<double>(d);
  
  // run through likelihood calculation once for each parameter to initialise lookup tables
  for (int i = 0; i < d; ++i) {
    double dummy = get_loglike(theta, i, true);
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
      double dummy = get_loglike(theta, i, false);
      
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
  
  // transition cubic spline nodes
  int pi = 0;
  for (int i = 0; i < s_ptr->p_AI_noden; ++i) {
    p_AI_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->p_AD_noden; ++i) {
    p_AD_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->p_ID_noden; ++i) {
    p_ID_node[i] = theta[pi++];
  }
  
  // mean durations
  double m_AI = theta[pi++];
  double m_AD = theta[pi++];
  for (int i = 0; i < s_ptr->m_AC_noden; ++i) {
    m_AC_node[i] = theta[pi++];
  }
  //double m_AC = theta[pi++];
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
  
  // censoring parameters
  double c_AD = theta[pi++];
  double c_AC = theta[pi++];
  double c_ID = theta[pi++];
  double c_IS = theta[pi++];
  
  // admission to labtest
  double m_AL = theta[pi++];
  
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
  int update_region = s_ptr->update_region[theta_i];
  
  // initialise loglikelihood
  double ret = 0.0;
  
  // ----------------------------------------------------------------
  // calculate cubic plines
  
  // get cubic spline for p_AI and transform to [0,1] interval
  cubic_spline(s_ptr->p_AI_nodex, p_AI_node, age_seq, p_AI);
  for (unsigned int i = 0; i < p_AI.size(); ++i) {
    p_AI[i] = 1.0 / (1.0 + exp(-p_AI[i]));
  }
  
  // get cubic spline for p_AD and transform to [0,1] interval
  cubic_spline(s_ptr->p_AD_nodex, p_AD_node, age_seq, p_AD);
  for (unsigned int i = 0; i < p_AD.size(); ++i) {
    p_AD[i] = 1.0 / (1.0 + exp(-p_AD[i]));
  }
  
  // get cubic spline for p_ID and transform to [0,1] interval
  cubic_spline(s_ptr->p_ID_nodex, p_ID_node, age_seq, p_ID);
  for (unsigned int i = 0; i < p_ID.size(); ++i) {
    p_ID[i] = 1.0 / (1.0 + exp(-p_ID[i]));
  }
  
  // get cubic spline for m_AC and transform to [0,1] interval
  cubic_spline(s_ptr->m_AC_nodex, m_AC_node, age_seq, m_AC);
  for (unsigned int i = 0; i < m_AC.size(); ++i) {
    m_AC[i] = 20.0 / (1.0 + exp(-m_AC[i]));
  }
  
  // ----------------------------------------------------------------
  // update objects used when calculating sitrep likelihood
  
  // calculate lab test weights
  for (unsigned int i = 0; i < neg_by_day.size(); ++i) {
    pos_on_day[i] = 0;//get_delay_density(i, m_AL, 1.0);
    neg_by_day[i] = 0;//get_delay_tail(i, m_AL, 1.0);
    pos_by_day[i] = 1.0;// - neg_by_day[i];
  }
  pos_on_day[0] = 1.0;
  
  // loop through sitrep regions
  for (int region_i = 0; region_i < s_ptr->n_region; ++region_i) {
    //continue;
    
    // skip based on update rules
    if (update_region != (region_i + 1) && update_region != 0) {
      continue;
    }
    
    // loop through sitrep age groups
    for (int age_sitrep = 0; age_sitrep < s_ptr->n_age_sitrep; ++age_sitrep) {
      
      // reset progression objects
      fill(admission_incidence[region_i][age_sitrep].begin(), admission_incidence[region_i][age_sitrep].end(), 0.0);
      fill(deaths_incidence[region_i][age_sitrep].begin(), deaths_incidence[region_i][age_sitrep].end(), 0.0);
      fill(discharges_incidence[region_i][age_sitrep].begin(), discharges_incidence[region_i][age_sitrep].end(), 0.0);
      fill(general_prevalence[region_i][age_sitrep].begin(), general_prevalence[region_i][age_sitrep].end(), 0.0);
      fill(critical_prevalence[region_i][age_sitrep].begin(), critical_prevalence[region_i][age_sitrep].end(), 0.0);
      
      // loop through all integer ages inside this sitrep age group
      for (unsigned int age2 = 0; age2 < s_ptr->age_weights[age_sitrep].size(); ++age2) {
        
        // get integer age
        int age = s_ptr->age_values[age_sitrep][age2];
        double weight = s_ptr->age_weights[age_sitrep][age2];
        
        // get age-specific parameters
        double p_AI_j = p_AI[age] * scale_p_AI[region_i];
        double p_AD_j = p_AD[age] * scale_p_AD[region_i];
        double p_ID_j = p_ID[age];// * scale_p_ID[region_i];
        double m_AC_j = m_AC[age];
        
        // check scaling doesn't push probabilities above 1
        if (p_AI_j > 1.0 || p_AD_j > 1.0 || p_ID_j > 1.0) {
          return -DBL_MAX/100.0;
        }
        
        // calculate unscaled progression vectors
        for (int i = 0; i < s_ptr->n_date_sitrep; ++i) {
          delta_stepup[i] = p_AI_j * get_delay_density(i, m_AI, s_AI);
          delta_stepdown[i] = (1 - p_ID_j) * get_delay_density(i, m_IS, s_IS);
          
          delta_deaths_general[i] = (1 - p_AI_j) * p_AD_j * get_delay_density(i, m_AD, s_AD);
          delta_discharges_general[i] =  (1 - p_AI_j) * (1 - p_AD_j) * get_delay_density(i, m_AC_j, s_AC);
          delta_open_general[i] =  (1 - p_AI_j) * p_AD_j * get_delay_tail(i, m_AD, s_AD) +
                                   (1 - p_AI_j) * (1 - p_AD_j) * get_delay_tail(i, m_AC_j, s_AC) +
                                   p_AI_j * get_delay_tail(i, m_AI, s_AI);
          delta_deaths_critical[i] = p_ID_j * get_delay_density(i, m_ID, s_ID);
          delta_open_critical[i] = p_ID_j * get_delay_tail(i, m_ID, s_ID) +
                                   (1 - p_ID_j) * get_delay_tail(i, m_IS, s_IS);
        }
        
        //print_vector(delta_open_general);
        
        // store incidence in stepup and stepdown care
        vector<double> stepup(s_ptr->n_date_sitrep);
        vector<double> stepdown(s_ptr->n_date_sitrep);
        
        // loop through each spline day in turn
        for (int i = 0; i < s_ptr->n_date_sitrep; ++i) {
          
          // get true admissions on this day from spline
          //double true_admissions = weight * admissions_spline[region_i][i];
          double true_admissions = weight * s_ptr->daily_influx[region_i][age_sitrep][i];
          
          //print(i, weight, true_admissions);
          
          // project true admissions out over future days
          for (int j = i; j < s_ptr->n_date_sitrep; ++j) {
            
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
            double discharges_stepdown = stepdown[i] * get_delay_density(j - i, m_SC, s_SC);
            double open_stepdown = stepdown[i] * get_delay_tail(j - i, m_SC, s_SC);
            
            // delay from admission to testing
            double pos_by_day_j = pos_by_day[j - i];
            double pos_on_day_j = pos_on_day[j - i];
            
            //print(i, j, pos_on_day_j, open_general, open_stepdown, open_critical);
            
            // incidence of observed admission
            //admission_incidence[region_i][age_sitrep][j] += pos_on_day_j * (open_general + open_critical + open_stepdown);
            admission_incidence[region_i][age_sitrep][j] += pos_on_day_j * true_admissions;
            
            // incidence of death and discharge
            deaths_incidence[region_i][age_sitrep][j] += pos_by_day_j * (deaths_general + deaths_critical);
            discharges_incidence[region_i][age_sitrep][j] += pos_by_day_j * (discharges_general + discharges_stepdown);
            
            // prevalence in general and critical beds
            general_prevalence[region_i][age_sitrep][j] += pos_by_day_j * (open_general + open_stepdown);
            critical_prevalence[region_i][age_sitrep][j] += pos_by_day_j * open_critical;
            
          }  // end j loop through spline days
          
        }  // end i loop through spline days
        
      }  // end age2 loop
      
    }  // end age_sitrep loop
    
  }  // end region_i loop
  
  
  // option to return at this stage
  if (quick_exit) {
    return 0.0;
  }
  
  // ----------------------------------------------------------------
  // individual-level component of likelihood
  
  // sum log-likelihood over individual-level data
  if (true) {  // (optionally skip this section without ever going into loop)
  for (int i = 0; i < s_ptr->n_ind; ++i) {
    
    // get age group of this individual
    int age_i = s_ptr->age[i];
    
    // get age-specific parameters
    double q_AI = p_AI[age_i];
    double q_AD = p_AD[age_i];
    double q_ID = p_ID[age_i];
    
    // if no recorded ICU stay
    if (s_ptr->icu[i] == 0) {
      
      // admission to death
      if (s_ptr->final_outcome[i] == 1) {
        int delta = s_ptr->date_final_outcome[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * q_AD * (1 - c_AD) * get_delay_density(delta, m_AD, s_AD) );
      }
      
      // admission to discharge
      if (s_ptr->final_outcome[i] == 2) {
        int delta = s_ptr->date_final_outcome[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * (1.0 - q_AD) * (1 - c_AC) * get_delay_density(delta, m_AC[age_i], s_AC) );
      }
      
      // open case in general ward
      if (s_ptr->final_outcome[i] == -1) {
        int delta = s_ptr->date_censor[i] - s_ptr->date_admission[i];
        ret += log( (1.0 - q_AI) * q_AD * (c_AD + (1 - c_AD) * get_delay_tail(delta, m_AD, s_AD)) +
                    (1.0 - q_AI) * (1.0 - q_AD) * (c_AC + (1 - c_AC) * get_delay_tail(delta, m_AC[age_i], s_AC)) +
                    q_AI * get_delay_tail(delta, m_AI, s_AI) );
      }
      
    }
    
    // if recorded ICU stay
    if (s_ptr->icu[i] == 1) {
      //continue;
      
      // admission to icu
      {
        int delta = s_ptr->date_icu[i] - s_ptr->date_admission[i];
        ret += log( q_AI * get_delay_density(delta, m_AI, s_AI) );
      }
      
      //continue;
      
      // if no recorded stepdown stay
      if (s_ptr->stepdown[i] == 0) {
        
        // icu to death
        if (s_ptr->final_outcome[i] == 1) {
          int delta = s_ptr->date_final_outcome[i] - s_ptr->date_icu[i];
          ret += log( q_ID * (1 - c_ID) * get_delay_density(delta, m_ID, s_ID) );
        }
        
        // open case in icu
        if (s_ptr->final_outcome[i] == -1) {
          int delta = s_ptr->date_censor[i] - s_ptr->date_icu[i];
          ret += log( q_ID * (c_ID + (1 - c_ID) * get_delay_tail(delta, m_ID, s_ID)) +
                      (1.0 - q_ID) * (c_IS + (1 - c_IS) * get_delay_tail(delta, m_IS, s_IS)) );
        }
        
      }
      
      // if recorded stepdown stay
      if (s_ptr->stepdown[i] == 1) {
        
        // icu to stepdown
        {
          int delta = s_ptr->date_stepdown[i] - s_ptr->date_icu[i];
          ret += log( (1.0 - q_ID) * (1 - c_IS) * get_delay_density(delta, m_IS, s_IS) );
        }
        
        // stepdown to discharge
        if (s_ptr->final_outcome[i] == 2) {
          int delta = s_ptr->date_final_outcome[i] - s_ptr->date_stepdown[i];
          ret += log( get_delay_density(delta, m_SC, s_SC) );
        }
        
        // open case in stepdown
        if (s_ptr->final_outcome[i] == -1) {
          int delta = s_ptr->date_censor[i] - s_ptr->date_stepdown[i];
          ret += log( get_delay_tail(delta, m_SC, s_SC) );
        }
        
      }  // end stepdown pathway
      
    }  // end icu pathway
    
  }  // end i loop
  }
  
  // ----------------------------------------------------------------
  // SitRep component of likelihood
  
  //print(ret);
  
  // loop through regions
  for (int region_i = 0; region_i < s_ptr->n_region; ++region_i) {
    //continue;
    
    double eps = 0.001;
    
    // loop through sitrep age groups
    for (int age_i = 0; age_i < s_ptr->n_age_sitrep; ++age_i) {
      
      // sum log-likelihood over all data
      for (int i = 0; i < s_ptr->n_date_sitrep; ++i) {
        
        //print(age_i, i);
        
        // adjust for spline offset
        int j = i;// + 15;
        
        //double k = 100.0;
        
        // new admissions
        if (s_ptr->daily_influx[region_i][age_i][i] != -1) {
          //ret += R::dpois(s_ptr->daily_influx[region_i][age_i][i], admission_incidence[region_i][age_i][j], true);
          //double p = k / (k + admission_incidence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->daily_influx[region_i][age_i][i], k, p, true);
        }
        
        // new deaths
        if (s_ptr->new_deaths[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->new_deaths[region_i][age_i][i], deaths_incidence[region_i][age_i][j] + eps, true);
          //double p = k / (k + deaths_incidence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->new_deaths[region_i][age_i][i], k, p, true);
        }
        //print("foo1", ret);
        
        // new discharges
        if (s_ptr->new_discharges[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->new_discharges[region_i][age_i][i], discharges_incidence[region_i][age_i][j] + eps, true);
          //double p = k / (k + discharges_incidence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->new_discharges[region_i][age_i][i], k, p, true);
        }
        //print("foo2", ret);
        
        // prevalence in general beds
        if (s_ptr->total_general[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->total_general[region_i][age_i][i], general_prevalence[region_i][age_i][j] + eps, true);
          //double p = k / (k + general_prevalence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->total_general[region_i][age_i][i], k, p, true);
        }
        //print("foo3", ret);
        //if (!isfinite(ret)) {
        //  print(s_ptr->total_general[region_i][age_i][i], general_prevalence[region_i][age_i][j]);
        //  Rcpp::stop("bar3");
        //}
        
        // prevalence in critical beds
        if (s_ptr->total_critical[region_i][age_i][i] != -1) {
          ret += R::dpois(s_ptr->total_critical[region_i][age_i][i], critical_prevalence[region_i][age_i][j] + eps, true);
          //double p = k / (k + critical_prevalence[region_i][age_i][j]);
          //ret += R::dnbinom(s_ptr->total_critical[region_i][age_i][i], k, p, true);
        }
        //print("foo4", ret);
        
      }
      
    }  // end age_i loop
    
  }  // end region_i loop
  
  // ----------------------------------------------------------------
  // return
  
  //foobar();
  //print(ret);
  
  // catch underflow
  if (!std::isfinite(ret)) {
    //print("catch");
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
// get density of delay distribution on day x
double Particle::get_delay_density(int x, double m, double s) {
  double ret;
  bool use_lookup = false;
  if (use_lookup) {
    double m_index = floor(m * 100);
    double s_index = floor(s * 100);
    ret = s_ptr->gamma_density_lookup[m_index][s_index][x];
  } else {
    ret = R::pgamma(x + 1, 1.0/(s*s), m*s*s, true, false) -
          R::pgamma(x, 1.0/(s*s), m*s*s, true, false) + 0.0001;
  }
  if (ret > 1.0) {
    ret = 1.0;
  }
  return ret;
}

//------------------------------------------------
// get tail of delay distribution past day x
double Particle::get_delay_tail(int x, double m, double s) {
  double ret;
  bool use_lookup = false;
  if (use_lookup) {
    double m_index = floor(m * 100);
    double s_index = floor(s * 100);
    ret = s_ptr->gamma_tail_lookup[m_index][s_index][x];
  } else {
    ret = R::pgamma(x + 1, 1.0/(s*s), m*s*s, false, false);
  }
  if (ret > 1.0) {
    ret = 1.0;
  }
  return ret;
}

