
#include "Particle.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

//------------------------------------------------
// initialise/reset particle
void Particle::init(System &s) {
  
  // pointer to system object
  this->s_ptr = &s;
  
  // parameters
  d = s_ptr->d;
  
  // rescaling parameters
  scale_p_AI = vector<double>(s_ptr->n_region);
  scale_p_AD = vector<double>(s_ptr->n_region);
  scale_p_ID = vector<double>(s_ptr->n_region);
  
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
  m_AI_node = vector<double>(s_ptr->m_AI_noden);
  m_AI = vector<double>(s_ptr->max_indlevel_age + 1);
  m_AD_node = vector<double>(s_ptr->m_AD_noden);
  m_AD = vector<double>(s_ptr->max_indlevel_age + 1);
  m_AC_node = vector<double>(s_ptr->m_AC_noden);
  m_AC = vector<double>(s_ptr->max_indlevel_age + 1);
  m_ID_node = vector<double>(s_ptr->m_ID_noden);
  m_ID = vector<double>(s_ptr->max_indlevel_age + 1);
  m_IS_node = vector<double>(s_ptr->m_IS_noden);
  m_IS = vector<double>(s_ptr->max_indlevel_age + 1);
  m_SC_node = vector<double>(s_ptr->m_SC_noden);
  m_SC = vector<double>(s_ptr->max_indlevel_age + 1);
  
  // objects for storing progression
  deaths_incidence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  discharges_incidence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  general_prevalence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  critical_prevalence = vector<vector<vector<double>>>(s_ptr->n_region, vector<vector<double>>(s_ptr->n_age_sitrep, vector<double>(s_ptr->n_date_sitrep)));
  
  delta_stepup = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  delta_stepdown = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  delta_deaths_general = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  delta_discharges_general = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  delta_open_general = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  delta_deaths_critical = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  delta_open_critical = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  
  stepup = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  stepdown = std::vector<std::vector<double>>(s_ptr->n_threads, std::vector<double>(s_ptr->n_date_sitrep));
  
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
    logprior_prop = get_logprior(theta_prop, i);
    
    // calculate Metropolis-Hastings ratio
    double MH = beta*(loglike_prop - loglike) + (logprior_prop - logprior) + adj;
    
    // accept or reject move
    bool MH_accept = (log(runif_0_1()) < MH);
    
    //if (i == 0) {
    //  print(i, theta[0], theta_prop[0], loglike, loglike_prop, logprior, logprior_prop, MH, MH_accept);
    //}
    
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
  
  //return(0.0);
  
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
  for (int i = 0; i < s_ptr->m_AI_noden; ++i) {
    m_AI_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_AD_noden; ++i) {
    m_AD_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_AC_noden; ++i) {
    m_AC_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_ID_noden; ++i) {
    m_ID_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_IS_noden; ++i) {
    m_IS_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_SC_noden; ++i) {
    m_SC_node[i] = theta[pi++];
  }
  
  // coefficients of variation of durations
  s_AI = theta[pi++];
  s_AD = theta[pi++];
  s_AC = theta[pi++];
  s_ID = theta[pi++];
  s_IS = theta[pi++];
  s_SC = theta[pi++];
  
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
  
  // get cubic spline for m_AI and transform to [0,20] interval
  cubic_spline(s_ptr->m_AI_nodex, m_AI_node, age_seq, m_AI);
  for (unsigned int i = 0; i < m_AI.size(); ++i) {
    m_AI[i] = 20.0 / (1.0 + exp(-m_AI[i]));
  }
  
  // get cubic spline for m_AD and transform to [0,20] interval
  cubic_spline(s_ptr->m_AD_nodex, m_AD_node, age_seq, m_AD);
  for (unsigned int i = 0; i < m_AD.size(); ++i) {
    m_AD[i] = 20.0 / (1.0 + exp(-m_AD[i]));
  }
  
  // get cubic spline for m_AC and transform to [0,20] interval
  cubic_spline(s_ptr->m_AC_nodex, m_AC_node, age_seq, m_AC);
  for (unsigned int i = 0; i < m_AC.size(); ++i) {
    m_AC[i] = 20.0 / (1.0 + exp(-m_AC[i]));
  }
  
  // get cubic spline for m_ID and transform to [0,20] interval
  cubic_spline(s_ptr->m_ID_nodex, m_ID_node, age_seq, m_ID);
  for (unsigned int i = 0; i < m_ID.size(); ++i) {
    m_ID[i] = 20.0 / (1.0 + exp(-m_ID[i]));
  }
  
  // get cubic spline for m_IS and transform to [0,20] interval
  cubic_spline(s_ptr->m_IS_nodex, m_IS_node, age_seq, m_IS);
  for (unsigned int i = 0; i < m_IS.size(); ++i) {
    m_IS[i] = 20.0 / (1.0 + exp(-m_IS[i]));
  }
  
  // get cubic spline for m_SC and transform to [0,20] interval
  cubic_spline(s_ptr->m_SC_nodex, m_SC_node, age_seq, m_SC);
  for (unsigned int i = 0; i < m_SC.size(); ++i) {
    m_SC[i] = 20.0 / (1.0 + exp(-m_SC[i]));
  }
  
  
  // ----------------------------------------------------------------
  // update objects used when calculating sitrep likelihood
  
  if (s_ptr->sitrep_loglike) {
    
    // loop through sitrep regions
    const size_t n_threads = s_ptr->n_threads;
  #pragma omp parallel for schedule(static) num_threads(n_threads)
    for (int region_i = 0; region_i < s_ptr->n_region; ++region_i) {
      update_region(region_i);
    }
    
  }
  
  // option to return at this stage
  if (quick_exit) {
    return 0.0;
  }
  
  
  // ----------------------------------------------------------------
  // individual-level component of likelihood
  
  // sum log-likelihood over individual-level data
  int n_age_indlevel = s_ptr->p_AI_numer.size();
  for (int i = 0; i < n_age_indlevel; ++i) {
    
    // transition probabilities
    ret += R::dbinom(s_ptr->p_AI_numer[i], s_ptr->p_AI_denom[i], p_AI[i], true);
    ret += R::dbinom(s_ptr->p_AD_numer[i], s_ptr->p_AD_denom[i], p_AD[i], true);
    ret += R::dbinom(s_ptr->p_ID_numer[i], s_ptr->p_ID_denom[i], p_ID[i], true);
    
    // duration m_AI
    for (unsigned int j = 0; j < s_ptr->m_AI_count[i].size(); ++j) {
      int tmp = s_ptr->m_AI_count[i][j];
      if (tmp > 0) {
        ret += tmp * log(get_delay_density(j, m_AI[i], s_AI));
      }
    }
    
    // duration m_AD
    for (unsigned int j = 0; j < s_ptr->m_AD_count[i].size(); ++j) {
      int tmp = s_ptr->m_AD_count[i][j];
      if (tmp > 0) {
        ret += tmp * log(get_delay_density(j, m_AD[i], s_AD));
      }
    }
    
    // duration m_AC
    for (unsigned int j = 0; j < s_ptr->m_AC_count[i].size(); ++j) {
      int tmp = s_ptr->m_AC_count[i][j];
      if (tmp > 0) {
        ret += tmp * log(get_delay_density(j, m_AC[i], s_AC));
      }
    }
    
    // duration m_ID
    for (unsigned int j = 0; j < s_ptr->m_ID_count[i].size(); ++j) {
      int tmp = s_ptr->m_ID_count[i][j];
      if (tmp > 0) {
        ret += tmp * log(get_delay_density(j, m_ID[i], s_ID));
      }
    }
    
    // duration m_IS
    for (unsigned int j = 0; j < s_ptr->m_IS_count[i].size(); ++j) {
      int tmp = s_ptr->m_IS_count[i][j];
      if (tmp > 0) {
        ret += tmp * log(get_delay_density(j, m_IS[i], s_IS));
      }
    }
    
    // duration m_SC
    for (unsigned int j = 0; j < s_ptr->m_SC_count[i].size(); ++j) {
      int tmp = s_ptr->m_SC_count[i][j];
      if (tmp > 0) {
        ret += tmp * log(get_delay_density(j, m_SC[i], s_SC));
      }
    }
    
    if (!isfinite(ret)) {
      //print(i, s_ptr->p_AI_denom[i], s_ptr->p_AI_numer[i], p_AI[i], ret);
      Rcpp::stop("ret non finite");
    }
  }
  
  // ----------------------------------------------------------------
  // SitRep component of likelihood
  
  if (s_ptr->sitrep_loglike) {
    
    // loop through regions
    for (int region_i = 0; region_i < s_ptr->n_region; ++region_i) {
      
      // loop through sitrep age groups
      for (int age_i = 0; age_i < s_ptr->n_age_sitrep; ++age_i) {
        
        // sum log-likelihood over all data
        for (int i = 0; i < s_ptr->n_date_sitrep; ++i) {
          
          // new deaths
          if (s_ptr->new_deaths[region_i][age_i][i] != -1) {
            double lambda = deaths_incidence[region_i][age_i][i];
            if (lambda < 0) {
              lambda = 0.001;
            }
            ret += R::dpois(s_ptr->new_deaths[region_i][age_i][i], lambda, true);
          }
          
          // new discharges
          if (s_ptr->new_discharges[region_i][age_i][i] != -1) {
            double lambda = discharges_incidence[region_i][age_i][i];
            if (lambda < 0) {
              lambda = 0.001;
            }
            ret += R::dpois(s_ptr->new_discharges[region_i][age_i][i], lambda, true);
          }
          
          // prevalence in general beds
          if (s_ptr->total_general[region_i][age_i][i] != -1) {
            double lambda = general_prevalence[region_i][age_i][i];
            if (lambda < 0) {
              lambda = 0.001;
            }
            ret += R::dpois(s_ptr->total_general[region_i][age_i][i], lambda, true);
          }
          
          // prevalence in critical beds
          if (s_ptr->total_critical[region_i][age_i][i] != -1) {
            double lambda = critical_prevalence[region_i][age_i][i];
            if (lambda < 0) {
              lambda = 0.001;
            }
            ret += R::dpois(s_ptr->total_critical[region_i][age_i][i], lambda, true);
          }
          
        }
        
      }  // end age_i loop
      
    }  // end region_i loop
  }
  
  
  // ----------------------------------------------------------------
  // return
  
  // catch underflow
  if (!std::isfinite(ret)) {
    ret = -DBL_MAX/100.0;
  }
  
  return ret;
}

//------------------------------------------------
// update objects required to calculate likelihood for a specific region
void Particle::update_region(int region_i) {
  
  //get thread index
  int thread_idx = omp_get_thread_num();
  //print(thread_idx);
  
  // loop through sitrep age groups
  for (int age_sitrep = 0; age_sitrep < s_ptr->n_age_sitrep; ++age_sitrep) {
    
    // reset progression objects
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
      double p_ID_j = p_ID[age] * scale_p_ID[region_i];
      double m_AI_j = m_AI[age];
      double m_AD_j = m_AD[age];
      double m_AC_j = m_AC[age];
      double m_ID_j = m_ID[age];
      double m_IS_j = m_IS[age];
      double m_SC_j = m_SC[age];
      
      // check scaling doesn't push probabilities above 1
      // to be compatible with openmp, this needs to be done outside of this loop.
      // if (p_AI_j > 1.0 || p_AD_j > 1.0 || p_ID_j > 1.0) {
      //   return -DBL_MAX/100.0;
      // }
      
      // calculate unscaled progression vectors
      for (int i = 0; i < s_ptr->n_date_sitrep; ++i) {
        delta_stepup[thread_idx][i] = p_AI_j * get_delay_density(i, m_AI_j, s_AI);
        delta_stepdown[thread_idx][i] = (1 - p_ID_j) * get_delay_density(i, m_IS_j, s_IS);
        
        delta_deaths_general[thread_idx][i] = (1 - p_AI_j) * p_AD_j * get_delay_density(i, m_AD_j, s_AD);
        delta_discharges_general[thread_idx][i] =  (1 - p_AI_j) * (1 - p_AD_j) * get_delay_density(i, m_AC_j, s_AC);
        delta_open_general[thread_idx][i] =  (1 - p_AI_j) * p_AD_j * get_delay_tail(i, m_AD_j, s_AD) +
          (1 - p_AI_j) * (1 - p_AD_j) * get_delay_tail(i, m_AC_j, s_AC) +
          p_AI_j * get_delay_tail(i, m_AI_j, s_AI);
        delta_deaths_critical[thread_idx][i] = p_ID_j * get_delay_density(i, m_ID_j, s_ID);
        delta_open_critical[thread_idx][i] = p_ID_j * get_delay_tail(i, m_ID_j, s_ID) +
          (1 - p_ID_j) * get_delay_tail(i, m_IS_j, s_IS);
      }
      
      fill(stepup[thread_idx].begin(), stepup[thread_idx].end(), 0.0);
      fill(stepdown[thread_idx].begin(), stepdown[thread_idx].end(), 0.0);
      
      // loop through each spline day in turn
      for (int i = 0; i < s_ptr->n_date_sitrep; ++i) {
        
        // get true admissions on this day from spline
        double true_admissions = weight * s_ptr->daily_influx[region_i][age_sitrep][i];
        
        // project true admissions out over future days
        for (int j = i; j < s_ptr->n_date_sitrep; ++j) {
          
          // update incidence in stepup and stepdown
          stepup[thread_idx][j] += true_admissions * delta_stepup[thread_idx][j - i];
          stepdown[thread_idx][j] += stepup[thread_idx][i] * delta_stepdown[thread_idx][j - i];
          
          // deaths, discharges and open cases in general ward
          double deaths_general = true_admissions * delta_deaths_general[thread_idx][j - i];
          double discharges_general = true_admissions * delta_discharges_general[thread_idx][j - i];
          double open_general = true_admissions * delta_open_general[thread_idx][j - i];
          
          // deaths and open cases in critical care
          double deaths_critical = stepup[thread_idx][i] * delta_deaths_critical[thread_idx][j - i];
          double open_critical = stepup[thread_idx][i] * delta_open_critical[thread_idx][j - i];
          
          // discharges and open cases in stepdown
          double discharges_stepdown = stepdown[thread_idx][i] * get_delay_density(j - i, m_SC_j, s_SC);
          double open_stepdown = stepdown[thread_idx][i] * get_delay_tail(j - i, m_SC_j, s_SC);
          
          // incidence of death and discharge
          deaths_incidence[region_i][age_sitrep][j] += deaths_general + deaths_critical;
          discharges_incidence[region_i][age_sitrep][j] += discharges_general + discharges_stepdown;
          
          // prevalence in general and critical beds
          general_prevalence[region_i][age_sitrep][j] += open_general + open_stepdown;
          critical_prevalence[region_i][age_sitrep][j] += open_critical;
          
        }  // end j loop through spline days
        
      }  // end i loop through spline days
      
    }  // end age2 loop
    
  }  // end age_sitrep loop
  
}

//------------------------------------------------
// define cpp logprior function
double Particle::get_logprior(vector<double> &theta, int theta_i) {
  
  // ----------------------------------------------------------------
  // unpack parameters
  
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
  for (int i = 0; i < s_ptr->m_AI_noden; ++i) {
    m_AI_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_AD_noden; ++i) {
    m_AD_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_AC_noden; ++i) {
    m_AC_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_ID_noden; ++i) {
    m_ID_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_IS_noden; ++i) {
    m_IS_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->m_SC_noden; ++i) {
    m_SC_node[i] = theta[pi++];
  }
  
  // ----------------------------------------------------------------
  // apply transformations and priors
  
  double k = 0.5;  // smoothing parameter
  double ret = 0.0;
  
  for (int i = 0; i < s_ptr->p_AI_noden; ++i) {
    if (i == 0) {
      ret += -p_AI_node[i] -2*log(1 + exp(-p_AI_node[i]));
    } else {
      ret += R::dnorm(p_AI_node[i], p_AI_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->p_AD_noden; ++i) {
    if (i == 0) {
      ret += -p_AD_node[i] -2*log(1 + exp(-p_AD_node[i]));
    } else {
      ret += R::dnorm(p_AD_node[i], p_AD_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->p_ID_noden; ++i) {
    if (i == 0) {
      ret += -p_ID_node[i] -2*log(1 + exp(-p_ID_node[i]));
    } else {
      ret += R::dnorm(p_ID_node[i], p_ID_node[i-1], k, true);
    }
  }
  
  // mean durations
  for (int i = 0; i < s_ptr->m_AI_noden; ++i) {
    if (i == 0) {
      ret += -m_AI_node[i] -2*log(1 + exp(-m_AI_node[i]));
    } else {
      ret += R::dnorm(m_AI_node[i], m_AI_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->m_AD_noden; ++i) {
    if (i == 0) {
      ret += -m_AD_node[i] -2*log(1 + exp(-m_AD_node[i]));
    } else {
      ret += R::dnorm(m_AD_node[i], m_AD_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->m_AC_noden; ++i) {
    if (i == 0) {
      ret += -m_AC_node[i] -2*log(1 + exp(-m_AC_node[i]));
    } else {
      ret += R::dnorm(m_AC_node[i], m_AC_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->m_ID_noden; ++i) {
    if (i == 0) {
      ret += -m_ID_node[i] -2*log(1 + exp(-m_ID_node[i]));
    } else {
      ret += R::dnorm(m_ID_node[i], m_ID_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->m_IS_noden; ++i) {
    if (i == 0) {
      ret += -m_IS_node[i] -2*log(1 + exp(-m_IS_node[i]));
    } else {
      ret += R::dnorm(m_IS_node[i], m_IS_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->m_SC_noden; ++i) {
    if (i == 0) {
      ret += -m_SC_node[i] -2*log(1 + exp(-m_SC_node[i]));
    } else {
      ret += R::dnorm(m_SC_node[i], m_SC_node[i-1], k, true);
    }
  }
  
  return ret;
}

//------------------------------------------------
// get density of delay distribution on day x
double Particle::get_delay_density(int x, double m, double s) {
#define USE_LOOKUP
#ifdef USE_LOOKUP
  double m_index = floor(m * 100);
  double s_index = floor(s * 100);
  if ((m_index < 0) || (m_index > 2000) || (s_index < 0) || (s_index > 100) || (x < 0)) {
    print("get_delay_density outside lookup range");
    print(x, m, s, m_index, s_index);
    Rcpp::stop("");
  }
  if (x > 100) {
    return 1e-200;
  }
  return s_ptr->gamma_density_lookup[m_index][s_index][x];
#else
  return 1.0;
#endif
}

//------------------------------------------------
// get tail of delay distribution past day x
double Particle::get_delay_tail(int x, double m, double s) {
#define USE_LOOKUP
#ifdef USE_LOOKUP
  double m_index = floor(m * 100);
  double s_index = floor(s * 100);
  if ((m_index < 0) || (m_index > 2000) || (s_index < 0) || (s_index > 100) || (x < 0)) {
    print("get_delay_tail outside lookup range");
    print(x, m, s, m_index, s_index);
    Rcpp::stop("");
  }
  if (x > 100) {
    return 1e-200;
  }
  return s_ptr->gamma_tail_lookup[m_index][s_index][x];
#else
  return 1.0;
#endif
}

