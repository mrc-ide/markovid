
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
  
  // vector over ages for cubic splines
  age_seq = vector<double>(s_ptr->max_indlevel_age + 1);
  for (unsigned int i = 0; i < age_seq.size(); ++i) {
    age_seq[i] = i;
  }
  
  // transition probabilities
  p_AI_node = vector<double>(s_ptr->n_node);
  p_AI = vector<double>(s_ptr->max_indlevel_age + 1);
  p_AD_node = vector<double>(s_ptr->n_node);
  p_AD = vector<double>(s_ptr->max_indlevel_age + 1);
  p_ID_node = vector<double>(s_ptr->n_node);
  p_ID = vector<double>(s_ptr->max_indlevel_age + 1);
  p_SD_node = vector<double>(s_ptr->n_node);
  p_SD = vector<double>(s_ptr->max_indlevel_age + 1);
  
  // theta is the parameter vector in natural space
  theta = s_ptr->theta_init;
  theta_prop = vector<double>(d);
  
  // phi is a vector of transformed parameters
  phi = vector<double>(d);
  theta_to_phi();
  phi_prop = vector<double>(d);
  
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
    loglike_prop = get_loglike(theta_prop, i);
    logprior_prop = get_logprior(theta_prop, i);
    
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
      
      // Robbins-Monro negative update (on the log scale)
      bw[i] = exp(log(bw[i]) - bw_stepsize*0.234/sqrt(bw_index[i]));
      bw_index[i]++;
      
    } // end MH step
    
  }  // end loop over parameters
    
}  // end update_univar function

//------------------------------------------------
// define cpp loglike function
double Particle::get_loglike(vector<double> &theta, int theta_i) {
  
  //return(0.0);
  
  // ----------------------------------------------------------------
  // unpack parameters and define fixed/derived parameters
  
  // transition cubic spline nodes
  int pi = 0;
  for (int i = 0; i < s_ptr->n_node; ++i) {
    p_AI_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    p_AD_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    p_ID_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    p_SD_node[i] = theta[pi++];
  }
  
  // mean durations
  m_AI = theta[pi++];
  m_AD = theta[pi++];
  m_AC = theta[pi++];
  m_ID = theta[pi++];
  m_I1S = theta[pi++];
  m_I2S = theta[pi++];
  m_SD = theta[pi++];
  m_SC = theta[pi++];
  
  // Erlang shape of durations
  s_AI = theta[pi++];
  s_AD = theta[pi++];
  s_AC = theta[pi++];
  s_ID = theta[pi++];
  s_I1S = theta[pi++];
  s_I2S = theta[pi++];
  s_SD = theta[pi++];
  s_SC = theta[pi++];
  
  // initialise loglikelihood
  double ret = 0.0;
  
  
  // ----------------------------------------------------------------
  // calculate cubic plines
  
  // get cubic spline for p_AI and transform to [0,1] interval
  cubic_spline(s_ptr->node_x, p_AI_node, age_seq, p_AI);
  for (unsigned int i = 0; i < p_AI.size(); ++i) {
    p_AI[i] = 1.0 / (1.0 + exp(-p_AI[i]));
  }
  
  // get cubic spline for p_AD and transform to [0,1] interval
  cubic_spline(s_ptr->node_x, p_AD_node, age_seq, p_AD);
  for (unsigned int i = 0; i < p_AD.size(); ++i) {
    p_AD[i] = 1.0 / (1.0 + exp(-p_AD[i]));
  }
  
  // get cubic spline for p_ID and transform to [0,1] interval
  cubic_spline(s_ptr->node_x, p_ID_node, age_seq, p_ID);
  for (unsigned int i = 0; i < p_ID.size(); ++i) {
    p_ID[i] = 1.0 / (1.0 + exp(-p_ID[i]));
  }
  
  // get cubic spline for p_SD and transform to [0,1] interval
  cubic_spline(s_ptr->node_x, p_SD_node, age_seq, p_SD);
  for (unsigned int i = 0; i < p_SD.size(); ++i) {
    p_SD[i] = 1.0 / (1.0 + exp(-p_SD[i]));
  }
  
  
  // ----------------------------------------------------------------
  // individual-level component of likelihood
  
  // transition probabilities
  for (int i = 0; i < (s_ptr->max_indlevel_age + 1); ++i) {
    ret += R::dbinom(s_ptr->p_AI_numer[i], s_ptr->p_AI_denom[i], p_AI[i], true);
    ret += R::dbinom(s_ptr->p_AD_numer[i], s_ptr->p_AD_denom[i], p_AD[i], true);
    ret += R::dbinom(s_ptr->p_ID_numer[i], s_ptr->p_ID_denom[i], p_ID[i], true);
    ret += R::dbinom(s_ptr->p_SD_numer[i], s_ptr->p_SD_denom[i], p_SD[i], true);
  }
  
  // duration m_AI
  for (unsigned int j = 0; j < s_ptr->m_AI_count.size(); ++j) {
    int tmp = s_ptr->m_AI_count[j];
    if (tmp > 0) {
      ret += tmp * log(get_delay_density(j, m_AI, s_AI));
    }
  }
  
  // duration m_AD
  for (unsigned int j = 0; j < s_ptr->m_AD_count.size(); ++j) {
    int tmp = s_ptr->m_AD_count[j];
    if (tmp > 0) {
      ret += tmp * log(get_delay_density(j, m_AD, s_AD));
    }
  }
  
  // duration m_AC
  for (unsigned int j = 0; j < s_ptr->m_AC_count.size(); ++j) {
    int tmp = s_ptr->m_AC_count[j];
    if (tmp > 0) {
      ret += tmp * log(get_delay_density(j, m_AC, s_AC));
    }
  }
  
  // duration m_ID
  for (unsigned int j = 0; j < s_ptr->m_ID_count.size(); ++j) {
    int tmp = s_ptr->m_ID_count[j];
    if (tmp > 0) {
      ret += tmp * log(get_delay_density(j, m_ID, s_ID));
    }
  }
  
  // duration m_I1S
  for (unsigned int j = 0; j < s_ptr->m_I1S_count.size(); ++j) {
    int tmp = s_ptr->m_I1S_count[j];
    if (tmp > 0) {
      ret += tmp * log(get_delay_density(j, m_I1S, s_I1S));
    }
  }
  
  // duration m_I2S
  for (unsigned int j = 0; j < s_ptr->m_I2S_count.size(); ++j) {
    int tmp = s_ptr->m_I2S_count[j];
    if (tmp > 0) {
      ret += tmp * log(get_delay_density(j, m_I2S, s_I2S));
    }
  }
  
  // duration m_SD
  for (unsigned int j = 0; j < s_ptr->m_SD_count.size(); ++j) {
    int tmp = s_ptr->m_SD_count[j];
    if (tmp > 0) {
      ret += tmp * log(get_delay_density(j, m_SD, s_SD));
    }
  }
  
  // duration m_SC
  for (unsigned int j = 0; j < s_ptr->m_SC_count.size(); ++j) {
    int tmp = s_ptr->m_SC_count[j];
    if (tmp > 0) {
      ret += tmp * log(get_delay_density(j, m_SC, s_SC));
    }
  }
  
  if (!isfinite(ret)) {
    //print(i, s_ptr->p_AI_denom[i], s_ptr->p_AI_numer[i], p_AI[i], ret);
    Rcpp::stop("ret non finite");
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
// define cpp logprior function
double Particle::get_logprior(vector<double> &theta, int theta_i) {
  
  // ----------------------------------------------------------------
  // unpack parameters
  
  // transition cubic spline nodes
  int pi = 0;
  for (int i = 0; i < s_ptr->n_node; ++i) {
    p_AI_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    p_AD_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    p_ID_node[i] = theta[pi++];
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    p_SD_node[i] = theta[pi++];
  }
  
  
  // ----------------------------------------------------------------
  // apply transformations and priors
  
  double k = 0.5;  // smoothing parameter
  double ret = 0.0;
  
  for (int i = 0; i < s_ptr->n_node; ++i) {
    if (i == 0) {
      ret += -p_AI_node[i] -2*log(1 + exp(-p_AI_node[i]));
    } else {
      ret += R::dnorm(p_AI_node[i], p_AI_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    if (i == 0) {
      ret += -p_AD_node[i] -2*log(1 + exp(-p_AD_node[i]));
    } else {
      ret += R::dnorm(p_AD_node[i], p_AD_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    if (i == 0) {
      ret += -p_ID_node[i] -2*log(1 + exp(-p_ID_node[i]));
    } else {
      ret += R::dnorm(p_ID_node[i], p_ID_node[i-1], k, true);
    }
  }
  for (int i = 0; i < s_ptr->n_node; ++i) {
    if (i == 0) {
      ret += -p_SD_node[i] -2*log(1 + exp(-p_SD_node[i]));
    } else {
      ret += R::dnorm(p_SD_node[i], p_SD_node[i-1], k, true);
    }
  }
  
  return ret;
}

//------------------------------------------------
// get density of delay distribution on day x
double Particle::get_delay_density(int x, double m, double s) {
#define USE_LOOKUP
#ifdef USE_LOOKUP
  int m_index = floor(m * 100);
  int s_index = floor(s);
  if ((m_index < 0) || (m_index > 2000) || (s_index < 0) || (s_index > 9) || (x < 0)) {
    print("get_delay_density outside lookup range");
    print(x, m, s, m_index, s_index);
    Rcpp::stop("");
  }
  if (x > 100) {
    return 1e-200;
  }
  return s_ptr->gamma_density_lookup[m_index][s_index][x];
#else
  double ret = R::pgamma(x + 1, s, m/s, true, false) - R::pgamma(x, s, m/s, true, false);
  if (ret < 1e-200) {
    ret = 1e-200;
  }
  return ret;
#endif
}

