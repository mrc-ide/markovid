
#pragma once

#include "System.h"
#include "misc_v10.h"
#include "probability_v10.h"

#include <Rcpp.h>

//------------------------------------------------
// class defining MCMC particle
class Particle {
  
public:
  // PUBLIC OBJECTS
  
  // pointer to system object
  System * s_ptr;
  
  // local copies of some parameters for convenience
  int d;
  
  // vector over ages for cubic splines
  std::vector<double> age_seq;
  
  // transition probabilities
  std::vector<double> p_AI_node;
  std::vector<double> p_AI;
  std::vector<double> p_AD_node;
  std::vector<double> p_AD;
  std::vector<double> p_ID_node;
  std::vector<double> p_ID;
  std::vector<double> p_SD_node;
  std::vector<double> p_SD;
  
  // mean durations
  double m_AI;
  double m_AD;
  double m_AC;
  double m_ID;
  double m_I1S;
  double m_I2S;
  double m_SD;
  double m_SC;
  
  // coefficients of variation of durations
  double s_AI;
  double s_AD;
  double s_AC;
  double s_ID;
  double s_I1S;
  double s_I2S;
  double s_SD;
  double s_SC;
  
  // theta is the parameter vector in natural space
  std::vector<double> theta;
  std::vector<double> theta_prop;
  
  // phi is a vector of transformed parameters
  std::vector<double> phi;
  std::vector<double> phi_prop;
  
  // proposal parameters
  std::vector<double> bw;
  std::vector<int> bw_index;
  double bw_stepsize;
  
  // likelihoods and priors
  double loglike;
  double loglike_prop;
  double logprior;
  double logprior_prop;
  
  // store acceptance rates
  int accept_count;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Particle() {};
  
  // initialise
  void init(System &s);
  
  // update theta[i] via univariate Metropolis-Hastings
  void update(double beta);
  
  // loglikelihood and logprior
  double get_loglike(std::vector<double> &theta, int theta_i);
  double get_logprior(std::vector<double> &theta, int theta_i);
  
  // other public methods
  double get_delay_density(int x, double m, double s);
  void phi_prop_to_theta_prop(int i);
  void theta_to_phi();
  double get_adjustment(int i);
  
};
