
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
  
  // mean durations
  std::vector<double> m_AI_node;
  std::vector<double> m_AI;
  std::vector<double> m_AD_node;
  std::vector<double> m_AD;
  std::vector<double> m_AC_node;
  std::vector<double> m_AC;
  std::vector<double> m_ID_node;
  std::vector<double> m_ID;
  std::vector<double> m_IS_node;
  std::vector<double> m_IS;
  std::vector<double> m_SC_node;
  std::vector<double> m_SC;
  
  // coefficients of variation of durations
  double s_AI;
  double s_AD;
  double s_AC;
  double s_ID;
  double s_IS;
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
  double get_delay_tail(int x, double m, double s);
  void phi_prop_to_theta_prop(int i);
  void theta_to_phi();
  double get_adjustment(int i);
  
};
