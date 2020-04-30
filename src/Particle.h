
#pragma once

#include "System.h"
#include "misc_v7.h"
#include "probability_v3.h"

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
  
  // beta_raised stores values of beta (the thermodynamic power), raised to the
  // power GTI_pow
  double beta_raised;
  
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
  
  // initialise everything EXCEPT FOR likelihood and prior values
  void init(System &s, double beta_raised);
  
  // initialise likelihood and prior values
  void init_like();
  
  // update theta[i] via univariate Metropolis-Hastings
  void update();
  
  double get_loglike();
  double get_logprior();
  
  // other public methods
  void propose_phi(int i);
  void phi_prop_to_theta_prop(int i);
  void theta_to_phi();
  double get_adjustment(int i);
  
};
