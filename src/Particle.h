
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
  
  // thermodynamic power
  double beta;
  
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
  void init(System &s, double beta);
  
  // update theta[i] via univariate Metropolis-Hastings
  void update();
  
  // loglikelihood and logprior
  double get_loglike(std::vector<double> &theta);
  double get_logprior(std::vector<double> &theta);
  
  // other public methods
  void phi_prop_to_theta_prop(int i);
  void theta_to_phi();
  double get_adjustment(int i);
  
};
