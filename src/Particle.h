
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
  
  // spline parameters
  std::vector<std::vector<double>> node_y;
  
  // population proportions
  std::vector<double> scale_rel_prop;
  
  // rescaling parameters
  std::vector<double> scale_p_AI;
  std::vector<double> scale_p_AD;
  std::vector<double> scale_p_ID;
  
  // lab test weights
  std::vector<double> pos_on_day;
  std::vector<double> neg_by_day;
  std::vector<double> pos_by_day;
  
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
  std::vector<double> m_AC_node;
  std::vector<double> m_AC;
  
  // dynamic lookup tables for interval distributions
  std::vector<double> density_AL;
  std::vector<double> density_AI;
  std::vector<double> density_AD;
  std::vector<std::vector<double>> density_AC;
  //std::vector<double> density_AC;
  std::vector<double> density_ID;
  std::vector<double> density_IS;
  std::vector<double> density_SC;
  
  // dynamic lookup tables for complementary cumulative density (ccdf) distributions
  std::vector<double> tail_AI;
  std::vector<double> tail_AD;
  std::vector<std::vector<double>> tail_AC;
  //std::vector<double> tail_AC;
  std::vector<double> tail_ID;
  std::vector<double> tail_IS;
  std::vector<double> tail_SC;
  
  // objects for storing progression over all stratification
  std::vector<std::vector<std::vector<double>>> admission_incidence;
  std::vector<std::vector<std::vector<double>>> deaths_incidence;
  std::vector<std::vector<std::vector<double>>> discharges_incidence;
  std::vector<std::vector<std::vector<double>>> general_prevalence;
  std::vector<std::vector<std::vector<double>>> critical_prevalence;
  
  std::vector<double> delta_stepup;
  std::vector<double> delta_stepdown;
  std::vector<double> delta_deaths_general;
  std::vector<double> delta_discharges_general;
  std::vector<double> delta_open_general;
  std::vector<double> delta_deaths_critical;
  std::vector<double> delta_open_critical;
  
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
  double get_loglike(std::vector<double> &theta, int theta_i, bool quick_exit = false);
  double get_logprior(std::vector<double> &theta, int theta_i);
  
  // other public methods
  double get_delay_density(int x, double m, double s);
  double get_delay_tail(int x, double m, double s);
  void update_gamma_density(std::vector<double> &density_vec, double m, double s);
  void update_gamma_tail(std::vector<double> &tail_vec, double m, double s);
  void phi_prop_to_theta_prop(int i);
  void theta_to_phi();
  double get_adjustment(int i);
  
};
