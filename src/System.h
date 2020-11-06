
#pragma once

#include <Rcpp.h>

#include <vector>

//------------------------------------------------
// class holding all data, parameters and functions
class System {
  
public:
  // PUBLIC OBJECTS
  
  // misc data
  int max_indlevel_age;
  
  // age splines
  std::vector<double> node_x;
  int n_node;
  
  // individual-level data
  std::vector<int> p_AI_numer;
  std::vector<int> p_AI_denom;
  std::vector<int> p_AD_numer;
  std::vector<int> p_AD_denom;
  std::vector<int> p_ID_numer;
  std::vector<int> p_ID_denom;
  std::vector<int> m_AI_count;
  std::vector<int> m_AD_count;
  std::vector<int> m_AC_count;
  std::vector<int> m_ID_count;
  std::vector<int> m_IS_count;
  std::vector<int> m_SC_count;
  
  // model parameters
  std::vector<double> theta_min;
  std::vector<double> theta_max;
  std::vector<double> theta_init;
  std::vector<int> trans_type;
  std::vector<bool> skip_param;
  int d;
  
  // MCMC parameters
  int burnin;
  int samples;
  std::vector<double> beta_vec;
  int rungs;
  int chain;
  
  // misc parameters
  bool pb_markdown;
  bool silent;
  
  // lookup tables
  std::vector<std::vector<std::vector<double>>> gamma_density_lookup;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args);
  
};
