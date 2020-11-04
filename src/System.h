
#pragma once

#include <Rcpp.h>

#include <vector>

//------------------------------------------------
// class holding all data, parameters and functions
class System {
  
public:
  // PUBLIC OBJECTS
  
  // option to return model fit
  bool return_fit;
  
  // misc data
  int lookup_max;
  int n_region;
  int n_age_sitrep;
  int n_date_sitrep;
  int max_indlevel_age;
  
  // age weights
  std::vector< std::vector<double> > age_weights;
  std::vector< std::vector<int> > age_values;
  
  // age splines
  std::vector<double> p_AI_nodex;
  int p_AI_noden;
  std::vector<double> p_AD_nodex;
  int p_AD_noden;
  std::vector<double> p_ID_nodex;
  int p_ID_noden;
  
  std::vector<double> m_AI_nodex;
  int m_AI_noden;
  std::vector<double> m_AD_nodex;
  int m_AD_noden;
  std::vector<double> m_AC_nodex;
  int m_AC_noden;
  std::vector<double> m_ID_nodex;
  int m_ID_noden;
  std::vector<double> m_IS_nodex;
  int m_IS_noden;
  std::vector<double> m_SC_nodex;
  int m_SC_noden;
  
  // individual-level data
  std::vector<int> p_AI_numer;
  std::vector<int> p_AI_denom;
  std::vector<int> p_AD_numer;
  std::vector<int> p_AD_denom;
  std::vector<int> p_ID_numer;
  std::vector<int> p_ID_denom;
  std::vector<std::vector<int>> m_AI_count;
  std::vector<std::vector<int>> m_AD_count;
  std::vector<std::vector<int>> m_AC_count;
  std::vector<std::vector<int>> m_ID_count;
  std::vector<std::vector<int>> m_IS_count;
  std::vector<std::vector<int>> m_SC_count;
  
  std::vector<int> age;
  std::vector<int> icu;
  std::vector<int> stepdown;
  std::vector<int> date_admission;
  std::vector<int> date_icu;
  std::vector<int> date_stepdown;
  std::vector<int> date_final_outcome;
  std::vector<int> final_outcome;
  std::vector<int> date_censor;
  int n_ind;
  
  // sitrep data
  std::vector<std::vector<std::vector<int>>> daily_influx;
  std::vector<std::vector<std::vector<int>>> new_deaths;
  std::vector<std::vector<std::vector<int>>> new_discharges;
  std::vector<std::vector<std::vector<int>>> total_general;
  std::vector<std::vector<std::vector<int>>> total_critical;
  
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
  bool sitrep_loglike;
  
  // lookup tables
  std::vector<std::vector<std::vector<double>>> gamma_density_lookup;
  std::vector<std::vector<std::vector<double>>> gamma_tail_lookup;
  
  std::vector<std::vector<double>> pgamma_lookup;
  
  size_t n_threads;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args);
  
};
