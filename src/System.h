
#pragma once

#include <Rcpp.h>

#include <vector>

//------------------------------------------------
// class holding all data, parameters and functions
class System {
  
public:
  // PUBLIC OBJECTS
  
  // individual-level data
  std::vector<int> age_group;
  std::vector<int> icu;
  std::vector<int> stepdown;
  std::vector<int> date_admission;
  std::vector<int> date_icu;
  std::vector<int> date_stepdown;
  std::vector<int> date_final_outcome;
  std::vector<int> final_outcome;
  std::vector<int> date_censor;
  int n_ind;
  
  // misc data
  int n_sitrep;
  std::vector<int> node_x;
  int n_node;
  int lookup_max;
  int n_date_sitrep;
  int n_age_sitrep;
  int n_age_indlevel;
  std::vector<double> rel_prop;
  std::vector<int> map_age_indlevel;
  std::vector<std::vector<int>> map_age_sitrep;
  std::vector<int> update_param;
  std::vector<int> update_age;
  
  // lookup tables
  int lookup_n_m;
  int lookup_n_s;
  std::vector<std::vector<std::vector<double>>> density_gamma;
  std::vector<std::vector<std::vector<double>>> tail_gamma;
  
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
  int rungs;
  int chain;
  
  // misc parameters
  bool pb_markdown;
  bool silent;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  System() {};
  
  // public methods
  void load(Rcpp::List args);
  
};
