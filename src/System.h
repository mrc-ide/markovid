
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
  int n_region;
  std::vector<int> node_x;
  int n_node;
  int n_spline;
  int lookup_max;
  int n_date_sitrep;
  int n_age_sitrep;
  int n_age_indlevel;
  std::vector<double> rel_prop;
  std::vector<int> map_age_indlevel;
  std::vector<std::vector<int>> map_age_sitrep;
  
  // update rules
  std::vector<int> update_density;
  std::vector<int> update_region;
  std::vector<int> update_indlevel_age;
  std::vector<int> update_sitrep_age;
  
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
