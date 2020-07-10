
#include "System.h"
#include "misc_v10.h"

using namespace std;

void System::load(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_params = args["args_params"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  Rcpp::List args_progress_burnin = args_progress["pb_burnin"];
  
  // option to return model fit
  return_fit = rcpp_to_bool(args_params["return_fit"]);
  
  // data list
  Rcpp::List data_list = args_params["data_list"];
  
  // misc data
  lookup_max = rcpp_to_int(data_list["lookup_max"]);
  n_region = rcpp_to_int(data_list["n_region"]);
  n_age_sitrep = rcpp_to_int(data_list["n_age_sitrep"]);
  n_date_sitrep = rcpp_to_int(data_list["n_date_sitrep"]);
  max_indlevel_age = rcpp_to_int(data_list["max_indlevel_age"]);
  
  // age weights
  age_weights = rcpp_to_matrix_double(data_list["age_weights"]);
  age_values = rcpp_to_matrix_int(data_list["age_values"]);
  
  // age splines
  p_AI_nodex = rcpp_to_vector_double(data_list["p_AI_nodex"]);
  p_AI_noden = p_AI_nodex.size();
  p_AD_nodex = rcpp_to_vector_double(data_list["p_AD_nodex"]);
  p_AD_noden = p_AI_nodex.size();
  p_ID_nodex = rcpp_to_vector_double(data_list["p_ID_nodex"]);
  p_ID_noden = p_AI_nodex.size();
  
  m_AC_nodex = rcpp_to_vector_double(data_list["m_AC_nodex"]);
  m_AC_noden = m_AC_nodex.size();
  
  // individual-level data
  Rcpp::List indlevel_list = data_list["indlevel"];
  age = rcpp_to_vector_int(indlevel_list["age"]);
  icu = rcpp_to_vector_int(indlevel_list["icu"]);
  stepdown = rcpp_to_vector_int(indlevel_list["stepdown"]);
  date_admission = rcpp_to_vector_int(indlevel_list["date_admission"]);
  date_icu = rcpp_to_vector_int(indlevel_list["date_icu"]);
  date_stepdown = rcpp_to_vector_int(indlevel_list["date_stepdown"]);
  date_final_outcome = rcpp_to_vector_int(indlevel_list["date_final_outcome"]);
  final_outcome = rcpp_to_vector_int(indlevel_list["final_outcome_numeric"]);
  date_censor = rcpp_to_vector_int(indlevel_list["date_censor"]);
  n_ind = age.size();
  
  // sitrep data
  Rcpp::List sitrep_list = data_list["sitrep"];
  daily_influx = vector<vector<vector<int>>>(n_region);
  daily_influx = vector<vector<vector<int>>>(n_region, vector<vector<int>>(n_age_sitrep));
  new_deaths = vector<vector<vector<int>>>(n_region, vector<vector<int>>(n_age_sitrep, vector<int>(n_date_sitrep)));
  new_discharges = vector<vector<vector<int>>>(n_region, vector<vector<int>>(n_age_sitrep, vector<int>(n_date_sitrep)));
  total_general = vector<vector<vector<int>>>(n_region, vector<vector<int>>(n_age_sitrep, vector<int>(n_date_sitrep)));
  total_critical = vector<vector<vector<int>>>(n_region, vector<vector<int>>(n_age_sitrep, vector<int>(n_date_sitrep)));
  
  for (int i = 0; i < n_region; ++i) {
    Rcpp::List sitrep_i = sitrep_list[i];
    for (int j = 0; j < n_age_sitrep; ++j) {
      Rcpp::List sitrep_j = sitrep_i[j];
      
      daily_influx[i][j] = rcpp_to_vector_int(sitrep_j["daily_influx"]);
      new_deaths[i][j] = rcpp_to_vector_int(sitrep_j["deaths"]);
      new_discharges[i][j] = rcpp_to_vector_int(sitrep_j["new_discharges"]);
      total_general[i][j] = rcpp_to_vector_int(sitrep_j["total_general"]);
      total_critical[i][j] = rcpp_to_vector_int(sitrep_j["total_hdu_icu"]);
    }
  }
  
  // model parameters
  theta_min = rcpp_to_vector_double(args_params["theta_min"]);
  theta_max = rcpp_to_vector_double(args_params["theta_max"]);
  theta_init = rcpp_to_vector_double(args_params["theta_init"]);
  trans_type = rcpp_to_vector_int(args_params["trans_type"]);
  skip_param = rcpp_to_vector_bool(args_params["skip_param"]);
  d = int(theta_min.size());
  
  // MCMC parameters
  burnin = rcpp_to_int(args_params["burnin"]);
  samples = rcpp_to_int(args_params["samples"]);
  beta_vec = rcpp_to_vector_double(args_params["beta_vec"]);
  rungs = beta_vec.size();
  chain = rcpp_to_int(args_params["chain"]);
  
  // misc parameters
  pb_markdown = rcpp_to_bool(args_params["pb_markdown"]);
  silent = rcpp_to_bool(args_params["silent"]);
  
  // populate lookup tables
  int n_m = 2001;
  int n_s = 101;
  gamma_density_lookup = std::vector<std::vector<std::vector<double>>>(n_m, std::vector<std::vector<double>>(n_s, std::vector<double>(lookup_max)));
  gamma_tail_lookup = std::vector<std::vector<std::vector<double>>>(n_m, std::vector<std::vector<double>>(n_s, std::vector<double>(lookup_max)));
  print("initializing lookup tables");
  for (int i = 0; i < n_m; ++i) {
    double m = double(i) / 100.0;
    for (int j = 0; j < n_s; ++j) {
      double s = double(j) / 100.0;
      for (int k = 0; k < lookup_max; ++k) {
        gamma_density_lookup[i][j][k] = R::pgamma(k + 1, 1.0/(s*s), m*s*s, true, false) -
                                        R::pgamma(k, 1.0/(s*s), m*s*s, true, false);
        gamma_tail_lookup[i][j][k] = R::pgamma(k + 1, 1.0/(s*s), m*s*s, false, false);
      }
    }
  }
  
}
