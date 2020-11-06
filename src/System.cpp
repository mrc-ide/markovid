
#include "System.h"
#include "misc_v10.h"

using namespace std;

void System::load(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_params = args["args_params"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  Rcpp::List args_progress_burnin = args_progress["pb_burnin"];
  
#define USE_LOOKUP
#ifdef USE_LOOKUP
  Rcpp::List args_lookup_density = args["args_lookup_density"];
#endif
  
  // data list
  Rcpp::List data_list = args_params["data_list"];
  
  // misc data
  max_indlevel_age = rcpp_to_int(data_list["max_indlevel_age"]);
  
  // age splines
  node_x = rcpp_to_vector_double(data_list["node_x"]);
  n_node = node_x.size();
  
  // individual-level data
  Rcpp::List indlevel_list = data_list["indlevel"];
  
  p_AI_numer = rcpp_to_vector_int(indlevel_list["p_AI_numer"]);
  p_AI_denom = rcpp_to_vector_int(indlevel_list["p_AI_denom"]);
  p_AD_numer = rcpp_to_vector_int(indlevel_list["p_AD_numer"]);
  p_AD_denom = rcpp_to_vector_int(indlevel_list["p_AD_denom"]);
  p_ID_numer = rcpp_to_vector_int(indlevel_list["p_ID_numer"]);
  p_ID_denom = rcpp_to_vector_int(indlevel_list["p_ID_denom"]);
  m_AI_count = rcpp_to_matrix_int(indlevel_list["m_AI_count"]);
  m_AD_count = rcpp_to_matrix_int(indlevel_list["m_AD_count"]);
  m_AC_count = rcpp_to_matrix_int(indlevel_list["m_AC_count"]);
  m_ID_count = rcpp_to_matrix_int(indlevel_list["m_ID_count"]);
  m_IS_count = rcpp_to_matrix_int(indlevel_list["m_IS_count"]);
  m_SC_count = rcpp_to_matrix_int(indlevel_list["m_SC_count"]);
  
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
#ifdef USE_LOOKUP
  int n_m = args_lookup_density.size();
  gamma_density_lookup = std::vector<std::vector<std::vector<double>>>(n_m);
  for (int i = 0; i < n_m; ++i) {
    gamma_density_lookup[i] = rcpp_to_matrix_double(args_lookup_density[i]);
  }
#endif
  
}
