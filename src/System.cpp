
#include "System.h"
#include "misc_v7.h"

using namespace std;

void System::load(Rcpp::List args) {
  
  // split argument lists
  Rcpp::List args_params = args["args_params"];
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::List args_progress = args["args_progress"];
  Rcpp::List args_progress_burnin = args_progress["pb_burnin"];
  
  // data list
  Rcpp::List data_list = args_params["data_list"];
  
  // individual-level data
  Rcpp::List indlevel_list = data_list["indlevel"];
  age_group = rcpp_to_vector_int(indlevel_list["age_group"]);
  icu = rcpp_to_vector_int(indlevel_list["icu"]);
  stepdown = rcpp_to_vector_int(indlevel_list["stepdown"]);
  date_admission = rcpp_to_vector_int(indlevel_list["date_admission"]);
  date_icu = rcpp_to_vector_int(indlevel_list["date_icu"]);
  date_stepdown = rcpp_to_vector_int(indlevel_list["date_stepdown"]);
  date_final_outcome = rcpp_to_vector_int(indlevel_list["date_final_outcome"]);
  final_outcome = rcpp_to_vector_int(indlevel_list["final_outcome_numeric"]);
  date_censor = rcpp_to_vector_int(indlevel_list["date_censor"]);
  n_ind = age_group.size();
  
  // misc data
  n_sitrep = rcpp_to_int(data_list["n_sitrep"]);
  node_x = rcpp_to_vector_int(data_list["node_x"]);
  n_node = node_x.size();
  lookup_max = rcpp_to_int(data_list["lookup_max"]);
  n_date_sitrep = rcpp_to_int(data_list["n_date_sitrep"]);
  n_age_sitrep = rcpp_to_int(data_list["n_age_sitrep"]);
  n_age_indlevel = rcpp_to_int(data_list["n_age_indlevel"]);
  rel_prop = rcpp_to_vector_double(data_list["rel_prop"]);
  map_age_indlevel = rcpp_to_vector_int(data_list["map_age_indlevel"]);
  map_age_sitrep = rcpp_to_matrix_int(data_list["map_age_sitrep"]);
  update_param = rcpp_to_vector_int(data_list["update_param"]);
  update_age = rcpp_to_vector_int(data_list["update_age"]);
  
  // lookup tables
  lookup_n_m = 100;
  lookup_n_s = 100;
  density_gamma = vector<vector<vector<double>>>(lookup_n_m, vector<vector<double>>(lookup_n_s, vector<double>(lookup_max)));
  tail_gamma = vector<vector<vector<double>>>(lookup_n_m, vector<vector<double>>(lookup_n_s, vector<double>(lookup_max)));
  for (int i = 0; i < lookup_n_m; ++i) {
    double m = 30.0 * i / double(lookup_n_m - 1);
    for (int j = 0; j < lookup_n_s; ++j) {
      double s = j / double(lookup_n_s - 1);
      for (int k = 0; k < lookup_max; ++k) {
        density_gamma[i][j][k] = R::pgamma(k + 1, 1.0/(s*s), m*s*s, true, false) -
                                 R::pgamma(k, 1.0/(s*s), m*s*s, true, false);
        tail_gamma[i][j][k] = R::pgamma(k + 1, 1.0/(s*s), m*s*s, false, false);
      }
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
  rungs = rcpp_to_int(args_params["rungs"]);
  chain = rcpp_to_int(args_params["chain"]);
  
  // misc parameters
  pb_markdown = rcpp_to_bool(args_params["pb_markdown"]);
  silent = rcpp_to_bool(args_params["silent"]);
  
}
