
#include "main.h"
#include "misc_v7.h"
#include "probability_v10.h"
#include "System.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// run MCMC
Rcpp::List run_mcmc_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // create sytem object and load args
  System s;
  s.load(args);
  
  // extract R utility functions that will be called from within MCMC
  Rcpp::List args_functions = args["args_functions"];
  Rcpp::Function test_convergence = args_functions["test_convergence"];
  Rcpp::Function update_progress = args_functions["update_progress"];
  
  // extract progress bar objects
  Rcpp::List args_progress = args["args_progress"];
  
  // local copies of some parameters for convenience
  int d = s.d;
  vector<double> beta_vec = s.beta_vec;
  int rungs = s.rungs;
  
  // initialise vector of particles
  vector<Particle> particle_vec(rungs);
  for (int r = 0; r < rungs; ++r) {
    particle_vec[r].init(s);
  }
  
  // specify rung order
  vector<int> rung_order = seq_int(0, rungs-1);
  
  // option to return model fit
  if (s.return_fit) {
    
    // return model fit
    return Rcpp::List::create(Rcpp::Named("admission_incidence") = particle_vec[0].admission_incidence,
                              Rcpp::Named("deaths_incidence") = particle_vec[0].deaths_incidence,
                              Rcpp::Named("discharges_incidence") = particle_vec[0].discharges_incidence,
                              Rcpp::Named("general_prevalence") = particle_vec[0].general_prevalence,
                              Rcpp::Named("critical_prevalence") = particle_vec[0].critical_prevalence);
  }
  
  // objects for storing loglikelihood and theta values over iterations
  vector<vector<double>> loglike_burnin(rungs, vector<double>(s.burnin));
  vector<vector<double>> logprior_burnin(rungs, vector<double>(s.burnin));
  vector<vector<vector<double>>> theta_burnin(rungs, vector<vector<double>>(s.burnin, vector<double>(d)));
  vector<vector<double>> loglike_sampling(rungs, vector<double>(s.samples));
  vector<vector<double>> logprior_sampling(rungs, vector<double>(s.samples));
  vector<vector<vector<double>>> theta_sampling(rungs, vector<vector<double>>(s.samples, vector<double>(d)));
  
  // specify stored values at first iteration. Ensures that user-defined initial
  // values are the first stored values
  for (int r = 0; r < rungs; ++r) {
    loglike_burnin[r][0] = particle_vec[r].loglike;
    logprior_burnin[r][0] = particle_vec[r].logprior;
    theta_burnin[r][0] = particle_vec[r].theta;
  }
  
  // store Metropolis coupling acceptance rates
  vector<int> mc_accept_burnin(rungs - 1);
  vector<int> mc_accept_sampling(rungs - 1);
  
  
  // ---------- burn-in MCMC ----------
  
  // print message to console
  if (!s.silent) {
    print("MCMC chain", s.chain);
    print("burn-in");
  }
  
  // loop through burn-in iterations
  for (int rep = 1; rep < s.burnin; ++rep) {
    
    // loop through rungs
    for (int r = 0; r < rungs; ++r) {
      
      // update particles
      particle_vec[rung_order[r]].update(beta_vec[rung_order[r]]);
      
      // store results
      loglike_burnin[r][rep] = particle_vec[rung_order[r]].loglike;
      logprior_burnin[r][rep] = particle_vec[rung_order[r]].logprior;
      theta_burnin[r][rep] = particle_vec[rung_order[r]].theta;
    }
    
    // perform Metropolis coupling
    coupling(particle_vec, mc_accept_burnin, true, beta_vec, rung_order);
    
    // update progress bars
    if (!s.silent) {
      int remainder = rep % int(ceil(double(s.burnin)/100));
      if ((remainder == 0 && !s.pb_markdown) || ((rep+1) == s.burnin)) {
        update_progress(args_progress, "pb_burnin", rep+1, s.burnin, false);
        if ((rep+1) == s.burnin) {
          print("");
        }
      }
    }
    
  }  // end burn-in MCMC loop
  
  // print phase diagnostics
  if (!s.silent) {
    double accept_rate = particle_vec[rungs-1].accept_count/double(s.burnin*d);
    Rcpp::Rcout << "acceptance rate: " << round(accept_rate*1000)/10.0 << "%\n";
  }
  
  
  // ---------- sampling MCMC ----------
  
  // print message to console
  if (!s.silent) {
    print("sampling phase");
  }
  
  // reset acceptance count of all rungs
  for (int r = 0; r < rungs; ++r) {
    particle_vec[r].accept_count = 0;
  }
  
  // loop through sampling iterations
  for (int rep = 0; rep < s.samples; ++rep) {
    
    // loop through rungs
    for (int r = 0; r < rungs; ++r) {
      
      // update particles
      particle_vec[rung_order[r]].update(beta_vec[rung_order[r]]);
      
      // store results
      loglike_sampling[r][rep] = particle_vec[rung_order[r]].loglike;
      logprior_sampling[r][rep] = particle_vec[rung_order[r]].logprior;
      theta_sampling[r][rep] = particle_vec[rung_order[r]].theta;
    }
    
    // perform Metropolis coupling
    coupling(particle_vec, mc_accept_sampling, false, beta_vec, rung_order);
    
    // update progress bars
    if (!s.silent) {
      int remainder = rep % int(ceil(double(s.samples)/100));
      if ((remainder == 0 && !s.pb_markdown) || ((rep+1) == s.samples)) {
        update_progress(args_progress, "pb_samples", rep+1, s.samples, false);
        if ((rep+1) == s.samples) {
          print("");
        }
      }
    }
    
  }  // end sampling MCMC loop
  
  // print final diagnostics
  if (!s.silent) {
    double accept_rate = particle_vec[rung_order[rungs-1]].accept_count/double(s.samples*d);
    Rcpp::Rcout << "acceptance rate: " << round(accept_rate*1000)/10.0 << "%\n";
  }
  
  
  // ---------- return ----------
  
  // end timer
  if (!s.silent) {
    print("");
    chrono_timer(t1);
  }
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("loglike_burnin") = loglike_burnin,
                                      Rcpp::Named("logprior_burnin") = logprior_burnin,
                                      Rcpp::Named("theta_burnin") = theta_burnin,
                                      Rcpp::Named("loglike_sampling") = loglike_sampling,
                                      Rcpp::Named("logprior_sampling") = logprior_sampling,
                                      Rcpp::Named("theta_sampling") = theta_sampling,
                                      Rcpp::Named("beta_vec") = beta_vec,
                                      Rcpp::Named("mc_accept_burnin") = mc_accept_burnin,
                                      Rcpp::Named("mc_accept_sampling") = mc_accept_sampling);
  return ret;
}

//------------------------------------------------
// Metropolis-coupling over temperature rungs
void coupling(vector<Particle> &particle_vec, vector<int> &mc_accept, bool adaptive,
              vector<double> &beta_vec, vector<int> &rung_order) {
  
  // get number of rungs
  int rungs = int(beta_vec.size());
  
  // return if single rung
  if (rungs == 1) {
    return;
  }
  
  // loop over rungs, starting with the hottest chain and moving to the cold
  // chain. Each time propose a swap with the next rung up
  for (int i = 1; i < rungs; ++i) {
    
    // define rungs of interest
    int rung1 = rung_order[i-1];
    int rung2 = rung_order[i];
    
    // get log-likelihoods and beta values of two chains in the comparison
    double loglike1 = particle_vec[rung1].loglike;
    double loglike2 = particle_vec[rung2].loglike;
    
    double beta1 = beta_vec[rung1];
    double beta2 = beta_vec[rung2];
    
    // calculate acceptance ratio (still in log space)
    double acceptance = (loglike2*beta1 + loglike1*beta2) - (loglike1*beta1 + loglike2*beta2);
    
    // accept or reject move
    bool accept_move = (log(runif_0_1()) < acceptance);
    
    // implement swap
    if (accept_move) {
      
      // swap beta values
      beta_vec[rung1] = beta2;
      beta_vec[rung2] = beta1;
      
      // swap rung order
      int tmp = rung_order[i-1];
      rung_order[i-1] = rung_order[i];
      rung_order[i] = tmp;
      
      // update acceptance rates
      mc_accept[i-1]++;
      
      // adaptive update
      if (adaptive) {
        
        
        
      }  // end adaptive update
      
    } else {  // if reject move
      
      // adaptive update
      if (adaptive) {
        
      }  // end adaptive update
      
    }
    
    //print_vector(rung_order);
    
  }  // end loop over rungs
  
}


