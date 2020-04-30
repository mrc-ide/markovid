
#include "System.h"
#include "Particle.h"

#include <Rcpp.h>

//------------------------------------------------
// run MCMC
// [[Rcpp::export]]
Rcpp::List run_mcmc_cpp(Rcpp::List args);

//------------------------------------------------
// Metropolis-coupling over temperature rungs
void coupling(std::vector<Particle> &particle_vec, std::vector<int> &mc_accept);
