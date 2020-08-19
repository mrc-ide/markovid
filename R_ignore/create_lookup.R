
# create_lookup.R
#
# Author: Bob Verity
# Date: 2020-08-18
#
# Purpose:
# Create gamma daily density and tail lookup tables to be read in by markovid.

# ------------------------------------------------------------------

# define vectors over which to create lookup
mvec <- seq(0, 20, 0.01)
svec <- seq(0, 1, 0.01)
kvec <- 0:100

# create matrices from marginal vectors
kmat <- matrix(rep(kvec, each = length(svec)), length(svec))
smat <- output <- matrix(rep(svec, length(kvec)), length(svec))

# save values into list
gamma_density <- gamma_tail <- list()
for (i in seq_along(mvec)) {
  m <- mvec[i]
  
  # calculate gamma density and gamma tail
  gamma_density[[i]] <- pgamma(kmat + 1, shape = 1/smat^2, scale = m*smat^2, lower.tail = TRUE) -
                        pgamma(kmat, shape = 1/smat^2, scale = m*smat^2, lower.tail = TRUE)
  gamma_tail[[i]] <- pgamma(kmat + 1, shape = 1/smat^2, scale = m*smat^2, lower.tail = FALSE)
  
  # replace NaN with 0
  gamma_density[[i]][is.na(gamma_density[[i]])] <- 0
  gamma_tail[[i]][is.na(gamma_tail[[i]])] <- 0
  
  # add tiny value to buffer against underflow
  gamma_density[[i]] <- gamma_density[[i]] + 1e-200
  gamma_tail[[i]] <- gamma_tail[[i]] + 1e-200
  
  # split into nested lists
  gamma_density[[i]] <- split(gamma_density[[i]], f = row(gamma_density[[i]]))
  gamma_tail[[i]] <- split(gamma_tail[[i]], f = row(gamma_tail[[i]]))
  
}

# save list objects to file
saveRDS(gamma_density, "inst/extdata/gamma_density.rds")
saveRDS(gamma_tail, "inst/extdata/gamma_tail.rds")
