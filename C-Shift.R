if (!requireNamespace("trustOptim", quietly = TRUE))  install.packages("trustOptim")
library(trustOptim)

################### functions ###################
# Trace objective function
objFunct_T <- function(alpha, Cov_pert, geneNo){
  D <- diag(alpha) %*% matrix(1, ncol = geneNo, nrow = geneNo)
  M <- Cov_pert + D + t(D)
  M_inv_hat <- solve(M, rep(1, geneNo), tol = 10 ^-25)
  V <- 1 / sum(M_inv_hat)
  C_alpha <- M - matrix(V, ncol = geneNo, nrow = geneNo)
  sum(diag(C_alpha))
}

# Trace gradient function
grad.objFunct_T <- function(alpha, Cov_pert, geneNo){
  D <- diag(alpha) %*% matrix(1, ncol = geneNo, nrow = geneNo)
  M <- Cov_pert + D + t(D)
  a <- sum(alpha)
  odin <- rep(1, geneNo)
  M_inv_hat <- solve(M, rep(1, geneNo), tol = 10 ^-25)
  V <- 1 / sum(M_inv_hat)
  grad <- 2 * (odin - geneNo * V * M_inv_hat)
  grad
}

# Frob Norm objective function
objFunct_F <- function(alpha, Cov_pert, geneNo){
  D <- diag(alpha) %*% matrix(1, ncol = geneNo, nrow = geneNo)
  M <- Cov_pert + D + t(D)
  M_inv_hat <- solve(M, rep(1, geneNo), tol = 10 ^-25)
  V <- 1 / sum(M_inv_hat)
  C_alpha <- M - matrix(V, ncol = geneNo, nrow = geneNo)
  (norm(C_alpha, "F")) ^ 2
}

# Frob Norm gradient function
grad.objFunct_F <- function(alpha, Cov_pert, geneNo){
  D <- diag(alpha) %*% matrix(1, ncol = geneNo, nrow = geneNo)
  M <- Cov_pert + D + t(D)
  a <- sum(alpha)
  odin <- rep(1, geneNo)
  M_inv_hat <- solve(M, rep(1, geneNo), tol = 10 ^-25)
  V <- 1 / sum(M_inv_hat)
  grad <- geneNo * alpha + ((geneNo * V) ^ 2 - sum(Cov_pert) * V - 2 * geneNo * V * a) * M_inv_hat + 
    Cov_pert %*% odin + (a - geneNo * V) * odin
  grad <- 4*grad
  grad
}

#################### Calculate OUTPUT matrices: Cov_new #################
CShift_fn <- function(Cov_obs, use_trace = T){
  geneNo <- dim(Cov_obs)[1] #extract number of genes
  
  # Make matrix pos def
  if (qr(Cov_obs)$rank == geneNo) {
    Fm <- matrix(0, ncol = geneNo, nrow = geneNo)
  }else{
    Fm <- diag(runif(geneNo))
  }
  Cov_pert <- Cov_obs + Fm
  
  # Run optimization
  alpha <- rep(0, geneNo)
  if(use_trace == T) ooo <- trust.optim(alpha, objFunct_T, grad.objFunct_T, method = "SR1", 
                                       control = list(maxit = 2000, 
                                                      stop.trust.radius = 0.000000000001, 
                                                      prec = 0.0001),
                                       Cov_pert = Cov_pert, geneNo = geneNo)
  if(use_trace == F) ooo <- trust.optim(alpha, objFunct_F, grad.objFunct_F, method = "SR1", 
                                       control = list(maxit = 2000, 
                                                      stop.trust.radius = 0.000000000001, 
                                                      prec = 0.0001),
                                       Cov_pert = Cov_pert, geneNo = geneNo)
  alpha <- ooo$solution
  
  # Use final alpha to calculate new covariance matrix
  D <- diag(alpha) %*% matrix(1, ncol = geneNo, nrow = geneNo)
  M <- Cov_pert + D + t(D)
  M_inv_hat <- solve(M, tol = 10 ^-25)
  V <- 1 / sum(M_inv_hat)
  C_alpha <- M - V * matrix(1, ncol = geneNo, nrow = geneNo)
  Cov_new <- C_alpha - Fm # Remove previously added variances
  
  # because of approximation errors the resulting matrix can be non-symmetric, thus we 
  # make these corrections
  C <- Cov_new * upper.tri(Cov_new)
  Cov_new <- C + t(C) + diag(diag(Cov_new))
  Cov_new
}
