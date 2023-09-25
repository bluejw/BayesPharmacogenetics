##########################################################################################
#                                   MCMC Main                                            #
##########################################################################################

library(Rcpp)
library(RcppArmadillo)
library(data.tree)
library(Matrix)
library(MCMCpack)
library(GIGrvg)

sourceCpp("MCMC_Rcpp_Functions.cpp") 
source("MCMC_R_Functions.R")
source("MCMC_Update_Tree.R")
load("Treatment_History_Data.Rdata")
source("Data_Generate.R") 


# MCMC Setup
Nit <- 5000 # number of MCMC iterations
burn.in <- 2500 # burn-in period
thin.fac <- 5 # thinning factor 
post_index <- seq(burn.in+1, Nit, by=thin.fac) # index of posterior samples
post_num <- (Nit-burn.in)/thin.fac # number of posterior samples

mcmc <- NULL # list of MCMC samples
mcmc$theta <- array(NA, dim=c(Nit, L_max, S_tilde))
mcmc$etheta <- matrix(0, nrow=Nit, ncol=S_tilde)
mcmc$Btheta <- array(0, dim=c(Nit, S_tilde, S_tilde))
mcmc$lambda <- matrix(NA, nrow=Nit, ncol=S_minus)
mcmc$tau <- rep(NA, Nit)
mcmc$nu <- matrix(NA, nrow=Nit, ncol=S_minus)
mcmc$psi <- rep(NA, Nit)
mcmc$omega <- array(NA, dim=c(Nit, n, max(J)))
mcmc$sigma2_omega <- matrix(NA, nrow=Nit, ncol=L_max)
mcmc$rho <- matrix(NA, nrow=Nit, ncol=L_max)
mcmc$sigma2 <- rep(NA, Nit)
mcmc$tree <- c()
mcmc$L <- rep(NA, Nit)
mcmc$R_tilde <- matrix(NA, nrow=Nit, ncol=n*max(J))
mcmc$logpost <- rep(NA, Nit) 


# Initialize
set.seed(124)
initial <- init()
mcmc$tree <- c(initial$tree)
mcmc$L[1] <- initial$L
mcmc$R_tilde[1,] <- initial$R_tilde
mcmc$theta[1,1:mcmc$L[1],] <- initial$theta
mcmc$etheta[1,] <- initial$etheta
mcmc$Btheta[1,,] <- initial$Btheta
mcmc$lambda[1,] <- initial$lambda
mcmc$tau[1] <- initial$tau
mcmc$nu[1,] <- initial$nu
mcmc$psi[1] <- initial$psi
mcmc$omega[1,,] <- initial$omega
omega_tilde <- as.vector(mcmc$omega[1,,])
mcmc$sigma2_omega[1,1:mcmc$L[1]] <- initial$sigma2_omega
mcmc$rho[1,1:mcmc$L[1]] <- initial$rho
mcmc$sigma2[1] <- initial$sigma2


# Start of the Chain
start.time = proc.time()

for (nit in 2:Nit){
  
  print(nit)
  
  # Update tree
  tree_update <- update_tree(Clone(mcmc$tree[[nit-1]]), mcmc$sigma2_omega[nit-1,1:mcmc$L[nit-1]], 
                             mcmc$rho[nit-1,1:mcmc$L[nit-1]],mcmc$sigma2[nit-1], 
                             mcmc$etheta[nit-1,], mcmc$Btheta[nit-1,,])
  print(tree_update, "split", "rule", "prior", "node_prior", "set")
  mcmc$tree <- c(mcmc$tree, tree_update)
  mcmc$R_tilde[nit,] <- tree_update$index
  mcmc$L[nit] <- length(unique(mcmc$R_tilde[nit,]))-1
  mcmc$logpost[nit] <- tree_update$logpost 
  mcmc$sigma2_omega[nit,1:mcmc$L[nit]] <- tree_update$sigma2_omega
  mcmc$rho[nit,1:mcmc$L[nit]] <- tree_update$rho
  
  # Update theta
  R_tilde <- mcmc$R_tilde[nit,]
  R_tilde_index <- sort(unique(R_tilde[R_tilde!=0])) # unique index of tree terminal nodes
  R_tilde_num <- as.numeric(factor(R_tilde, levels=R_tilde_index)) # numerical index of the tree terminal node
  mcmc$theta[nit,1:mcmc$L[nit],] <- update_theta(mcmc$L[nit], mcmc$R_tilde[nit,], 
                                                 mcmc$sigma2_omega[nit,1:mcmc$L[nit]], mcmc$rho[nit,1:mcmc$L[nit]],
                                                 mcmc$sigma2[nit-1], mcmc$etheta[nit-1,], mcmc$Btheta[nit-1,,])
  
  # Update lambda, tau, nu, psi
  mcmc$lambda[nit,] <- update_lambda(mcmc$L[nit], mcmc$theta[nit,1:mcmc$L[nit],2:S_tilde], mcmc$tau[nit-1], mcmc$nu[nit-1,])
  mcmc$tau[nit] <- update_tau(mcmc$L[nit], mcmc$theta[nit,1:mcmc$L[nit],2:S_tilde], mcmc$lambda[nit,], mcmc$psi[nit-1])
  mcmc$nu[nit,] <- update_nu(mcmc$lambda[nit,])
  mcmc$psi[nit] <- update_psi(mcmc$tau[nit])
  mcmc$Btheta[nit,1,1] <- 10; mcmc$Btheta[nit,2:S_tilde,2:S_tilde] <- mcmc$tau[nit]^2*diag(mcmc$lambda[nit,]^2)
  
  # Update omega 
  mcmc$omega[nit,,] <- update_omega(mcmc$L[nit], mcmc$R_tilde[nit,], 
                                    mcmc$theta[nit,1:mcmc$L[nit],], mcmc$sigma2[nit-1],
                                    mcmc$sigma2_omega[nit,1:mcmc$L[nit]], mcmc$rho[nit,1:mcmc$L[nit]])
  omega_tilde <- as.vector(mcmc$omega[nit,,])
  
  # Udpate sigma2_omega
  mcmc$sigma2_omega[nit,1:mcmc$L[nit]] <- update_sigma2_omega(mcmc$L[nit], mcmc$R_tilde[nit,], 
                                                              mcmc$omega[nit,,], mcmc$rho[nit,1:mcmc$L[nit]])
  
  # Update rho
  mcmc$rho[nit,1:mcmc$L[nit]] <- update_rho(mcmc$L[nit], mcmc$R_tilde[nit,], mcmc$omega[nit,,], 
                                            mcmc$sigma2_omega[nit,1:mcmc$L[nit]], mcmc$rho[nit,1:mcmc$L[nit]])
  
  # Update sigma2
  mcmc$sigma2[nit] <- update_sigma2(mcmc$L[nit], mcmc$R_tilde[nit,], mcmc$theta[nit,1:mcmc$L[nit],], omega_tilde)
}

duration = proc.time()-start.time
# End of the Chain


# Posterior Samples
post <- NULL
post$tree <- mcmc$tree[post_index] # tree
post$L <- mcmc$L[post_index] # number of terminal nodes
post$R_tilde <- mcmc$R_tilde[post_index,] # index of clusters
post$theta <- mcmc$theta[post_index,1:max(post$L),] # fixed effect coefficients
post$etheta <- mcmc$etheta[post_index,] # prior mean for fixed effect coefficients
post$Btheta <- mcmc$Btheta[post_index,,] # prior covariance for fixed effect coefficients
post$lambda <- mcmc$lambda[post_index,] # local shrinkage parameter
post$tau <- mcmc$tau[post_index] # global shrinkage parameter
post$sigma2_omega <- mcmc$sigma2_omega[post_index,1:max(post$L)] # GP variance
post$rho <- mcmc$rho[post_index,1:max(post$L)] # GP length scale
post$sigma2 <- mcmc$sigma2[post_index] # residual variance
