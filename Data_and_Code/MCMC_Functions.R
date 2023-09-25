################################################## MCMC Functions ##############################################



init <- function(){
  
  # MCMC initialization 
  # initialize parameters need to be estimated in MCMC
  # returns: init_list
  
  L_init <- 1
  R_tilde_init <- rep("0-1", n*max(J))
  R_tilde_init[na_index] <- 0 
  
  tree_init <- Node$new("0-1") # tree root node
  tree_init$index <- R_tilde_init # index of tree terminal nodes
  tree_init$split <- c() # split position of root node 
  tree_init$rule <- c() # split rule of root node 
  tree_init$prior <- 1-a # prior of the tree
  tree_init$node_prior <- 1-a # prior of root node
  tree_init$set <- available_set # available covariate split position 
  
  theta_init <- rmvn_rcpp(L_init, rep(0,S_tilde), diag(1,S_tilde))
  sigma2_omega_init <- 1
  rho_init <- 5
  omega_init <- matrix(NA, nrow=n, ncol=max(J))
  for (i in 1:n){ omega_init[i,1:J[i]] <- rmvn_rcpp(1, rep(0,J[i]), diag(1,J[i])) }
  sigma2_init <- 1
  
  etheta_init <- rep(0, S_tilde) # etheta_init <- rmvn_rcpp(1, e_0, E_0)
  Btheta_init <- diag(1, S_tilde) # Btheta_init <- riwish_rcpp(b_0, B_0)
  lambda_init <- rep(1, S_minus); tau_init <- 1
  nu_init <- rep(1, S_minus); psi_init <- 1
  
  init_list <- list(tree=tree_init, L=L_init, R_tilde=R_tilde_init, theta=theta_init, 
                    omega=omega_init, sigma2_omega=sigma2_omega_init, rho=rho_init, sigma2=sigma2_init,
                    etheta=etheta_init, Btheta=Btheta_init,
                    lambda=lambda_init, tau=tau_init, nu=nu_init, psi=psi_init)
  return(init_list)
}



update_theta <- function(L, R_tilde, sigma2_omega, rho, sigma2, etheta, Btheta){
  
  # update theta
  # args: L: number of tree terminal nodes
  #       R_tilde: index of tree terminal nodes 
  #       sigma2_omega, rho, sigma2: parameters
  #       etheta, Btheta: hyper-parameters
  # returns: theta_update
  
  theta_update <- matrix(NA, nrow=L, ncol=S_tilde)
  
  Btheta_inv <- chol2inv(chol(Btheta))
  Btheta_inv_etheta <- Btheta_inv %*% etheta
  R <- matrix(R_tilde, nrow=n, ncol=max(J))[,1] # index of tree terminal nodes for each individuals
  R_index <- sort(unique(R)) # unique index of tree terminal nodes
  
  for (l in 1:L){
    R_l <- which(R==R_index[l]) # index of individuals in l-th terminal nodes 
    X_tilde_sum <- matrix(0, nrow=S_tilde, ncol=S_tilde)
    Xy_tilde_sum <- rep(0, S_tilde)
    
    for (i in R_l){
      # covariance matrix for i-th individual
      Kappa_i <- Kappa[i,1:J[i],1:J[i]]
      Rho_i <- matrix(NA, nrow=J[i], ncol=J[i])
      for (j in 1:J[i]){ Rho_i[j,1:J[i]] <- exp(-abs(t[i,j]-t[i,1:J[i]])/rho[l]) }
      Sigma_i <- sigma2_omega[l]*Kappa_i*Rho_i + diag(sigma2,J[i])
      
      # calculate X_tilde_sum, Xy_tilde_sum
      X_tilde_Sigma_inv <- t(bX[i,1:J[i],]) %*% chol2inv(chol(Sigma_i))
      X_tilde_sum <- X_tilde_sum + X_tilde_Sigma_inv %*% bX[i,1:J[i],]
      Xy_tilde_sum <- Xy_tilde_sum + X_tilde_Sigma_inv %*% y[i,1:J[i]]
    }
    
    # posterior distribution for theta_l, l = 1,2,...,L
    V_n <- chol2inv(chol(X_tilde_sum + Btheta_inv))
    mu_n <- V_n%*%(Xy_tilde_sum+Btheta_inv_etheta)
    
    # update theta_l, l = 1,2,...,L
    theta_update[l,] <- as.vector(rmvn_rcpp(1, mu_n, V_n))
  }
  
  return(theta_update)
}



update_lambda <- function(L, theta, tau, nu){
  
  # update lambda
  # args: L: number of tree terminal nodes
  #       theta: S_minus dimensional parameters with shrinkage priors
  #       tau, nu: hyper-parameters
  # returns: lambda_update
  
  lambda_update <- rep(NA, S_minus)
  if (length(dim(theta))!=2) { theta <- array(theta, dim=c(1, S_minus)) }
  
  for (s in 1:S_minus){
    a_star <- (L+1)/2
    b_star <- 1/nu[s] + sum(theta[,s]^2)/(2*tau^2)
    lambda_update[s] <- sqrt(rinvgamma(1, a_star, b_star))
  }
  return(lambda_update)
}



update_tau <- function(L, theta, lambda, psi){
  
  # update tau
  # args: L: number of tree terminal nodes
  #       theta: S_minus dimensional parameters with shrinkage priors
  #       lambda, psi: hyper-parameters
  # returns: tau_update
  
  tau_update <- rep(NA, 1)
  if (length(dim(theta))!=2) { theta <- array(theta, dim=c(1, S_minus)) }
  
  theta_star <- 0
  for (s in 1:S_minus){ 
    theta_star <- theta_star + sum((theta[,s]/lambda[s])^2) 
  }
  a_star <- (L*S_minus+1)/2
  b_star <- 1/psi + theta_star/2
  tau_update <- sqrt(rinvgamma(1, a_star, b_star))
  return(tau_update)
}



update_nu <- function(lambda){
  
  # update nu
  # args: lambda: hyper-parameters
  # returns: nu_update
  
  nu_update <- rep(NA, S_minus)
  for (s in 1:S_minus){
    nu_update[s] <- rinvgamma(1, 1, 1+1/lambda[s]^2)
  }
  return(nu_update)
}



update_psi <- function(tau){
  
  # update psi
  # args: tau: hyper-parameters
  # returns: psi_update
  
  psi_update <- rinvgamma(1, 1+1/tau^2)
  return(psi_update)
}



update_omega <- function(L, R_tilde, theta, sigma2, sigma2_omega, rho){
  
  # update omega
  # args: L: number of tree terminal nodes
  #       R_tilde: index of tree terminal nodes 
  #       theta, sigma2: parameters
  #       sigma2_omega, rho: hyper-parameters
  # returns: omega_update
  
  omega_update <- matrix(NA, nrow=n, ncol=max(J))
  
  y_tilde_vec <- rep(NA, n*max(J))
  R_tilde_index <- sort(unique(R_tilde[R_tilde!=0])) # unique index of tree terminal nodes
  if (length(dim(theta))!=2) { theta <- array(theta, dim=c(1, S_tilde)) }
  for (l in 1:L){ 
    R_l <- which(R_tilde==R_tilde_index[l]) # index of data points in l-th terminal node
    y_tilde_vec[R_l] <- y_tilde[R_l] - bX_tilde[R_l,]%*%theta[l,]
  }
  y_tilde_mat <- matrix(y_tilde_vec, nrow=n, ncol=max(J))
  
  R <- matrix(R_tilde, nrow=n, ncol=max(J))[,1]
  for (i in 1:n){
    R_i <- which(R_tilde_index==R[i]) # terminal node index of i-th individual
    # covariance matrix for i-th individual
    Kappa_i <- Kappa[i,1:J[i],1:J[i]]
    Rho_i <- matrix(NA, nrow=J[i], ncol=J[i])
    for (j in 1:J[i]){ Rho_i[j,1:J[i]] <- exp(-abs(t[i,j]-t[i,1:J[i]])/rho[R_i]) }
    # update omega_i, i = 1,2,..,n
    V_n <- chol2inv(chol(diag(1/sigma2,J[i]) + 1/sigma2_omega[R_i]*chol2inv(chol(Kappa_i*Rho_i))))
    mu_n <- V_n %*% y_tilde_mat[i,1:J[i]]/sigma2
    omega_update[i,1:J[i]] <- rmvn_rcpp(1, mu_n, V_n)
  }
  
  return(omega_update)
}



update_sigma2_omega <- function(L, R_tilde, omega, rho){
  
  # update sigma2_omega
  # args: L: number of tree terminal nodes
  #       R_tilde: index of tree terminal nodes 
  #       omega, rho: parameters
  # returns: sigma2_omega_update
  
  sigma2_omega_update <- rep(NA, L)
  R <- matrix(R_tilde, nrow=n, ncol=max(J))[,1]
  R_index <- sort(unique(R)) # unique index of tree terminal nodes
  
  for (l in 1:L){
    omega_sum <- 0
    R_l <- which(R==R_index[l]) # index of individuals in l-th terminal nodes 
    for (i in R_l){
      # covariance matrix for i-th individual
      Kappa_i <- Kappa[i,1:J[i],1:J[i]]
      Rho_i <- matrix(NA, nrow=J[i], ncol=J[i])
      for (j in 1:J[i]){ Rho_i[j,1:J[i]] <- exp(-abs(t[i,j]-t[i,1:J[i]])/rho[l]) }
      omega_sum <- omega_sum + t(omega[i,1:J[i]])%*%chol2inv(chol(Kappa_i*Rho_i))%*%omega[i,1:J[i]]
    }
    # update sigma2_omega
    sigma2_omega_update[l] <- rgig(n=1, lambda=f_1-sum(J[R_l])/2, chi=omega_sum, psi=2*f_2)
  }
  return(sigma2_omega_update)
}



update_rho <- function(L, R_tilde, omega, sigma2_omega, rho){
  
  # update rho
  # args: L: number of tree terminal nodes
  #       R_tilde: index of tree terminal nodes 
  #       omega, sigma2_omega, rho: parameters
  # returns: rho_update
  
  rho_update <- rep(NA, L)

  # step size of proposal distribution
  step <- 0.2 
  # avoid numerical issue
  lower_bound <- 0.1; upper_bound <- 100
  R <- matrix(R_tilde, nrow=n, ncol=max(J))[,1]
  R_index <- sort(unique(R)) # unique index of tree terminal nodes
  
  for (l in 1:L){
    R_l <- which(R==R_index[l]) # index of individuals in l-th terminal nodes 
    lower <- max(rho[l]-step, lower_bound); upper <- min(rho[l]+step, upper_bound)
    rho_new <- runif(1, lower, upper) # uniform random walk
    ratio <- logpost_rho(R_l, omega, sigma2_omega[l], rho_new) - logpost_rho(R_l, omega, sigma2_omega[l], rho[l])
    if (log(runif(1)) < ratio){
      rho_update[l] <- rho_new
    }else { rho_update[l] <- rho[l] }
  }
  return(rho_update)
}



logpost_rho <- function(R_l, omega, sigma2_omega, rho){
  
  # calculate log posterior of rho
  # args: R_l: index of individuals in l-th terminal nodes 
  #       omega, sigma2_omega, rho: parameters
  # returns: logpost
  
  logll <- 0
  for (i in R_l){
    # covariance matrix for i-th individual
    Kappa_i <- Kappa[i,1:J[i],1:J[i]]
    Rho_i <- matrix(NA, nrow=J[i], ncol=J[i])
    for (j in 1:J[i]){ Rho_i[j,1:J[i]] <- exp(-abs(t[i,j]-t[i,1:J[i]])/rho) }
    logll <- logll + dmvn_rcpp(omega[i,1:J[i]], rep(0,J[i]), sigma2_omega*Kappa_i*Rho_i, logd=T)
  }
  
  # logpost <- logll + dbeta(rho, a_rho, b_rho, log=T)
  logpost <- logll + log(dinvgamma(rho, a_rho, b_rho))
  return(logpost)
}



update_sigma2 <- function(L, R_tilde, theta, omega_tilde){
  
  # update sigma2
  # args: L: number of tree terminal nodes
  #       R_tilde: index of tree terminal nodes 
  #       theta, omega_tilde: parameters
  # returns: sigma2_update
  
  y_tilde_sum <- 0
  if (length(dim(theta))!=2) { theta <- array(theta, dim=c(1, S_tilde)) }
  R_tilde_index <- sort(unique(R_tilde[R_tilde!=0])) # unique index of tree terminal nodes
  for (l in 1:L){ 
    R_l <- which(R_tilde==R_tilde_index[l]) # index of data points in l-th terminal node
    y_tilde_sum <- y_tilde_sum + sum((y_tilde[R_l] - bX_tilde[R_l,]%*%theta[l,] - omega_tilde[R_l])^2) 
  }
  
  # update sigma2
  g_1_star <- g_1 + sum(J)/2
  g_2_star <- g_2 + y_tilde_sum/2
  sigma2_udpate <- rinvgamma(1, g_1_star, g_2_star)
  return(sigma2_udpate)
}