############################################### Simulated Data #################################################

# load the treatment history data
n <- data$n # number of individuals
J <- data$J # number of visits
kappa <- data$kappa # ST kernel similarity matrix
index_regimens <- data$index_regimens # index of unique cART regimen

# start to generate simulated data 
set.seed(123) 

S <- 10 # dimension of covariates 
C_1 <- 4 # number of categories for the first categorical covariate
C_2 <- 3 # number of categories for the second categorical covariate
X <- array(NA, dim=c(n, max(J), S)) # the covaraite matrix
Xcat_1 <- matrix(NA, nrow=n, ncol=max(J)) # time-invariant categorical covariate
Xcat_2 <- matrix(NA, nrow=n, ncol=max(J)) # time-invariant categorical covariate

# generate the baseline covariates
for (i in 1:n){
  X[i,1:J[i],1] <- rep(1, J[i]) # intercept term
  X[i,1:J[i],2] <- sample(0:1, J[i], replace=T) # time-varying binary covariate
  X[i,1:J[i],3] <- rnorm(J[i], 0, 1) # time-varying continuous covariate
  # first time-invariant categorical covariate
  Xcat_1[i,1:J[i]] <- rep(sample(0:(C_1-1), 1, replace=F), J[i])
  if (Xcat_1[i,1] == 1){
    X[i,1:J[i],4] <- rep(1, J[i])
    X[i,1:J[i],5] <- rep(0, J[i])
    X[i,1:J[i],6] <- rep(0, J[i])
  }else if (Xcat_1[i,1] == 2){
    X[i,1:J[i],4] <- rep(0, J[i])
    X[i,1:J[i],5] <- rep(1, J[i])
    X[i,1:J[i],6] <- rep(0, J[i])
  }else if (Xcat_1[i,1] == 3){
    X[i,1:J[i],4] <- rep(0, J[i])
    X[i,1:J[i],5] <- rep(0, J[i])
    X[i,1:J[i],6] <- rep(1, J[i])
  }else { 
    X[i,1:J[i],4] <- rep(0, J[i])
    X[i,1:J[i],5] <- rep(0, J[i])
    X[i,1:J[i],6] <- rep(0, J[i])
  }
  X[i,1:J[i],7] <- sample(0:1, J[i], replace=T) # time-varying binary covariate
  X[i,1:J[i],8] <- rnorm(J[i], 0, 1) # time-varying continuous covariate
  # second time-invariant categorical covariate
  Xcat_2[i,1:J[i]] <- rep(sample(0:(C_2-1), 1, replace=F), J[i])
  if (Xcat_2[i,1] == 1){
    X[i,1:J[i],9] <- rep(1, J[i])
    X[i,1:J[i],10] <- rep(0, J[i])
  }else if (Xcat_2[i,1] == 2){
    X[i,1:J[i],9] <- rep(0, J[i])
    X[i,1:J[i],10] <- rep(1, J[i])
  }else {
    X[i,1:J[i],9] <- rep(0, J[i])
    X[i,1:J[i],10] <- rep(0, J[i])
  }
}

# transfer the range of continuous variable to be between 0 and 1
x3min <- range(X[,,3][!is.na(X[,,3])])[1] 
x3max <- range(X[,,3][!is.na(X[,,3])])[2] 
for (i in 1:n){ X[i,1:J[i],3] <- (X[i,1:J[i],3]-x3min)/(x3max-x3min) }
x8min <- range(X[,,8][!is.na(X[,,8])])[1] 
x8max <- range(X[,,8][!is.na(X[,,8])])[2] 
for (i in 1:n){ X[i,1:J[i],8] <- (X[i,1:J[i],8]-x8min)/(x8max-x8min) }

M <- 3 # dimension of genetics variable 
G <- array(NA, dim=c(n, max(J), M))
for (i in 1:n){
  G_ij <- sample(0:1, M, replace=T)
  G[i,1:J[i],1:M] <- matrix(G_ij, nrow=J[i], ncol=M, byrow=T)
}

D <- 4 # number of representative ART regimens
index_kernel <- c(47, 12, 66, 149) # selected representative cART regimens for kernel regression
H <- array(NA, dim=c(n, max(J), D)) # drug kernel regression design matrix
for (i in 1:n){ 
  for (j in 1:J[i]){
    if (sum(kappa[index_regimens[i,j],index_kernel])!=0){
      H[i,j,1:D] <- kappa[index_regimens[i,j],index_kernel] / sum(kappa[index_regimens[i,j],index_kernel]) 
    }else{
      H[i,j,1:D] <- 0 
    }
  }
}

Q <- array(NA, dim=c(n, max(J), M*D)) # drug-genetics interaction design matrix         
for (i in 1:n){
  for (j in 1:J[i]){
    Q[i,j,] <- H[i,j,] %x% G[i,j,]
  }
}

# normalize the design matrix 
H_mean <- rep(NA, D); H_sd <- rep(NA, D)
for (d in 1:D){ 
  H_hat <- H[,,d][!is.na(H[,,d])]
  H[,,d] <- (H[,,d] - mean(H_hat)) / sd(H_hat) 
  H_mean[d] <- mean(H_hat); H_sd[d] <- sd(H_hat)
}

Q_mean <- rep(NA, M*D); Q_sd <- rep(NA, M*D)
for (q in 1:(M*D)){ 
  Q_hat <- Q[,,q][!is.na(Q[,,q])]
  Q[,,q] <- (Q[,,q] - mean(Q_hat)) / sd(Q_hat) 
  Q_mean[q] <- mean(Q_hat); Q_sd[q] <- sd(Q_hat)
}

# regression tree
L <- 4 # number of tree terminal nodes 
# L <- 1 # number of tree terminal nodes 
R <- array(NA, dim=c(n, max(J))) # name index of tree terminal node for each (i,j) data points
Rindex <- array(NA, dim=c(n, max(J))) # number index of tree terminal node for each (i,j) data points
for (i in 1:n){ 
  # tree node name has the form "depth-index"
  # index is from 1 to 2^{depth}, e.g., the root node is "0-1"
  for (j in 1:J[i]){
    if (sum(X[i,1,4:6]) == 0 | X[i,1,4] == 1){
      if (X[i,1,2] == 0){
        R[i,j] <- "2-1"; Rindex[i,j] <- 1
      }else { R[i,j] <- "2-2"; Rindex[i,j] <- 2 }
    }else {
      if (X[i,1,3] <= 0.5){
        R[i,j] <- "2-3"; Rindex[i,j] <- 3
      }else { R[i,j] <- "2-4"; Rindex[i,j] <- 4 }
    }
    # R[i,j] <- "0-1"; Rindex[i,j] <- 1
  }
}

beta <- matrix(0, nrow=L, ncol=S)
gamma <- matrix(0, nrow=L, ncol=D)
delta <- matrix(0, nrow=L, ncol=M)
alpha <- matrix(0, nrow=L, ncol=M*D)

beta[1,1:3] <- c(1, -1, 1.5)
beta[2,1:3] <- c(-1, 1, -1.5)
beta[3,1:3] <- c(1.5, 1.5, -1)
beta[4,1:3] <- c(1.5, 1.5, 1)
gamma[1:4,1] <- c(1.5, -1.5, 1, -1)
delta[1:4,1] <- c(1, 1.5, 2, 2.5)
alpha[1,1:2] <- c(1, -1)
alpha[2,1:2] <- c(-1.5, -1.5)
alpha[3,1:2] <- c(-1, 1.5)
alpha[4,1:2] <- c(1, 1)

# gaussian process 
omega <- array(NA, dim=c(n, max(J)))
sigma2_omega <- c(0.5, 0.75, 1, 1.25)
rho <- c(4, 6, 8, 10)
K <- array(NA, dim=c(n, max(J), max(J))) # kernel matrix for Gaussian process
Kappa <- array(NA, dim=c(n, max(J), max(J))) # drug similarity matrix for Gaussian process
Rho <- array(NA, dim=c(n, max(J), max(J))) # time discount matrix for Gaussian process
t <- array(NA, dim=c(n, max(J))) 
for (i in 1:n){
  t[i,1:J[i]] <- sort(rnorm(J[i], 0, 1))
  R_i <- Rindex[i,1] # terminal node index of i-th individual
  for (j in 1:J[i]){
    Rho[i,j,1:J[i]] <- exp(-abs(t[i,j]-t[i,1:J[i]])/rho[R_i])
    Kappa[i,j,1:J[i]] <- kappa[index_regimens[i,j],index_regimens[i,1:J[i]]] / 
      sqrt(kappa[index_regimens[i,j],index_regimens[i,j]] * diag(kappa[index_regimens[i,1:J[i]],index_regimens[i,1:J[i]]]))
  }
  K[i,1:J[i],1:J[i]] <- sigma2_omega[R_i] * Kappa[i,1:J[i],1:J[i]] * Rho[i,1:J[i],1:J[i]] 
  omega[i,1:J[i]] <- rmvn_rcpp(1, rep(0,J[i]), K[i,1:J[i],1:J[i]])
}


# simulated continuous outcomes
sigma2 <- 1 # variance of independent Gaussian noise
y <- array(NA, dim=c(n, max(J))) 
for (i in 1:n){
  for (j in 1:J[i]){
    y[i,j] <- rnorm(1, X[i,j,]%*%beta[Rindex[i,j],] + H[i,j,]%*%gamma[Rindex[i,j],] + 
        G[i,j,]%*%delta[Rindex[i,j],] + Q[i,j,]%*%alpha[Rindex[i,j],] + omega[i,j], sigma2)
  }
}

# prior specification: hyper-parameters
S_tilde <- S+D+M+M*D # dimension of all the parameters
S_minus <- S_tilde - 1 # dimension of shrinkage parameters 
theta <- cbind(beta, gamma, delta, alpha)
etheta <- rep(0, S_tilde); Btheta <- diag(1, S_tilde)
Btheta_inv <- chol2inv(chol(Btheta)); Btheta_inv_etheta <- Btheta_inv %*% etheta
f_1 <- 1; f_2 <- 10 # sigma2_omega \sim gamma(f_1, f_2)
g_1 <- 1; g_2 <- 1 # sigma2 \sim inverse-gamma(g_1, g_2)
a_rho <- 1; b_rho <- 1 # rho \sim inverse-gamma(a_rho, b_rho)
a <- 0.95; b <- 1 # tree split prior parameters
xi <- 0.2 # power for the fractional bayes factor

# pre-calculated data for saving running time  
R_tilde <- as.vector(R); R_tilde[is.na(R_tilde)] <- 0
Rindex_tilde <- as.vector(Rindex); Rindex_tilde[is.na(Rindex_tilde)] <- 0
y_tilde <- as.vector(y); y_tilde[is.na(y_tilde)] <- 0
omega_tilde <- as.vector(omega); omega_tilde[is.na(omega_tilde)] <- 0
X_tilde <- matrix(X, nrow=n*max(J), ncol=S); X_tilde[is.na(X_tilde)] <- 0
H_tilde <- matrix(H, nrow=n*max(J), ncol=D); H_tilde[is.na(H_tilde)] <- 0
G_tilde <- matrix(G, nrow=n*max(J), ncol=M); G_tilde[is.na(G_tilde)] <- 0
Q_tilde <- matrix(Q, nrow=n*max(J), ncol=M*D); Q_tilde[is.na(Q_tilde)] <- 0
bX_tilde <- cbind(X_tilde, H_tilde, G_tilde, Q_tilde)
bX_tilde_sum <- array(0, dim=c(n*max(J), S_tilde, S_tilde))
for (i in 1:(n*max(J))){ bX_tilde_sum[i,,] <- bX_tilde[i,] %*% t(bX_tilde[i,]) }
bX_tilde_trans <- t(bX_tilde)
bX_tilde_sum_trans <- aperm(bX_tilde_sum, c(2,3,1))
bX <- array(NA, dim=c(n, max(J), S_tilde))
for (i in 1:n){ for (j in 1:J[i]){ bX[i,j,] <- c(X[i,j,], H[i,j,], G[i,j,], Q[i,j,]) }}

# tree related pre-calculated data 
depth_max <- 3 # maximum tree depth 
L_max <- 2^depth_max # maximum number of tree terminal nodes
Xcat_1_tilde <- matrix(Xcat_1, nrow=n*max(J), ncol=1); Xcat_1_tilde[is.na(Xcat_1_tilde)] <- 0
Xcat_2_tilde <- matrix(Xcat_2, nrow=n*max(J), ncol=1); Xcat_2_tilde[is.na(Xcat_2_tilde)] <- 0
available_set <- c(1, 2, 3, 4, 5, 6) # available covariate split position 
data_type <- c("binary", "continuous", "ordinal", "binary", "continuous", "categorical") # data type of the covariate
na_index <- which(R_tilde==0); data_index <- which(R_tilde!=0) # index of NA and non-NA data
min_num <- 10 # minimum number of data points in each terminal nodes to split further 
S_tree <- length(available_set) 
X_tree_array <- array(NA, dim=c(n, max(J), S_tree))
for (i in 1:n){ for (j in 1:J[i]){ X_tree_array[i,j,1:S_tree] <- c(X[i,1,2:3], Xcat_1[i,1], X[i,1,7:8], Xcat_2[i,1]) }}
Xtree <- matrix(X_tree_array, nrow=n*max(J), ncol=S_tree) # available split covariates for regression tree

# tree simulation truth in "data.tree" type 
tree <- Node$new("0-1") # tree root node
tree$index <- R_tilde # index of tree terminal nodes
tree$split <- c(3) # split position of root node 
tree$rule <- c(1) # split rule of root node 
tree$prior <- a*((a*(1+1)^(-b)))^2*(1-(a*(1+2)^(-b)))^4 # prior of the tree
tree$node_prior <- a # prior of root node
tree$set <- available_set # available covariate split position 
tree$AddChild("1-1"); tree$`1-1`$node_prior <- a*(1+1)^(-b)
tree$`1-1`$split <- c(1)
tree$`1-1`$rule <- c(0)
tree$`1-1`$set <- available_set[available_set!=1] 
tree$`1-1`$AddChild("2-1"); tree$`1-1`$`2-1`$node_prior <- 1-(a*(1+2)^(-b))
tree$`1-1`$`2-1`$set <- available_set[available_set!=1] 
tree$`1-1`$AddChild("2-2"); tree$`1-1`$`2-2`$node_prior <- 1-(a*(1+2)^(-b))
tree$`1-1`$`2-2`$set <- available_set[available_set!=1] 
tree$AddChild("1-2"); tree$`1-2`$node_prior <- a*(1+1)^(-b)
tree$`1-2`$split <- c(2)
tree$`1-2`$rule <- c(0.5)
tree$`1-2`$set <- available_set
tree$`1-2`$AddChild("2-3"); tree$`1-2`$`2-3`$node_prior <- 1-(a*(1+2)^(-b))
tree$`1-2`$`2-3`$set <- available_set
tree$`1-2`$AddChild("2-4"); tree$`1-2`$`2-4`$node_prior <- 1-(a*(1+2)^(-b))
tree$`1-2`$`2-4`$set <- available_set
print(tree, "split", "rule", "prior", "node_prior", "set")