//[[Rcpp::depends(RcppArmadillo)]]
# include <RcppArmadillo.h>


using namespace std;
using namespace arma;
using namespace Rcpp; 


const double log2pi = log(2.0*M_PI);



// [[Rcpp::export]]
double dmvn_rcpp(arma::rowvec& x, arma::rowvec& mean, arma::mat& sigma, bool logd = false){ 
  
  // calculate density of multivariate normal distribution
  //
  // args: x: row vector data
  //      mean: row vector mean, sigma: covaraicne matrix  
  //      logd: true for taking log
  // returns: out: pdf (or log pdf) of multivariate normal distribution
  
  int xdim = x.size(); 
  arma::mat rooti = trans(inv(trimatu(chol(sigma))));
  double rootisum = sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
  
  arma::vec z = rooti*trans(x-mean);
  double out = constants-0.5*sum(z%z)+rootisum;
  
  if (logd == false){ out = exp(out); }
  return(out);
}



// [[Rcpp::export]]
arma::mat rmvn_rcpp(const int n, arma::vec& mean, arma::mat& sigma){
  
  // randomly generate samples from multivariate normal distribution
  //
  // args: n: number of data 
  //      mean: row vector mean, sigma: covaraicne matrix  
  // returns: out: random samples from multivariate normal distribution
  
  int k = sigma.n_cols; // dimension of the multivariate normal distribution
  arma::mat z = randn(n, k);
  arma::mat out = repmat(mean,1,n).t()+z*chol(sigma);
  return(out);
}



// [[Rcpp::export]]
arma::mat riwish_rcpp(const int df, arma::mat& S){
  
  // randomly generate matrix from inverse wishart distribution
  //
  // args: df: degrees of freedom
  //       S: inverse of scale matrix 
  // returns: out: random matrix from inverse wishart distribution
  
  S = S.i(); 
  int m = S.n_rows;
  
  arma::mat Z(m,m); // Bartlett decomposition
  for (int i = 0; i < m; i++){
    Z(i,i) = sqrt(R::rchisq(df-i)); // Fill the diagonal
  }
  for (int j = 0; j < m; j++){  
    for(int i = j+1; i < m; i++){    
      Z(i,j) = R::rnorm(0,1); // Fill the lower matrix 
    }
  }
  
  arma::mat C = trimatl(Z).t() * chol(S);
  // Random matrix from inverse wishart distribution
  arma::mat out = (C.t()*C).i();
  return(out);
}



// [[Rcpp::export]]
double logll_tree_rcpp(arma::vec& R_num, arma::vec& sigma2_omega, arma::vec& rho,
                       const double sigma2, arma::vec& etheta, arma::mat& Btheta,
                       const int L, const int S_tilde, arma::vec& J, arma::cube& Kappa, 
                       arma::mat& t, arma::cube& bX, arma::mat& y, double frac = 1.0){
  
  // calculate the log (power) integrated likelihood of the tree 
  // args: R_num: index of tree terminal nodes 
  //       L: number of tree terminal nodes
  //       sigma2_omega, rho, sigma2: parameters
  //       etheta, Btheta: hyper-parameters
  //       S_tilde, J, Kappa, t, bX, y: data
  //       frac: power of the fractional bayes factor
  // returns: logll
  
  double logll = 0;
  arma::mat Btheta_inv = inv_sympd(Btheta);
  arma::vec Btheta_inv_etheta = Btheta_inv * etheta;
  arma::uvec R_l;
  double y_tilde_sum, Sigma_det_logsum;
  arma::mat X_tilde_sum(S_tilde, S_tilde);
  arma::vec Xy_tilde_sum(S_tilde);
  double sigma2_omega_l, rho_l;
  arma::vec e_n(S_tilde);
  arma::mat B_n(S_tilde, S_tilde), B_n_inv(S_tilde, S_tilde);
  
  // tree log integrated likelihood
  for (int l=0; l<L; l++){
    R_l = find(R_num == (l+1)); // index of individuals in l-th terminal nodes 
    sigma2_omega_l = sigma2_omega(l); rho_l = rho(l); 
    X_tilde_sum.fill(0); Xy_tilde_sum.fill(0);
    y_tilde_sum = 0; Sigma_det_logsum = 0; 
    
    for (int i=0; i<R_l.size(); i++){
      
      // covariance matrix for i-th individual
      int index_i = R_l(i); int J_i = J(index_i);
      arma::mat Kappa_i(J_i, J_i), Rho_i(J_i, J_i);
      Kappa_i = Kappa.subcube(index_i, 0, 0, index_i, J_i-1, J_i-1);
      for (int j=0; j<J_i; j++){
        for (int q=0; q<J_i; q++){
          Rho_i(j,q) = exp(-abs(t(index_i,j) - t(index_i,q))/rho_l);
        }
      }
      arma::mat Sigma_i(J_i, J_i, fill::eye); Sigma_i *= sigma2;
      Sigma_i += sigma2_omega_l * Kappa_i % Rho_i;
      arma::mat Sigma_i_inv = inv_sympd(Sigma_i);
      
      // calculate X_tilde_sum, Xy_tilde_sum, and y_tilde_sum
      arma::mat bX_i(J_i, S_tilde);
      bX_i = bX.subcube(index_i, 0, 0, index_i, J_i-1, S_tilde-1);
      arma::vec y_i(J_i);
      y_i = arma::conv_to< arma::vec >::from(y.submat(index_i, 0, index_i, J_i-1));
      arma::mat X_tilde_Sigma_inv = bX_i.t() * Sigma_i_inv; 
      X_tilde_sum += X_tilde_Sigma_inv * bX_i;
      Xy_tilde_sum += X_tilde_Sigma_inv * y_i;
      y_tilde_sum += arma::conv_to< double >::from(y_i.t() * Sigma_i_inv * y_i);
      
      // calculate Sigma_det_logsum
      double logdet_Sigma_i, sign_Sigma_i;
      log_det(logdet_Sigma_i, sign_Sigma_i, Sigma_i);
      Sigma_det_logsum += logdet_Sigma_i;
    }
    
    // calculate the tree marginal likelihood 
    B_n_inv =  X_tilde_sum + Btheta_inv;
    B_n = inv_sympd(B_n_inv);
    e_n = B_n * (Xy_tilde_sum + Btheta_inv_etheta);
    double logdet_Btheta, sign_Btheta, logdet_B_n, sign_B_n;
    log_det(logdet_Btheta, sign_Btheta, Btheta);
    log_det(logdet_B_n, sign_B_n, B_n);
    logll += frac / 2 * (- sum(J.elem(R_l))*log2pi - Sigma_det_logsum - logdet_Btheta + logdet_B_n - y_tilde_sum - 
      arma::conv_to< double >::from(etheta.t()*Btheta_inv_etheta) + arma::conv_to< double >::from(e_n.t()*B_n_inv*e_n));
  }
  
  return(logll);
}