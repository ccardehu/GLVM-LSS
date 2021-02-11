#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector prb (NumericVector x) {
  
  const double eps = pow(2.220446e-16,2) ; 
  NumericVector res = Rcpp::plogis(x, 0.0, 1.0, true, false) ;
  
  for(int i = 0; i < x.size(); i++){
   if(res[i] == 1) {
    res[i] = 1 - eps;
   }
   if(res[i] == 0) {
    res[i] = eps;
   } 
  }
  return(res);
}

// [[Rcpp::export]]
NumericVector dZIPo (NumericVector Y,
                     NumericVector mu,
                     NumericVector sg,
                     bool rlog = true) {
  
  NumericVector u(Y.size());
  NumericVector lf(Y.size()); 
  
  u.fill(0) ;
  for(int i = 0; i < Y.size(); i++){
    if(Y[i] == 0) {
     u[i] = 1; 
    }
  }
  for(int j = 0; j < Y.size(); j++){
    lf[j] = u[j]*std::log(sg[j] + (1-sg[j])*exp(-mu[j])) + (1-u[j])*(log(1-sg[j]) - mu[j] + Y[j]*log(mu[j]) - std::lgamma(Y[j] + 1)) ;
  }
  if(rlog == true) return lf ;
  else  {
   return exp(lf);
  } 

}

// [[Rcpp::export]]
arma::mat Zreg(DataFrame df,
               Formula formula) {
   Rcpp::Environment stats_env("package:stats");
   Rcpp::Function model_matrix = stats_env["model.matrix"];
   arma::mat df_new = as<arma::mat>(model_matrix(Rcpp::_["object"] = formula, Rcpp::_["data"] = df));
   return(df_new);
}

//// [[Rcpp::export]]
 
// arma::mat dcF (arma::mat Y,
//                Rcpp::List Z,
//                Rcpp::List b,
//                CharacterVector fam) {
//    // Inputs:
//    
//    int n = Y.n_rows ;
//    int p = Y.n_cols ;
//    
//    // Containers:
//    
//    arma::mat R(n,p) ; 
//    R.fill(0)
//      
//    // Algorithm:
//    
//    for(int i = 0; n < p; i++){
//    
//    if(fam[i] == "normal"){
//      arma::vec Zmu = as<arma::mat>(Z["mu"]) ;
//      arma::vec Zsg = as<arma::mat>(Z["sigma"]) ;
//      arma::vec bmu = as<arma::mat>(b["mu"]) ;
//      arma::vec bsg = as<arma::mat>(b["sigma"]) ;
//      arma::vec mu = Z * bmu.row(i) ;
//      arma::vec sg = Z * bsg.row(i) ;
//      R.col(i) = R::dnorm(Y.col(i), mu, sg, log = true) ;
//    }
//    
//    if(fam[i] == "poisson"){
//      arma::vec Zmu = as<arma::mat>(Z["mu"]) ;
//      arma::vec bmu = as<arma::mat>(b["mu"]) ;
//      arma::vec mu = exp(Z * bmu.row(i)) ;
//      R.col(i) = R::dpois(Y.col(i), mu, log = true) ;
//    }
//    
//    if(fam[i] == "gamma"){
//      arma::vec Zmu = as<arma::mat>(Z["mu"]) ;
//      arma::vec Zsg = as<arma::mat>(Z["sigma"]) ;
//      arma::vec bmu = as<arma::mat>(b["mu"]) ;
//      arma::vec bsg = as<arma::mat>(b["sigma"]) ;
//      arma::vec mu = exp(Z * bmu.row(i)) ;
//      arma::vec sg = exp(Z * bsg.row(i)) ;
//      R.col(i) = R::dgamma(Y.col(i), shape = mu, scale = sg, log = true) ;
//    }
//    
//    if(fam[i] == "binom"){
//      arma::vec Zmu = as<arma::mat>(Z["mu"]) ;
//      arma::vec bmu = as<arma::mat>(b["mu"]) ;
//      arma::vec mu = prb(Z * bmu.row(i)) ;
//      R.col(i) = R::dbinom(Y.col(i), 1, mu, log = true) ;
//    }
//      
//    }
//    
//      
//    // Returns:
//    
// } // Rcpp::sourceCpp("C:/Users/carde/Dropbox/Camilo and Irini/Research/GitHub/SPLVM/C")


