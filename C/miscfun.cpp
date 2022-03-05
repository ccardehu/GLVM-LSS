// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <miscfun.h>

using namespace Rcpp;
using namespace arma;

arma::mat DFtoNM(Rcpp::DataFrame x) {
  int nRows = x.nrows();
  arma::mat y(nRows, x.size());
  for (int i=0; i<x.size(); i++){
    y.col(i)= as<arma::vec>(x[i]);
  }
  return y;
}

Rcpp::NumericVector cprobs (arma::mat x) {
 const double eps = pow(2.220446e-16,0.5) ; 
 int p = x.size() ; 
 NumericVector res(p); 
 for(int i = 0; i < p; i ++){
  res(i) = R::plogis(x(i), 0.0, 1.0, true, false) ;
  if(res(i) == 1) { res(i) = 1 - eps ; } 
  if(res(i) == 0) { res(i) = eps ; } 
 }
 return(res);
}

NumericVector dZIPo (NumericVector Y,
                     NumericVector mu,
                     NumericVector sg,
                     bool rlog = true) {
 NumericVector u(Y.size(), 0.0);
 NumericVector lf(Y.size());
 for(int i = 0; i < Y.size(); i++){ if(Y[i] == 0) u[i] = 1; }
 for(int j = 0; j < Y.size(); j++){
  lf[j] = u[j]*std::log(sg[j] + (1-sg[j])*exp(-mu[j])) + (1-u[j])*(log(1-sg[j]) - mu[j] + Y[j]*log(mu[j]) - std::lgamma(Y[j] + 1)) ;
 }
 if(rlog == true) return lf ;
 else { return exp(lf); }
}
