#ifndef MISCFUN_H
#define MISCFUN_H

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;

inline
arma::mat DFtoNM(Rcpp::DataFrame x) {
  int nRows = x.nrows();
  arma::mat y(nRows, x.size());
  for (int i=0; i<x.size(); i++){
    y.col(i)= as<arma::vec>(x[i]);
  }
  return y;
}

inline
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

inline
double dZIPo (double Y,
              double mu,
              double sg,
              bool rlog = true) {
 int u = 0;
 double lf;
 if(Y == 0) u = 1;
 lf = u*std::log(sg + (1-sg)*std::exp(-mu)) + (1-u)*(std::log(1-sg) - mu + Y*std::log(mu) - std::lgamma(Y + 1)) ;
 if(rlog == true) return lf ;
 else { return std::exp(lf); }
}

#endif //MISCFUN_H