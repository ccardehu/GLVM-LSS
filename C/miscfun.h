#ifndef MISCFUN_H
#define MISCFUN_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

arma::mat DFtoNM(Rcpp::DataFrame x);

Rcpp::NumericVector cprobs (arma::mat x);

Rcpp::NumericVector dZIPo (NumericVector Y,
                           NumericVector mu,
                           NumericVector sg,
                           bool rlog = true);

#endif //MISCFUN_H