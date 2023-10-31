#ifndef UTILS_ESTIM_H
#define UTILS_ESTIM_H

arma::vec d1ll(Rcpp::NumericMatrix& Yr,
               Rcpp::CharacterVector& famr,
               ghQ& Q,
               arma::cube& b,
               arma::mat& pD);

arma::mat d2ll_EM(Rcpp::NumericMatrix& Yr,
                  Rcpp::CharacterVector& famr,
                  ghQ& Q,
                  arma::cube& b,
                  arma::mat& pD,
                  const bool flagEIM);

arma::vec solve_EM(Rcpp::NumericMatrix& Yr,
                   Rcpp::CharacterVector& famr,
                   ghQ& Q,
                   arma::cube& b,
                   arma::mat& pD,
                   const bool flagEIM);

arma::mat SigmaGrad(ghQ& Q, fYZ& fyz);

double fZY(ghQ& Q, fYZ& fyz);

void update_sigm(ghQ& Q, fYZ& fyz);

Rcpp::List SEs(Rcpp::NumericMatrix& Yr,
               Rcpp::CharacterVector& famr,
               ghQ& Q,
               arma::cube& b,
               arma::mat& pD,
               bool FLAGCORLV);

void EM_step(Rcpp::NumericMatrix& Yr,
             Rcpp::CharacterVector& famr,
             arma::cube& b,
             ghQ& Q,
             fYZ& fyz,
             Rcpp::List& control);

void DM_step(Rcpp::NumericMatrix& Yr,
             Rcpp::CharacterVector& famr,
             arma::cube& b,
             ghQ& Q,
             fYZ& fyz,
             Rcpp::List& control);

#endif