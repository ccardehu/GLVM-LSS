#ifndef UTILS_BASIC_H
#define UTILS_BASIC_H

arma::uword ixd2ll(int& k1, int& k2, int& K);

arma::mat count_b(arma::cube& b, const int p);

arma::cube lb2Cb(Rcpp::List& b, const int p, const int q);

arma::mat Cb2Mb(arma::cube& b);

arma::vec Cb2Vb(arma::cube& b);

arma::cube Mb2Cb(arma::mat mb, const int q, const int K);

arma::cube Vb2Cb(arma::vec& vb, arma::cube& b, const int q);

arma::mat vec2mat(arma::vec& x, int nrow, int ncol);

void fixY(arma::vec& y);

arma::vec expitvec(arma::vec eta);

arma::vec logitvec(arma::vec mu);

arma::vec dMd1Evec(arma::vec eta);

arma::vec dMd2Evec(arma::vec eta);

arma::vec geta(arma::vec& mu, arma::vec&sg);

arma::vec getb(arma::vec& mu, arma::vec&sg);

arma::vec zeta1(arma::vec& y);

arma::vec zeta2(arma::vec& y);

arma::mat d1mu(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat d1sg(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat d1nu(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat d2mu(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat d2sg(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat d2nu(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat dcms(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat dcmn(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat dcsn(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat d2muF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat d2sgF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat d2nuF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat dcmsF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat dcmnF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

arma::mat dcsnF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu);

double dsn(double& x, double& mu, double& sigma, double& nu);

double rsn(double& mu, double& sigma, double& nu);

arma::vec nc2nd(arma::vec& nuc);

arma::vec sc2sd(arma::vec& sgc, arma::vec& nuc);

arma::vec mc2md(arma::vec& muc, arma::vec& sgc, arma::vec& nuc);

arma::vec p2nuc(arma::vec& p);

arma::vec dsddsc(arma::vec& nuc);

arma::vec dnddnc(arma::vec& nuc);

arma::vec dmddsc(arma::vec& nuc);

arma::vec dmddnc(arma::vec& sgc, arma::vec& nuc);

arma::vec dsddnc(arma::vec& sgc, arma::vec& nuc);

arma::mat Ezm(int q,
              arma::mat& Qpi,
              arma::vec& Qw,
              arma::mat& pD);

arma::cube Vzm(int q,
               arma::mat& Qpi,
               arma::vec& Qw,
               arma::mat& pD,
               arma::mat& EZM);

double EVzm(arma::mat& G,
            arma::mat& EZM,
            arma::cube& VZM);

arma::vec log_mvnN(arma::mat& y, arma::mat& S);

void fixL(arma::mat& L);

Rcpp::NumericVector chol2cor(Rcpp::NumericVector& vLi, int& q);

Rcpp::NumericMatrix Jchol(Rcpp::NumericVector& vLi,int q);

#endif
