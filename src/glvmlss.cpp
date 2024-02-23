#include <RcppArmadillo.h>
#include <RcppEnsmallen.h>
#include "ghq_class.h"
#include "dist_class.h"
#include "fyz_class.h"
#include "utils_basic.h"
#include "utils_estimation.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEnsmallen)]]

using namespace Rcpp;
using namespace arma;
using namespace ens;

// [[Rcpp::export]]
Rcpp::List glvmlss_rcpp(Rcpp::NumericMatrix& Yr,
                        Rcpp::CharacterVector& famr,
                        const int q,
                        arma::cube& b,
                        Rcpp::List& control){

  int nqp = control["nQP"];
  bool flagVERBO = control["verbose"], flagCORLV = control["corLV"] ;
  bool flagGD = control["EM.useGD"];
  arma::mat Si = control["start.R"];

  ghQ Q(q,nqp,Si);
  fYZ fyz(Yr,famr,Q,b);

  int EM_iL = control["EM.maxit"];
  int DM_iL = control["DM.maxit"];

  if(EM_iL != 0){
    if(flagGD){
      EM_stepGD(Yr,famr,b,Q,fyz,control);
    } else {
      EM_step(Yr,famr,b,Q,fyz,control);
    }
  }
  if(DM_iL != 0) DM_step(Yr,famr,b,Q,fyz,control);

  arma::mat pD = fyz.get_pD();
  if(flagVERBO) Rcpp::Rcout << "\n Computing standard errors ... " ;
  Rcpp::List se = SEs(Yr,famr,Q,b,pD,flagCORLV);
  if(flagVERBO) Rcpp::Rcout << "\r Computing standard errors ... Done!" ;

  arma::cube seb = se["seb"];

  arma::mat fsc = fscores(Q,pD);
  int n = Yr.nrow() ;
  int K = se["K"];
  double AIC = -2.0*fyz.get_ll() + 2.0*K;
  double BIC = -2.0*fyz.get_ll() + std::log(n)*K;

  if(flagCORLV){
    arma::vec seR = se["seR"];
    return Rcpp::List::create(Rcpp::Named("b") = b,
                              Rcpp::Named("R") = Q.get_sigm(),
                              Rcpp::Named("ll") = fyz.get_ll(),
                              Rcpp::Named("seb") = seb,
                              Rcpp::Named("seR") = seR,
                              Rcpp::Named("z") = fsc,
                              Rcpp::Named("AIC") = AIC,
                              Rcpp::Named("BIC") = BIC,
                              Rcpp::Named("K") = K) ;
  } else{
    return Rcpp::List::create(Rcpp::Named("b") = b,
                              Rcpp::Named("R") = Q.get_sigm(),
                              Rcpp::Named("ll") = fyz.get_ll(),
                              Rcpp::Named("seb") = seb,
                              Rcpp::Named("z") = fsc,
                              Rcpp::Named("AIC") = AIC,
                              Rcpp::Named("BIC") = BIC,
                              Rcpp::Named("K") = K) ;
  }
}

// [[Rcpp::export]]
Rcpp::List glvmlss_rcpp_sim(const int n,
                            Rcpp::CharacterVector& famr,
                            const int q,
                            arma::cube& b,
                            Rcpp::List& control){

  int nqp = control["nQP"];
  arma::mat Si = control["start.R"];

  int p = famr.length();
  arma::mat Yo(n,p);

  ghQ Q(q,nqp,Si);
  Rcpp::NumericMatrix Zr = Q.get_simZ(q,n);
  arma::mat Z(Zr.begin(),n,q+1,false,false);

  for(int i = 0; i < p; i++){
    std::string fam = Rcpp::as<std::string>(famr(i));
    Dist Y(fam) ;
    arma::vec Yi = Y.sim(n,i,b,Z);
    Yo.col(i) = Yi;
  }

  return Rcpp::List::create(Rcpp::Named("Y") = Yo,
                            Rcpp::Named("Z") = Z,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("R") = Si) ;
}
