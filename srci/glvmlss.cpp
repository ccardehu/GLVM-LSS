#include <RcppArmadillo.h>
#include <RcppEnsmallen.h>
#include "ghq_class.h"
#include "dists.h"
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

  arma::mat Si = arma::eye(q,q);
  if(flagCORLV){
    Si.fill(0.05);
    Si.diag().fill(1.0);
  }
  
  ghQ Q(q,nqp,Si);
  fYZ fyz(Yr,famr,Q,b);
  
  int EM_iL = control["EM.iter.lim"];
  int DM_iL = control["DM.iter.lim"];
  
  if(EM_iL != 0) EM_step(Yr,famr,b,Q,fyz,control);
  if(DM_iL != 0) DM_step(Yr,famr,b,Q,fyz,control);
  
  arma::mat pD = fyz.get_pD();
  if(flagVERBO) std::cout << "\n Computing standard errors ... " ;
  Rcpp::List se = SEs(Yr,famr,Q,b,pD,flagCORLV);
  if(flagVERBO) std::cout << "\r Computing standard errors ... Done!" ;
  
  arma::cube seb = se["seb"];
  
  if(flagCORLV){
    arma::vec seR = se["seR"];
    return Rcpp::List::create(Rcpp::Named("b") = b,
                              Rcpp::Named("Rz") = Q.get_sigm(),
                              Rcpp::Named("ll") = fyz.get_ll(),
                              Rcpp::Named("seb") = seb,
                              Rcpp::Named("seR") = seR) ; 
  } else{
    return Rcpp::List::create(Rcpp::Named("b") = b,
                              Rcpp::Named("Rz") = Q.get_sigm(),
                              Rcpp::Named("ll") = fyz.get_ll(),
                              Rcpp::Named("seb") = seb) ; 
  }
}

// [[Rcpp::export]]
arma::vec glvmlss_rcpp4sim(Rcpp::NumericMatrix& Yr,
                           Rcpp::CharacterVector& famr,
                           const int q,
                           arma::cube& b,
                           Rcpp::List& control){

  int nqp = control["nQP"];
  bool flagVERBO = control["verbose"], flagCORLV = control["corLV"] ;

  arma::mat Si = arma::eye(q,q);
  if(flagCORLV){
    Si.fill(0.05);
    Si.diag().fill(1.0);
  }
  
  ghQ Q(q,nqp,Si);
  fYZ fyz(Yr,famr,Q,b);
  
  int EM_iL = control["EM.iter.lim"];
  int DM_iL = control["DM.iter.lim"];
  
  if(EM_iL != 0) EM_step(Yr,famr,b,Q,fyz,control);
  if(DM_iL != 0) DM_step(Yr,famr,b,Q,fyz,control);
  
  arma::mat pD = fyz.get_pD();
  if(flagVERBO) std::cout << "\n Computing standard errors ... " ;
  Rcpp::List se = SEs(Yr,famr,Q,b,pD,flagCORLV);
  if(flagVERBO) std::cout << "\r Computing standard errors ... Done!" ;
  
  arma::cube seb = se["seb"];
  
  if(flagCORLV){
    arma::vec seR = se["seR"];
    arma::vec phis = Q.get_phis();
    arma::vec bvec = Cb2Vb(b);
    arma::vec bsev = Cb2Vb(seb);
    arma::vec out = arma::join_cols(bvec,phis,bsev,seR);
    return out;
  } else{
    arma::vec bvec = Cb2Vb(b);
    arma::vec bsev = Cb2Vb(seb);
    arma::vec out = arma::join_cols(bvec,bsev);
    return out;
  }
}

