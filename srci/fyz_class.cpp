#include <RcppArmadillo.h>
#include "utils_basic.h"
#include "ghq_class.h"
#include "dists.h"
#include "fyz_class.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void fYZ::update(Rcpp::NumericMatrix& Yr,
                 Rcpp::CharacterVector& famr,
                 ghQ& Q,
                 arma::cube& b){
  
  int p = Yr.ncol(), n = Yr.nrow() ;
  arma::mat Yc(Yr.begin(),n,p,false);
  arma::mat Qp = Q.get_grid() ;
  arma::vec Qw = Q.get_weig() ;
  int tqp = Qp.n_rows ;
  
  pD.fill(0.0);

  for(int i = 0; i < p; i++){
    std::string fam = Rcpp::as<std::string>(famr[i]);
    Dist Y(fam) ;
    arma::mat dYi = Y.dy(i,Yc,Qp,b);
    for(int ii = 0; ii < tqp; ii++){
      for(int jj = 0; jj < n; jj++){
        pD(jj,ii) += std::isnan(dYi(jj,ii)) ? 0.0 : dYi(jj,ii) ;
      }
    }
  }
   
  arma::vec llv = arma::exp(pD)*Qw;
  ll = arma::accu(arma::log(llv));
  pD = arma::exp(pD) ;
  pD.each_col() /= llv;
}

arma::mat fYZ::get_pD(){
    return pD;
}

double fYZ::get_ll(){
    return ll;
}