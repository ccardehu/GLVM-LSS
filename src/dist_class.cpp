#include <RcppArmadillo.h>
#include "utils_basic.h"
#include "dist_class.h"
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

using namespace Rcpp;
using namespace arma;
using namespace boost::math;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

arma::vec Dist::sim(const int n,
                    int i,
                    arma::cube& b,
                    arma::mat& Z){

  arma::colvec res(n);

  if(fam == "normal"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));

    for(int idx = 0; idx < n; idx++){
      res(idx) = R::rnorm(mu(idx),sg(idx));
    }
    return res;
  } else if(fam == "binomial"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));

    for(int idx = 0; idx < n; idx++){
      res(idx) = R::rbinom(1,mu(idx));
    }
    return res;
  } else if(fam == "beta"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));
    arma::vec sg = expitvec(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec a = geta(mu,sg);
    arma::vec b = getb(mu,sg);

    for(int idx = 0; idx < n; idx++){
      res(idx) = R::rbeta(a(idx),b(idx));
    }
    return res;
  } else if(fam == "sn"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::vec bnu = b.slice(2).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::uvec bnu0 = arma::find(bnu != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec nu = expitvec(Z.cols(bnu0)*bnu.elem(bnu0));

    arma::vec nuc = p2nuc(nu);
    arma::vec nud = nc2nd(nuc);
    arma::vec sgd = sc2sd(sg,nuc);
    arma::vec mud = mc2md(mu,sg,nuc);

    for(int idx = 0; idx < n; idx++){
      res(idx) = rsn(mud(idx),sgd(idx),nud(idx));
    }
    return res;
  } else return res;
}

arma::mat Dist::dy(const int i,
                   arma::mat& Yr,
                   arma::mat& QP,
                   arma::cube& b){

  arma::vec  Yc = Yr.col(i) ;
  arma::mat  Z(QP) ;

  int n = Yc.size() ;
  int tqp = Z.n_rows ;
  arma::mat pD(n,tqp, arma::fill::zeros) ;

  if(fam == "normal"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));

    for(int jj = 0; jj < n; jj++){
      if(std::isnan(Yc(jj))){
        pD.row(jj).fill(arma::datum::nan) ;
      } else {
        for(int ii = 0; ii < tqp; ii++){
          double mud = mu(ii);
          double sgd = sg(ii);
          pD(jj,ii) = R::dnorm(Yc(jj), mud, sgd, 1);
        }
      }
    }
    return pD;
  } else if(fam == "binomial"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));

    for(int jj = 0; jj < n; jj++){
      if(std::isnan(Yc(jj))){
        pD.row(jj).fill(arma::datum::nan) ;
      } else {
        for(int ii = 0; ii < tqp; ii++){
          double mud = mu(ii);
          pD(jj,ii) = R::dbinom(Yc(jj), 1, mud, 1);
        }
      }
    }
    return pD;
  } else if(fam == "beta"){
    fixY(Yc);
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));
    arma::vec sg = expitvec(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec al = geta(mu,sg);
    arma::vec be = getb(mu,sg);

    for(int jj = 0; jj < n; jj++){
      if(std::isnan(Yc(jj))){
        pD.row(jj).fill(arma::datum::nan) ;
      } else {
        for(int ii = 0; ii < tqp; ii++){
          double ad = al(ii);
          double bd = be(ii);
          pD(jj,ii) = R::dbeta(Yc(jj), ad, bd, 1);
        }
      }
    }
    return pD;
  } else if(fam == "sn"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::vec bnu = b.slice(2).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::uvec bnu0 = arma::find(bnu != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec nu = expitvec(Z.cols(bnu0)*bnu.elem(bnu0));

    arma::vec nuc = p2nuc(nu);
    arma::vec nud = nc2nd(nuc);
    arma::vec sgd = sc2sd(sg,nuc);
    arma::vec mud = mc2md(mu,sg,nuc);

    for(int jj = 0; jj < n; jj++){
      if(std::isnan(Yc(jj))){
        pD.row(jj).fill(arma::datum::nan) ;
      } else {
        for(int ii = 0; ii < tqp; ii++){
          pD(jj,ii) = dsn(Yc(jj), mud(ii), sgd(ii), nud(ii));
        }
      }
    }
    return pD;
  } else return pD;
}

arma::cube Dist::d1y(const int i,
                     arma::mat& Yr,
                     arma::mat& QP,
                     arma::cube& b){

  arma::vec  Yc = Yr.col(i) ;
  int n = Yc.size() ;
  arma::mat  Z(QP) ;
  int tqp = Z.n_rows ;

  arma::cube D ;

  if(fam == "normal"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));

    D.set_size(n,tqp,2);

    for(int jj = 0; jj < n; jj++){
      if(std::isnan(Yc(jj))){
        D.row(jj).fill(arma::datum::nan) ;
      } else {
        for(int ii = 0; ii < tqp; ii++){
          double mud = mu(ii);
          double sgd = sg(ii);
          // First derivatives of mu: d1fdmu * dmudeta
          D(jj,ii,0) = (Yc(jj) - mud)/std::pow(sgd,2); // * 1
          // First derivatives of sg: d1fdsg * dsgdeta
          D(jj,ii,1) = (std::pow(Yc(jj) - mud,2) - std::pow(sgd,2))/std::pow(sgd,3) * sgd;
        }
      }
    }
    return D;
  } else if(fam == "binomial"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));
    arma::vec dX = dMd1Evec(logitvec(mu));

    D.set_size(n,tqp,1);

    for(int jj = 0; jj < n; jj++){
      if(std::isnan(Yc(jj))){
        D.row(jj).fill(arma::datum::nan) ;
      } else {
        for(int ii = 0; ii < tqp; ii++){
          double mud = mu(ii);
          double dXd = dX(ii);
          // First derivatives of mu: d1fdmu * dmudeta
          D(jj,ii,0) = (Yc(jj) - mud)/(mud * (1 - mud)) * dXd ;
        }
      }
    }
    return D;
  } else if(fam == "beta"){
    fixY(Yc);
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));
    arma::vec sg = expitvec(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec al = geta(mu,sg);
    arma::vec be = getb(mu,sg);
    arma::vec dXm = dMd1Evec(logitvec(mu));
    arma::vec dXs = dMd1Evec(logitvec(sg));

    D.set_size(n,tqp,2);

    for(int jj = 0; jj < n; jj++){
      if(std::isnan(Yc(jj))){
        D.row(jj).fill(arma::datum::nan) ;
      } else {
        for(int ii = 0; ii < tqp; ii++){
          double mud = mu(ii);
          double sgd = sg(ii);
          double ad = al(ii);
          double bd = be(ii);
          double dXmd = dXm(ii);
          double dXsd = dXs(ii);
          // First derivatives of mu: d1fdmu * dmudeta
          D(jj,ii,0) = ((1.0 - std::pow(sgd,2))/std::pow(sgd,2)) * (-1.0*boost::math::digamma(ad) +
            boost::math::digamma(bd) + std::log(Yc(jj)) - std::log(1.0-Yc(jj))) * dXmd ;
          // First derivatives of sg: d1fdsg * dsgdeta
          D(jj,ii,1) = -1.0*(2/(std::pow(sgd,3))) * (mud * (-1.0*boost::math::digamma(ad) + boost::math::digamma(ad + bd) +
            std::log(Yc(jj))) + (1-mud)*(-1.0*boost::math::digamma(bd) + boost::math::digamma(ad + bd) +
            std::log(1.0 - Yc(jj))) ) * dXsd ;
        }
      }
    }
    return D;
  } else if(fam == "sn"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::vec bnu = b.slice(2).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::uvec bnu0 = arma::find(bnu != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec nu = expitvec(Z.cols(bnu0)*bnu.elem(bnu0));
    arma::vec dXn = dMd1Evec(logitvec(nu));

    arma::vec nuc = p2nuc(nu);
    arma::vec nud = nc2nd(nuc);
    arma::vec sgd = sc2sd(sg,nuc);
    arma::vec mud = mc2md(mu,sg,nuc);

    double dnucdp = (2*std::sqrt(2)*(4-M_PI)/std::pow(M_PI-2,1.5));

    D.set_size(n,tqp,3);

    arma::mat D1M = d1mu(Yc,mud,sgd,nud);
    arma::mat D1S = d1sg(Yc,mud,sgd,nud);
    arma::mat D1N = d1nu(Yc,mud,sgd,nud);

    arma::mat D1Ms = D1M * arma::diagmat(dmddsc(nuc)) ;
    arma::mat D1Mn = D1M * arma::diagmat(dmddnc(sg,nuc)) ;
    arma::mat D1Ss = D1S * arma::diagmat(dsddsc(nuc)) ;
    arma::mat D1Sn = D1S * arma::diagmat(dsddnc(sg,nuc)) ;
    D1N = D1N * arma::diagmat(dnddnc(nuc));

    D.slice(0) = D1M ;
    D.slice(1) = (D1Ms + D1Ss) * arma::diagmat(sg) ;
    D.slice(2) = (D1Mn + D1Sn + D1N) * arma::diagmat(dnucdp * dXn);
    return D;
  } else return D;
}

arma::cube Dist::d2y(const int i,
                     arma::mat& Yr,
                     arma::mat& QP,
                     arma::cube& b,
                     const bool flagEIM,
                     arma::cube& d1Yi){

  arma::vec  Yc = Yr.col(i) ;
  int n = Yc.size() ;
  arma::mat  Z(QP) ;
  int tqp = Z.n_rows ;

  arma::cube D ;

  if(fam == "normal"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));

    D.set_size(n,tqp,2);

    if(flagEIM){
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(Yc(jj))){
          D.row(jj).fill(arma::datum::nan) ;
        } else {
          for(int ii = 0; ii < tqp; ii++){
            // double mud = mu(ii);
            double sgd = sg(ii);
            // Second derivatives of mu: d2fdmu * dmudeta^2
            D(jj,ii,0) = -1/std::pow(sgd,2); // * 1
            // Second derivatives of sg: d2fdsg * dsgdeta^2
            D(jj,ii,1) = -2/std::pow(sgd,2) * std::pow(sgd,2); // exp(log(sgd)) ^ 2 = sgd^2;
          }
        }
      }
      return D;
    } else {
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(Yc(jj))){
          D.row(jj).fill(arma::datum::nan) ;
        } else {
          for(int ii = 0; ii < tqp; ii++){
            double mud = mu(ii);
            double sgd = sg(ii);
            // Second derivatives of mu: (d2fdmu * dmudeta^2) + d1fdmu
            D(jj,ii,0) = -1/std::pow(sgd,2); // * 1
            // Second derivatives of sg: (d2fdsg * dsgdeta^2) + d1fdsg
            D(jj,ii,1) = (-3*std::pow(Yc(jj) - mud,2) / std::pow(sgd,4) + 1/std::pow(sgd,2)) * std::pow(sgd,2) + d1Yi(jj,ii,1) ; // exp(log(sgd)) ^ 2 = sgd^2;
          }
        }
      }
      return D;
    }
  } else if(fam == "binomial"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));
    arma::vec d1X = dMd1Evec(logitvec(mu));
    arma::vec d2X = dMd2Evec(mu);

    D.set_size(n,tqp,1);

    if(flagEIM){
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(Yc(jj))){
          D.row(jj).fill(arma::datum::nan) ;
        } else {
          for(int ii = 0; ii < tqp; ii++){
            double mud = mu(ii);
            double d1Xd = d1X(ii);
            // Second derivatives of mu: d2fdmu * dmudeta^2
            D(jj,ii,0) = -1/(mud * (1-mud)) * std::pow(d1Xd,2) ;
          }
        }
      }
      return D;
    } else {
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(Yc(jj))){
          D.row(jj).fill(arma::datum::nan) ;
        } else {
          for(int ii = 0; ii < tqp; ii++){
            double mud = mu(ii);
            double d1Xd = d1X(ii);
            double d2Xd = d2X(ii);
            // Second derivatives of mu: (d2fdmu * dmudeta^2) + d1fdmu
            D(jj,ii,0) = (-(std::pow(mud,2) - 2*Yc(jj)*mud + Yc(jj)) / std::pow((mud-1)*mud,2) * std::pow(d1Xd,2) + d2Xd * d1Yi(jj,ii,0)) ;
          }
        }
      }
      return D;
    }
  } else if(fam == "beta"){
    fixY(Yc);
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));
    arma::vec sg = expitvec(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec al = geta(mu,sg);
    arma::vec be = getb(mu,sg);
    arma::vec d1Xm = dMd1Evec(logitvec(mu));
    arma::vec d1Xs = dMd1Evec(logitvec(sg));
    arma::vec d2Xm = dMd2Evec(mu);
    arma::vec d2Xs = dMd2Evec(sg);

    D.set_size(n,tqp,2);

    if(flagEIM){
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(Yc(jj))){
          D.row(jj).fill(arma::datum::nan) ;
        } else {
          for(int ii = 0; ii < tqp; ii++){
            double mud = mu(ii);
            double sgd = sg(ii);
            double ad = al(ii);
            double bd = be(ii);
            double d1Xmd = d1Xm(ii);
            double d1Xsd = d1Xs(ii);
            // Second derivatives of mu: d2fdmu * dmudeta^2
            D(jj,ii,0) = -1.0*(std::pow(1.0 - std::pow(sgd,2),2)/std::pow(sgd,4)) * (boost::math::trigamma(ad) +
                               boost::math::trigamma(bd)) * std::pow(d1Xmd,2) ;
            // Second derivatives of sg: d2fdsg * dsgdeta^2
            D(jj,ii,1) = -1.0*(4.0/std::pow(sgd,6)) * (std::pow(mud,2) * boost::math::trigamma(ad) +
                               std::pow(1.0 - mud,2) * boost::math::trigamma(bd) - boost::math::trigamma(ad + bd)) * std::pow(d1Xsd,2) ;
          }
        }
      }
      return D;
    } else {
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(Yc(jj))){
          D.row(jj).fill(arma::datum::nan) ;
        } else {
          for(int ii = 0; ii < tqp; ii++){
            double mud = mu(ii);
            double sgd = sg(ii);
            double ad = al(ii);
            double bd = be(ii);
            double d1Xmd = d1Xm(ii);
            double d1Xsd = d1Xs(ii);
            double d2Xmd = d2Xm(ii);
            double d2Xsd = d2Xs(ii);
            // Second derivatives of mu: (d2fdmu * dmudeta^2) + d1fdmu
            D(jj,ii,0) = -1.0*(std::pow(1.0 - std::pow(sgd,2),2)/std::pow(sgd,4)) * (boost::math::trigamma(ad) +
                               boost::math::trigamma(bd)) * std::pow(d1Xmd,2) + d2Xmd * d1Yi(jj,ii,0) ;
            // Second derivatives of sg: (d2fdsg * dsgdeta^2) + d1fdsg
            D(jj,ii,1) = -1.0*(4.0/std::pow(sgd,6)) * (std::pow(mud,2) * boost::math::trigamma(ad) +
                               std::pow(1.0 - mud,2) * boost::math::trigamma(bd) - boost::math::trigamma(ad + bd)) * std::pow(d1Xsd,2) + d2Xsd * d1Yi(jj,ii,1) ;
          }
        }
      }
      return D;
    }
  } else if(fam == "sn"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::vec bnu = b.slice(2).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::uvec bnu0 = arma::find(bnu != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec nu = expitvec(Z.cols(bnu0)*bnu.elem(bnu0));
    arma::vec d1Xn = dMd1Evec(logitvec(nu));
    arma::vec d2Xn = dMd2Evec(nu);

    arma::vec nuc = p2nuc(nu);
    arma::vec nud = nc2nd(nuc);
    arma::vec sgd = sc2sd(sg,nuc);
    arma::vec mud = mc2md(mu,sg,nuc);

    double dnucdp = (2*std::sqrt(2)*(4-M_PI)/std::pow(M_PI-2,1.5));

    D.set_size(n,tqp,3);

    if(flagEIM){
      arma::mat D2M = d2muF(Yc,mud,sgd,nud);
      arma::mat D2S = d2sgF(Yc,mud,sgd,nud);
      arma::mat D2N = d2nuF(Yc,mud,sgd,nud);
      arma::mat DcMS = dcmsF(Yc,mud,sgd,nud);
      arma::mat DcMN = dcmnF(Yc,mud,sgd,nud);
      arma::mat DcSN = dcsnF(Yc,mud,sgd,nud);

      arma::mat D2Ms = D2M * arma::diagmat(dmddsc(nuc)) ;
      arma::mat D2Mn = D2M * arma::diagmat(dmddnc(sg,nuc)) ;
      arma::mat DcMSs1 = DcMS * arma::diagmat(dsddsc(nuc));
      arma::mat DcMSs2 = DcMS * arma::diagmat(dmddsc(nuc));
      arma::mat DcMSn1 = DcMS * arma::diagmat(dsddnc(sg,nuc));
      arma::mat DcMSn2 = DcMS * arma::diagmat(dmddnc(sg,nuc));
      arma::mat DcMNn1 = DcMN * arma::diagmat(dnddnc(nuc));
      arma::mat DcMNn2 = DcMN * arma::diagmat(dmddnc(sg,nuc));

      arma::mat D2Ss = D2S * arma::diagmat(dsddsc(nuc)) ;
      arma::mat D2Sn = D2S * arma::diagmat(dsddnc(sg,nuc)) ;
      arma::mat DcSNn1 = DcSN * arma::diagmat(dnddnc(nuc));
      arma::mat DcSNn2 = DcSN * arma::diagmat(dsddnc(sg,nuc));

      D2N = D2N * arma::diagmat(dnddnc(nuc)) ;

      D.slice(0) = D2M ;
      D.slice(1) = ((D2Ms + DcMSs1) * arma::diagmat(dmddsc(nuc)) +
                    (DcMSs2 + D2Ss) * arma::diagmat(dsddsc(nuc))) * arma::diagmat(arma::pow(sg,2)) ;
      D.slice(2) = ((D2Mn + DcMSn1 + DcMNn1) * arma::diagmat(dmddnc(sg,nuc)) +
                    (DcMSn2 + D2Sn + DcSNn1) * arma::diagmat(dsddnc(sg,nuc)) +
                    (DcMNn2 + DcSNn2 + D2N) * arma::diagmat(dnddnc(nuc))) * arma::diagmat(arma::pow(dnucdp * d1Xn,2) ) ;
      return D;
    } else {
      arma::mat D2M = d2mu(Yc,mud,sgd,nud);
      arma::mat D2S = d2sg(Yc,mud,sgd,nud);
      arma::mat D2N = d2nu(Yc,mud,sgd,nud);
      arma::mat DcMS = dcms(Yc,mud,sgd,nud);
      arma::mat DcMN = dcmn(Yc,mud,sgd,nud);
      arma::mat DcSN = dcsn(Yc,mud,sgd,nud);

      arma::mat D2Ms = D2M * arma::diagmat(dmddsc(nuc)) ;
      arma::mat D2Mn = D2M * arma::diagmat(dmddnc(sg,nuc)) ;
      arma::mat DcMSs1 = DcMS * arma::diagmat(dsddsc(nuc));
      arma::mat DcMSs2 = DcMS * arma::diagmat(dmddsc(nuc));
      arma::mat DcMSn1 = DcMS * arma::diagmat(dsddnc(sg,nuc));
      arma::mat DcMSn2 = DcMS * arma::diagmat(dmddnc(sg,nuc));
      arma::mat DcMNn1 = DcMN * arma::diagmat(dnddnc(nuc));
      arma::mat DcMNn2 = DcMN * arma::diagmat(dmddnc(sg,nuc));

      arma::mat D2Ss = D2S * arma::diagmat(dsddsc(nuc)) ;
      arma::mat D2Sn = D2S * arma::diagmat(dsddnc(sg,nuc)) ;
      arma::mat DcSNn1 = DcSN * arma::diagmat(dnddnc(nuc));
      arma::mat DcSNn2 = DcSN * arma::diagmat(dsddnc(sg,nuc));

      D2N = D2N * arma::diagmat(dnddnc(nuc)) ;

      D.slice(0) = D2M ;
      D.slice(1) = ((D2Ms + DcMSs1) * arma::diagmat(dmddsc(nuc)) +
                    (DcMSs2 + D2Ss) * arma::diagmat(dsddsc(nuc))) * arma::diagmat(arma::pow(sg,2)) + d1Yi.slice(1) ;
      D.slice(2) = ((D2Mn + DcMSn1 + DcMNn1) * arma::diagmat(dmddnc(sg,nuc)) +
                    (DcMSn2 + D2Sn + DcSNn1) * arma::diagmat(dsddnc(sg,nuc)) +
                    (DcMNn2 + DcSNn2 + D2N) * arma::diagmat(dnddnc(nuc))) * arma::diagmat(arma::pow(dnucdp * d1Xn,2))
                    + (d1Yi.slice(2) * arma::diagmat(d2Xn)) ;
      return D;
    }
  } else return D;
}

arma::cube Dist::dCy(const int i,
                     arma::mat& Yr,
                     arma::mat& QP,
                     arma::cube& b,
                     const bool flagEIM){

  arma::vec  Yc = Yr.col(i) ;
  int n = Yc.size() ;
  arma::mat  Z(QP) ;
  int tqp = Z.n_rows ;

  arma::cube D ;

  if(fam == "normal"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));

    D.set_size(n,tqp,1);

    if(flagEIM){
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(Yc(jj))){
          D.row(jj).fill(arma::datum::nan) ;
        } else {
          D.row(jj).fill(0.0);
        }
      }
      return D;
    } else {
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(Yc(jj))){
          D.row(jj).fill(arma::datum::nan) ;
        } else {
          for(int ii = 0; ii < tqp; ii++){
            double mud = mu(ii);
            double sgd = sg(ii);
            // Cross derivatives of mu: (d2f/dmudsg * dmudeta * dsgdeta)
            D(jj,ii,0) = -2*(Yc(jj) - mud)/std::pow(sgd,3) * sgd; // * 1 ;
          }
        }
      }
      return D;
    }
  } else if(fam == "binomial"){
    return D;
  } else if(fam == "beta"){
    fixY(Yc);
    // Yc.clamp(arma::datum::eps,1.0-arma::datum::eps);
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::vec mu = expitvec(Z.cols(bmu0)*bmu.elem(bmu0));
    arma::vec sg = expitvec(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec al = geta(mu,sg);
    arma::vec be = getb(mu,sg);
    arma::vec d1Xm = dMd1Evec(logitvec(mu));
    arma::vec d1Xs = dMd1Evec(logitvec(sg));

    D.set_size(n,tqp,1);

    for(int jj = 0; jj < n; jj++){
      if(std::isnan(Yc(jj))){
        D.row(jj).fill(arma::datum::nan) ;
      } else {
        for(int ii = 0; ii < tqp; ii++){
          double mud = mu(ii);
          double sgd = sg(ii);
          double ad = al(ii);
          double bd = be(ii);
          double d1Xmd = d1Xm(ii);
          double d1Xsd = d1Xs(ii);
          // Cross derivatives of mu: (d2f/dmudsg * dmudeta * dsgdeta)
          D(jj,ii,0) = (2.0 * (1.0 - std::pow(sgd,2))/std::pow(sgd,5)) * (mud * boost::math::trigamma(ad) - (1.0 - mud)*boost::math::trigamma(bd)) * d1Xmd * d1Xsd ;
        }
      }
    }
    return D;
  } else if(fam == "sn"){
    arma::vec bmu = b.slice(0).row(i).t() ;
    arma::vec bsg = b.slice(1).row(i).t() ;
    arma::vec bnu = b.slice(2).row(i).t() ;
    arma::uvec bmu0 = arma::find(bmu != 0) ; // add here if != nan
    arma::uvec bsg0 = arma::find(bsg != 0) ;
    arma::uvec bnu0 = arma::find(bnu != 0) ;
    arma::vec mu = Z.cols(bmu0)*bmu.elem(bmu0);
    arma::vec sg = exp(Z.cols(bsg0)*bsg.elem(bsg0));
    arma::vec nu = expitvec(Z.cols(bnu0)*bnu.elem(bnu0));
    arma::vec dXn = dMd1Evec(logitvec(nu));

    arma::vec nuc = p2nuc(nu);
    arma::vec nud = nc2nd(nuc);
    arma::vec sgd = sc2sd(sg,nuc);
    arma::vec mud = mc2md(mu,sg,nuc);

    double dnucdp = (2*std::sqrt(2)*(4-M_PI)/std::pow(M_PI-2,1.5));

    D.set_size(n,tqp,3);

    if(flagEIM){
      arma::mat D2M = d2muF(Yc,mud,sgd,nud);
      arma::mat D2S = d2sgF(Yc,mud,sgd,nud);
      arma::mat DcMS = dcmsF(Yc,mud,sgd,nud);
      arma::mat DcMN = dcmnF(Yc,mud,sgd,nud);
      arma::mat DcSN = dcsnF(Yc,mud,sgd,nud);

      arma::mat D2Ms = D2M * arma::diagmat(dmddsc(nuc)) ;
      arma::mat D2Mn = D2M * arma::diagmat(dmddnc(sg,nuc)) ;
      arma::mat DcMSs1 = DcMS * arma::diagmat(dsddsc(nuc));
      arma::mat DcMSn1 = DcMS * arma::diagmat(dsddnc(sg,nuc));
      arma::mat DcMSn2 = DcMS * arma::diagmat(dmddnc(sg,nuc));
      arma::mat DcMNn1 = DcMN * arma::diagmat(dnddnc(nuc));

      arma::mat D2Sn = D2S * arma::diagmat(dsddnc(sg,nuc)) ;
      arma::mat DcSNn1 = DcSN * arma::diagmat(dnddnc(nuc));

      D.slice(0) = (D2Ms + DcMSs1) * arma::diagmat(sg) ;
      D.slice(1) = (D2Mn + DcMSn1 + DcMNn1) * arma::diagmat(dXn) ;
      D.slice(2) = ((D2Mn + DcMSn1 + DcMNn1) * arma::diagmat(dmddsc(nuc)) +
                    (DcMSn2 + D2Sn + DcSNn1) * arma::diagmat(dsddsc(nuc))) * arma::diagmat(sg) * arma::diagmat(dnucdp * dXn);
      return D;
    } else {
      arma::mat D2M = d2mu(Yc,mud,sgd,nud);
      arma::mat D2S = d2sg(Yc,mud,sgd,nud);
      arma::mat DcMS = dcms(Yc,mud,sgd,nud);
      arma::mat DcMN = dcmn(Yc,mud,sgd,nud);
      arma::mat DcSN = dcsn(Yc,mud,sgd,nud);

      arma::mat D2Ms = D2M * arma::diagmat(dmddsc(nuc)) ;
      arma::mat D2Mn = D2M * arma::diagmat(dmddnc(sg,nuc)) ;
      arma::mat DcMSs1 = DcMS * arma::diagmat(dsddsc(nuc));
      arma::mat DcMSn1 = DcMS * arma::diagmat(dsddnc(sg,nuc));
      arma::mat DcMSn2 = DcMS * arma::diagmat(dmddnc(sg,nuc));
      arma::mat DcMNn1 = DcMN * arma::diagmat(dnddnc(nuc));

      arma::mat D2Sn = D2S * arma::diagmat(dsddnc(sg,nuc)) ;
      arma::mat DcSNn1 = DcSN * arma::diagmat(dnddnc(nuc));

      D.slice(0) = (D2Ms + DcMSs1) * arma::diagmat(sg) ;
      D.slice(1) = (D2Mn + DcMSn1 + DcMNn1) * arma::diagmat(dXn) ;
      D.slice(2) = ((D2Mn + DcMSn1 + DcMNn1) * arma::diagmat(dmddsc(nuc)) +
                    (DcMSn2 + D2Sn + DcSNn1) * arma::diagmat(dsddsc(nuc))) * (arma::diagmat(sg) * arma::diagmat(dnucdp * dXn));
      return D;
    }
  } else return D;
}

