#include <RcppArmadillo.h>
#include "utils_basic.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

arma::uword ixd2ll(int& k1, int& k2, int& K){
  if(K == 2){
    if((k1 == 0) & (k2 == 1)){return 0;}
  } else if(K == 3){
    if((k1 == 0) & (k2 == 1)){return 0;}
    if((k1 == 0) & (k2 == 2)){return 1;}
    if((k1 == 1) & (k2 == 2)){return 2;}
  } else if(K == 4){
    if((k1 == 0) & (k2 == 1)){return 0;}
    if((k1 == 0) & (k2 == 2)){return 1;}
    if((k1 == 1) & (k2 == 2)){return 2;}
    if((k1 == 1) & (k2 == 2)){return 3;}
    if((k1 == 1) & (k2 == 3)){return 4;}
    if((k1 == 2) & (k2 == 3)){return 5;}
  }
  return 0;
}

arma::mat count_b(arma::cube& b, const int p){
  int K = b.n_slices;
  arma::mat out(p,K);
  for(int k = 0; k < K; k++){
    arma::mat btmp = b.slice(k);
    for(int i = 0; i < p; i++){
      arma::rowvec rtmp = btmp.row(i);
      int q = rtmp.size();
      for(int j = 0; j < q; j++){
        out(i,k) += (!std::isnan(btmp(i,j)) & (btmp(i,j) != 0)) ? 1.0 : 0.0;
      }
    }
  }
  return out;
}

arma::cube lb2Cb(Rcpp::List& b, const int p, const int q){
  int K = b.size() ;
  arma::cube bc(p,q+1,K,fill::zeros);
  for(int i = 0; i < K; i++){
    arma::mat Mc = Rcpp::as<arma::mat>(b(i));
    arma::uvec idx = arma::find(Mc != 0);
    bc.slice(i).elem(idx) = Mc.elem(idx);
  }
  return bc;
}

arma::mat Cb2Mb(arma::cube& b){
  int K = b.n_slices, p = b.n_rows, q = b.n_cols ;
  arma::mat bt(b.slice(0).begin(),p,q,false) ;
  for(int k = 1; k  < K; k++){
    arma::mat bt1(b.slice(k).begin(),p,q,false);
    bt = arma::join_rows(bt,bt1);
  }
  return bt;
}

arma::vec Cb2Vb(arma::cube& b){
  return arma::nonzeros(arma::vectorise(Cb2Mb(b), 1).t()) ;
}

arma::cube Mb2Cb(arma::mat mb,
                 const int q,
                 const int K){
  arma::cube b(mb.n_rows,q+1,K);
  int mem = 0;
  for(int k = 0; k < K; k++){
    arma::uvec id0 = arma::regspace<arma::uvec>(mem, (mem + q));
    b.slice(k) = mb.cols(id0);
    mem += (q+1);
  }
  return b;
}

arma::cube Vb2Cb(arma::vec& vb,
                 arma::cube& b,
                 const int q){
  arma::mat bM = Cb2Mb(b).t();
  arma::uvec id0 = arma::find(bM != 0);
  bM(id0) = vb;
  int k = b.n_slices;
  arma::cube bC = Mb2Cb(bM.t(),q,k);
  return bC;
}

arma::mat vec2mat(arma::vec& x, int nrow, int ncol){
  arma::mat y(x);
  y.reshape(nrow, ncol);
  return y;
}

void fixY(arma::vec& y){
  int n = y.n_elem;
  for(int i = 0; i < n; i++){
    if(y(i) == 1.0) y(i) = 1.0 - std::sqrt(arma::datum::eps);
    if(y(i) == 0.0) y(i) = std::sqrt(arma::datum::eps);
  }
}

arma::vec expitvec(arma::vec eta){
  arma::vec expx = arma::exp(eta);
  return expx / (1.0 + expx);
}

arma::vec logitvec(arma::vec mu){
  return arma::log(mu / (1.0 - mu));
}

arma::vec dMd1Evec(arma::vec eta){
  arma::vec expx = arma::exp(eta);
  return expx / arma::pow(1.0 + expx, 2);
}

arma::vec dMd2Evec(arma::vec mu){
  return 1.0 - 2.0*mu ;
}

arma::vec geta(arma::vec& mu, arma::vec& sg){
  arma::vec a(mu.n_elem);
  a = mu % (1.0 - arma::pow(sg,2)) / arma::pow(sg,2);
  return a;
}

arma::vec getb(arma::vec& mu, arma::vec& sg){
  arma::vec b(mu.n_elem);
  b = (1.0 - mu) % (1.0 - arma::pow(sg,2)) / arma::pow(sg,2);
  return b;
}

arma::vec zeta1(arma::vec& y){
  int n = y.n_elem;
  arma::vec out(n);
  for(int i = 0; i < n; i++){
    if(std::isnan(y(i))){
      out(i) = arma::datum::nan ;
    } else {
      out(i) = R::dnorm(y(i),0,1,1) - R::pnorm(y(i),0,1,1,1);
    }
  }
  return arma::exp(out) ;
}

arma::vec zeta2(arma::vec& y){
  arma::vec temp = zeta1(y);
  return -y % temp - arma::pow(temp,2) ;
}

arma::mat d1mu(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = z/sigma(i) - (nu(i)/sigma(i))*zeta1(w);
  }
  return out ;
}

arma::mat d1sg(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = (-zeta1(w)%w)/sigma(i) + (arma::pow(z,2)-1.0)/sigma(i);
  }
  return out ;
}

arma::mat d1nu(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = zeta1(w) % z;
  }
  return out ;
}

arma::mat d2mu(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = -std::pow(sigma(i),-2) + std::pow(nu(i)/sigma(i),2)*zeta2(w);
  }
  return out ;
}

arma::mat d2sg(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = std::pow(sigma(i),-2) * (1 - 3*arma::pow(z,2) + 2*w% zeta1(w) + arma::pow(w,2)%zeta2(w)) ;
  }
  return out ;
}

arma::mat d2nu(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = zeta2(w) % arma::pow(z,2) ;
  }
  return out ;
}

arma::mat dcms(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = std::pow(sigma(i),-2) * (-2*z + nu(i)*zeta1(w) + std::pow(nu(i),2)*zeta2(w) % z) ;
  }
  return out ;
}

arma::mat dcmn(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = -std::pow(sigma(i),-1) * zeta1(w) - w%zeta2(w)*std::pow(sigma(i),-1) ;
  }
  return out ;
}

arma::mat dcsn(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    out.col(i) = -zeta1(w)%z/sigma(i) - zeta2(w)%w%z/sigma(i) ;
  }
  return out ;
}

arma::mat d2muF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    double i1 = arma::as_scalar(arma::mean(arma::pow(zeta1(w),2)));
    out.col(i).fill(-std::pow(sigma(i),-2) - std::pow(nu(i)/sigma(i),2)*i1) ;
  }
  return out ;
}

arma::mat d2sgF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    double i1 = arma::as_scalar(arma::mean(arma::pow(z % zeta1(w),2)));
    out.col(i).fill(std::pow(sigma(i),-2) * (-2 - std::pow(nu(i),2)*i1)) ;
  }
  return out ;
}

arma::mat d2nuF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    double i1 = arma::as_scalar(arma::mean(arma::pow(z % zeta1(w),2)));
    out.col(i).fill(-i1);
  }
  return out ;
}

arma::mat dcmsF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  double log2pi = std::sqrt(2/M_PI);
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    double i1 = arma::as_scalar(arma::mean(z % arma::pow(zeta1(w),2)));
    double i2 = std::pow(sigma(i),-2) * (-log2pi*nu(i)*(1+
                         2*std::pow(nu(i),2))/std::pow(1+std::pow(nu(i),2),1.5) - std::pow(nu(i),2)*i1) ;
    out.col(i).fill(i2);
  }
  return out ;
}

arma::mat dcmnF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  double log2pi = std::sqrt(2/M_PI);
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    double i1 = arma::as_scalar(arma::mean(z % arma::pow(zeta1(w),2)));
    double i2 = 1/sigma(i)*(-log2pi/std::pow(1+std::pow(nu(i),2),1.5) + nu(i)*i1) ;
    out.col(i).fill(i2);
  }
  return out ;
}

arma::mat dcsnF(arma::vec& y, arma::vec& mu, arma::vec& sigma, arma::vec& nu){
  int n = y.n_elem, tqp = mu.n_elem;
  arma::mat out(n,tqp);
  arma::vec z(tqp), w(tqp);
  for(int i = 0; i < tqp; i++){
    z = (y - mu(i))/sigma(i);
    w = nu(i) * z;
    double i1 = arma::as_scalar(arma::mean(arma::pow(z % zeta1(w),2)));
    double i2 = nu(i)*i1/sigma(i) ;
    out.col(i).fill(i2);
  }
  return out ;
}

double dsn(double& x, double& mu, double& sigma, double& nu){
  double z = (x-mu)/sigma ;
  double d1 = -std::log(2*M_PI)/2 - std::log(sigma) - std::pow(z,2)/2 ;
  double d2 = R::pnorm(nu*z,0,1,1,1) + std::log(2) ;
  return d1 + d2 ;
}

double rsn(double& mu, double& sigma, double& nu){
  double delta = nu/std::sqrt(1+std::pow(nu,2));
  arma::vec tn(2);
  for(int i = 0; i < 2; i++){
    tn(i) = R::rnorm(0.0,1.0);
  }
  double z = delta * std::abs(tn(0)) + std::sqrt(1-std::pow(delta,2)) * tn(1) ;
  return mu + sigma*z  ;
}

arma::vec nc2nd(arma::vec& nuc){
  double log2pi = std::sqrt(2/M_PI);
  double OoT = 1.0/3.0 ;
  arma::vec R = arma::pow((2*arma::abs(nuc))/(4-M_PI),OoT) % arma::sign(nuc);
  return R / arma::sqrt(std::pow(log2pi,2) - (1-std::pow(log2pi,2)) * arma::pow(R,2));
}

arma::vec sc2sd(arma::vec& sgc, arma::vec& nuc){
  double log2pi = std::sqrt(2/M_PI);
  arma::vec nud = nc2nd(nuc) ;
  arma::vec R = nud / arma::sqrt(1+arma::pow(nud,2)) ;
  return sgc / arma::sqrt(1-std::pow(log2pi,2) * arma::pow(R,2)) ;
}

arma::vec mc2md(arma::vec& muc, arma::vec& sgc, arma::vec& nuc){
  double log2pi = std::pow(2/M_PI,0.5);
  arma::vec nud = nc2nd(nuc) ;
  arma::vec R = nud / arma::sqrt(1+arma::pow(nud,2)) ;
  arma::vec sgd = sc2sd(sgc,nuc);
  return muc - log2pi * sgd % R;
}

arma::vec p2nuc(arma::vec& p){
  double g1max = std::sqrt(2)*(4-M_PI)/std::pow(M_PI-2,1.5);
  return(g1max*(2*p - 1));
}

arma::vec dsddsc(arma::vec& nuc){
  double log2pi = std::pow(2/M_PI,0.5);
  arma::vec nud = nc2nd(nuc);
  arma::vec R = nud / arma::sqrt(1+arma::pow(nud,2)) ;
  return 1 / arma::sqrt(1-std::pow(log2pi,2) * arma::pow(R,2)) ;
}

arma::vec dnddnc(arma::vec& nuc){
  double a = 2/(3*(4-M_PI));
  double log2pi = std::pow(2/M_PI,0.5);
  double OoT = 1.0/3.0 ;
  arma::vec R = arma::pow((2*arma::abs(nuc))/(4-M_PI),OoT) % arma::sign(nuc);
  arma::vec t = arma::sqrt(std::pow(log2pi,2) - (1-std::pow(log2pi,2)) * arma::pow(R,2));
  return a * (1/(t % arma::pow(R,2)) + (1-std::pow(log2pi,2))/arma::pow(t,3));
}

arma::vec dmddsc(arma::vec& nuc){
  double log2pi = std::pow(2/M_PI,0.5);
  arma::vec nud = nc2nd(nuc);
  arma::vec R = nud / arma::sqrt(1+arma::pow(nud,2)) ;
  return - (log2pi * R) / arma::sqrt(1 - std::pow(log2pi,2) * arma::pow(R,2)) ;
}

arma::vec dmddnc(arma::vec& sgc, arma::vec& nuc){
  double log2pi = std::pow(2/M_PI,0.5);
  arma::vec nud = nc2nd(nuc);
  arma::vec R = nud / arma::sqrt(1+arma::pow(nud,2)) ;
  return - (sgc % (log2pi * R)) / (3*arma::sqrt(1 - std::pow(log2pi,2) * arma::pow(R,2)) % nuc) ;
}

arma::vec dsddnc(arma::vec& sgc, arma::vec& nuc){
  double log2pi = std::pow(2/M_PI,0.5);
  arma::vec nud = nc2nd(nuc);
  arma::vec R = nud / arma::sqrt(1+arma::pow(nud,2)) ;
  arma::vec mz = log2pi * R ;
  arma::vec sz = arma::sqrt(1 - std::pow(log2pi,2) * arma::pow(R,2)) ;
  arma::vec d1 = -sgc / arma::pow(sz,2);
  arma::vec d2 = -mz/sz % (log2pi / arma::pow(1+arma::pow(nud,2) ,1.5)) ;
  arma::vec d3 = dnddnc(nuc);
  return d1 % d2 % d3;
}

arma::mat Ezm(int q,
              arma::mat& Qpi,
              arma::vec& Qw,
              arma::mat& pD){

  int n = pD.n_rows;
  int tqp = pD.n_cols;

  arma::mat zm(n,q);
  for(int m = 0; m < n; m++){
    for(int i = 0; i < tqp; i++){
      zm.row(m) += Qpi.row(i) * arma::as_scalar(pD(m,i)) * arma::as_scalar(Qw(i)) ;
    }
  }
  return zm ;
}

arma::cube Vzm(int q,
               arma::mat& Qpi,
               arma::vec& Qw,
               arma::mat& pD,
               arma::mat& EZM){

  int n = pD.n_rows;
  int tqp = pD.n_cols;

  arma::cube Vm(q,q,n);
  for(int m = 0; m < n; m++){
    for(int i = 0; i < tqp; i++){
      arma::rowvec tmp = Qpi.row(i) - EZM.row(m) ;
      Vm.slice(m) += (tmp.t()*tmp * arma::as_scalar(pD(m,i)) * arma::as_scalar(Qw(i)));
    }
  }
  return Vm;
}

double EVzm(arma::mat& G,
            arma::mat& EZM,
            arma::cube& VZM){
  int n = EZM.n_rows;
  double out = 0.0 ;
  for(int m = 0; m < n; m++){
    out += arma::trace(G*VZM.slice(m)) + arma::as_scalar(EZM.row(m)*G*EZM.row(m).t());
  }
  return out ;
}

arma::vec log_mvnN(arma::mat& y,
                   arma::mat& S){

  int tqp = y.n_rows ;
  int q = y.n_cols ;
  double lS;
  double sign;
  bool ok = arma::log_det(lS, sign, S);

  arma::vec out(tqp);

  if(ok){
    for(int m = 0; m < tqp ; m++){
      out(m) = -0.5*(q*std::log(2*M_PI) + lS + arma::as_scalar(y.row(m)*(arma::inv_sympd(S)*y.row(m).t()))) ;
    }
  }
  return out;
}

void fixL(arma::mat& L){
  int q = L.n_rows;
  for(int i = 0; i < q; i++){
    double cross = arma::dot(L.row(i),L.row(i));
    L.row(i) = L.row(i)/std::sqrt(cross);
  }
}

Rcpp::NumericVector chol2cor(Rcpp::NumericVector& vLi, int& q){
  int ncor = vLi.length();
  arma::vec vL(vLi.begin(),ncor,false,false);
  arma::mat L = eye(q,q);
  arma::uvec Li = arma::trimatl_ind(arma::size(L));
  arma::uvec idx = arma::regspace<arma::uvec>(1,Li.n_elem-1);
  Li = Li(idx);
  L(Li) = vL;
  arma::mat R = L*L.t();
  arma::uvec Ri = arma::trimatl_ind(arma::size(R),-1);
  arma::vec out = R(Ri);
  return Rcpp::wrap(out);
}

Rcpp::NumericMatrix Jchol(Rcpp::NumericVector& vLi,
                          int q) {

  Rcpp::Environment numDeriv("package:numDeriv");
  Rcpp::Function jacobian = numDeriv["jacobian"];
  Rcpp::NumericMatrix out = jacobian(Rcpp::_["func"] = Rcpp::InternalFunction(chol2cor),
                                     Rcpp::_["x"] = vLi,
                                     Rcpp::_["q"] = q);

  return out;
}
