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

arma::vec d1ll(Rcpp::NumericMatrix& Yr,
               Rcpp::CharacterVector& famr,
               ghQ& Q,
               arma::cube& b,
               arma::mat& pD){

  int p = Yr.ncol(), n = Yr.nrow() ;
  arma::mat Yc(Yr.begin(),n,p,false);
  arma::mat Qp = Q.get_grid() ;
  arma::vec Qw = Q.get_weig() ;
  int tqp = Qp.n_rows ;

  arma::mat nbs = count_b(b,p);
  int npars = arma::accu(nbs);

  arma::vec sc(npars);
  int mem = 0;

  for(int i = 0; i < p; i++){
    std::string fam = Rcpp::as<std::string>(famr(i));
    Dist Y(fam) ;
    arma::cube d1Yi = Y.d1y(i,Yc,Qp,b);
    int K = d1Yi.n_slices;
    for(int kk = 0; kk < K; kk++){
      arma::mat tmpb = b.slice(kk);
      arma::mat tmpM = d1Yi.slice(kk);
      arma::rowvec vec1(tqp) ;
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(tmpM(jj,0))){
          continue;
        } else {
          vec1 += tmpM.row(jj) % pD.row(jj) ;
        }
      }
      arma::uvec idx1 = arma::find(tmpb.row(i) != 0);
      arma::vec vec2 = ((vec1.t() % Qw).t() * Qp.cols(idx1)).t() ;
      arma::uvec idx2 = arma::regspace<arma::uvec>(mem, (mem + nbs(i,kk) - 1));
      sc(idx2) = vec2 ;
      mem += nbs(i,kk);
    }
  }
  return sc;
}

arma::mat d2ll_EM(Rcpp::NumericMatrix& Yr,
                  Rcpp::CharacterVector& famr,
                  ghQ& Q,
                  arma::cube& b,
                  arma::mat& pD,
                  const bool flagEIM){
  
  int p = Yr.ncol(), n = Yr.nrow() ;
  arma::mat Yc(Yr.begin(),n,p,false);
  arma::mat Qp = Q.get_grid() ;
  arma::vec Qw = Q.get_weig() ;
  int tqp = Qp.n_rows ;

  arma::mat nbs = count_b(b,p);
  int npars = arma::accu(nbs);

  arma::mat hess(npars, npars);
  int mem1 = 0, mem2 = 0;

  for(int i = 0; i < p; i++){
    std::string fam = Rcpp::as<std::string>(famr(i));
    Dist Y(fam) ;
    arma::cube d1Yi, dCYi ;
    if(!flagEIM){
      d1Yi = Y.d1y(i,Yc,Qp,b);
    }
    arma::cube d2Yi = Y.d2y(i,Yc,Qp,b,flagEIM,d1Yi);
    int K = d2Yi.n_slices;
    if(K > 1){
      dCYi = Y.dCy(i,Yc,Qp,b,flagEIM);
    }
    
    for(int k1 = 0; k1 < K; k1++){
      arma::mat bk1 = b.slice(k1);
      arma::uvec idk1 = arma::find(bk1.row(i) != 0);
      arma::mat Qpk1 = Qp.cols(idk1) ;
      arma::uvec id1 = arma::regspace<arma::uvec>(mem1, (mem1 + nbs(i,k1) - 1));
      
      arma::mat dYk1 = d2Yi.slice(k1);
      arma::rowvec veck1a(tqp) ;
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(dYk1(jj,0))){
          continue;
        } else {
          veck1a += dYk1.row(jj) % pD.row(jj) ;
        }
      }
      arma::vec veck1b = (veck1a.t() % Qw) ;
      
      for(int k2 = k1; k2 < K; k2++){
        if(k2 == k1){
          arma::cube h1(idk1.n_rows,idk1.n_rows,tqp);
          for(int rr = 0; rr < tqp; rr++){
            h1.slice(rr) = (Qpk1.row(rr).t() * Qpk1.row(rr)) * arma::as_scalar(veck1b(rr)) ;
            hess(id1, id1) += h1.slice(rr);
          }
          mem1 += nbs(i,k2);
          mem2 = mem1;
        } else {
          arma::mat bk2 = b.slice(k2);
          arma::uvec idk2 = arma::find(bk2.row(i) != 0);
          arma::mat Qpk2 = Qp.cols(idk2) ;
          arma::uvec id2 = arma::regspace<arma::uvec>(mem2, (mem2 + nbs(i,k2) - 1));
          
          arma::uword idH = ixd2ll(k1,k2,K);
          arma::mat dYk2 = dCYi.slice(idH);
          arma::rowvec veck2a(tqp) ;
          for(int jj = 0; jj < n; jj++){
            if(std::isnan(dYk2(jj,0))){
              continue;
            } else {
              veck2a += dYk2.row(jj) % pD.row(jj) ;
            }
          }
          arma::vec veck2b = (veck2a.t() % Qw) ;
          
          arma::cube h2(idk1.n_rows,idk2.n_rows,tqp);
          for(int rr = 0; rr < tqp; rr++){
            h2.slice(rr) = (Qpk1.row(rr).t() * Qpk2.row(rr)) * arma::as_scalar(veck2b(rr)) ;
            hess(id1, id2) += h2.slice(rr);
            hess(id2, id1) += h2.slice(rr).t();
          }
          mem2 += nbs(i,k2);
        }
      }
    }
  }
  return hess;
}

arma::vec solve_EM(Rcpp::NumericMatrix& Yr,
                   Rcpp::CharacterVector& famr,
                   ghQ& Q,
                   arma::cube& b,
                   arma::mat& pD,
                   const bool flagEIM){
  
  int p = Yr.ncol(), n = Yr.nrow() ;
  arma::mat Yc(Yr.begin(),n,p,false);
  arma::mat Qp = Q.get_grid() ;
  arma::vec Qw = Q.get_weig() ;
  int tqp = Qp.n_rows ;

  arma::mat nbs = count_b(b,p);
  int npars = arma::accu(nbs);
  
  arma::vec sc(npars);
  arma::mat hess(npars, npars);
  int mem1 = 0, mem2 = 0;

  for(int i = 0; i < p; i++){
    std::string fam = Rcpp::as<std::string>(famr(i));
    Dist Y(fam) ;
    arma::cube d1Yi = Y.d1y(i,Yc,Qp,b);
    arma::cube d2Yi = Y.d2y(i,Yc,Qp,b,flagEIM,d1Yi);
    arma::cube dCYi ;
    int K = d1Yi.n_slices;
    if(K > 1){
      dCYi = Y.dCy(i,Yc,Qp,b,flagEIM);
    }
    
    for(int k1 = 0; k1 < K; k1++){
      arma::mat bk1 = b.slice(k1);
      arma::uvec idk1 = arma::find(bk1.row(i) != 0);
      arma::mat Qpk1 = Qp.cols(idk1) ;
      arma::uvec id1 = arma::regspace<arma::uvec>(mem1, (mem1 + nbs(i,k1) - 1));
      
      arma::mat dYk0 = d1Yi.slice(k1);
      arma::mat dYk1 = d2Yi.slice(k1);
      arma::rowvec veck0a(tqp), veck1a(tqp) ;
      for(int jj = 0; jj < n; jj++){
        if(std::isnan(dYk1(jj,0))){
          continue;
        } else {
          veck0a += dYk0.row(jj) % pD.row(jj) ;
          veck1a += dYk1.row(jj) % pD.row(jj) ;
        }
      }
      arma::vec veck0b = (veck0a.t() % Qw) ;
      arma::vec veck1b = (veck1a.t() % Qw) ;
      
      for(int k2 = k1; k2 < K; k2++){
        if(k2 == k1){
          for(int rr = 0; rr < tqp; rr++){
            // Update for score vector
            sc(id1) += Qpk1.row(rr).t() * arma::as_scalar(veck0b(rr)) ;
            // Update for information matrix (complete log-likelihood / EM)
            hess(id1, id1) += (Qpk1.row(rr).t() * Qpk1.row(rr)) * arma::as_scalar(veck1b(rr)) ;
          }
          mem1 += nbs(i,k2);
          mem2 = mem1;
        } else {
          arma::mat bk2 = b.slice(k2);
          arma::uvec idk2 = arma::find(bk2.row(i) != 0);
          arma::mat Qpk2 = Qp.cols(idk2) ;
          arma::uvec id2 = arma::regspace<arma::uvec>(mem2, (mem2 + nbs(i,k2) - 1));
          
          arma::uword idH = ixd2ll(k1,k2,K);
          arma::mat dYk2 = dCYi.slice(idH);
          arma::rowvec veck2a(tqp) ;
          for(int jj = 0; jj < n; jj++){
            if(std::isnan(dYk2(jj,0))){
              continue;
            } else {
              veck2a += dYk2.row(jj) % pD.row(jj) ;
            }
          }
          arma::vec veck2b = (veck2a.t() % Qw) ;
          
          for(int rr = 0; rr < tqp; rr++){
            hess(id1, id2) += (Qpk1.row(rr).t() * Qpk2.row(rr)) * arma::as_scalar(veck2b(rr)) ;
            hess(id2, id1) += (Qpk1.row(rr).t() * Qpk2.row(rr)).t() * arma::as_scalar(veck2b(rr)) ;
          }
          mem2 += nbs(i,k2);
        }
      }
    }
  }
  return arma::solve(hess,sc, solve_opts::fast + solve_opts::likely_sympd) ;
}

arma::mat SigmaGrad(ghQ& Q, fYZ& fyz){
  
  arma::mat Sig = Q.get_sigm();
  arma::mat Qpi = Q.get_poin();
  arma::vec Qw = Q.get_weig();
  arma::mat pD = fyz.get_pD();
  
  int q = Sig.n_rows;
  int n = pD.n_rows;
  
  arma::mat L = arma::chol(Sig,"lower");
  arma::uvec Li = arma::trimatl_ind(arma::size(L));
  arma::vec vL = L(Li);
  arma::uvec idx = arma::regspace<arma::uvec>(1,Li.n_elem-1);
  Li = Li(idx);
  vL = vL(idx);

  arma::mat EZm = Ezm(q,Qpi,Qw,pD);
  arma::cube VZm = Vzm(q,Qpi,Qw,pD,EZm);
  
  int ncor = vL.n_elem;

  arma::mat iS = arma::inv_sympd(Sig);
  arma::vec grad(ncor) ;
  
  for(int i = 0; i < ncor; i++){
    arma::uword ii = Li(i); 
    arma::mat D(q,q,arma::fill::zeros);
    D(ii) = 1.0 ;
    arma::mat SiD = arma::solve(Sig,D);
    arma::mat G = SiD*L.t()*iS;
    grad(i) = -n*arma::trace(L.t()*SiD) + EVzm(G,EZm,VZm);
  }
  
  arma::mat out(grad);
  return out ;
}

double fZY(ghQ& Q, fYZ& fyz){
  
  arma::mat Sig = Q.get_sigm();
  arma::mat Qpi = Q.get_poin();
  arma::vec Qw = Q.get_weig();
  arma::mat pD = fyz.get_pD();
  
  int n = pD.n_rows ;
  double out = 0.0 ;
  for(int m = 0; m < n; m++){
    out += arma::accu(log_mvnN(Qpi,Sig) % pD.row(m).t() % Qw);
  }
  return out;
}

class opL {
  // Private class members
private:
  ghQ& Q ;
  fYZ& fyz;
  double objf;
  arma::mat SI, L;
  arma::uvec Li, idx;
  
public:
  // Constructor
  opL(ghQ& Qi, fYZ& fyzi) : Q(Qi), fyz(fyzi), objf(0.0) {
    SI = Q.get_sigm();
    L = arma::chol(SI,"lower");
    Li = arma::trimatl_ind(arma::size(L));
    idx = arma::regspace<arma::uvec>(1,Li.n_elem-1);
    Li = Li(idx);
  }
  
  double EvaluateWithGradient(const arma::mat& vL, arma::mat& grad){
    
    arma::vec tvL(vL);
    L(Li) = tvL;
    fixL(L);
    SI = L*L.t() ;
    Q.set_sigm(SI);
    grad = -SigmaGrad(Q,fyz) ;
    objf = fZY(Q,fyz);
    return -1.0*objf;
  }
  
  double get_fzy(){
    return objf;
  }
  
};

void update_sigm(ghQ& Q, fYZ& fyz){
  
  arma::mat Sig = Q.get_sigm();
  arma::mat L = arma::chol(Sig,"lower");
  arma::uvec Li = arma::trimatl_ind(arma::size(L));
  arma::vec vL = L(Li);
  arma::uvec idx = arma::regspace<arma::uvec>(1,Li.n_elem-1);
  Li = Li(idx);
  vL = vL(idx);
  
  opL optL(Q,fyz);
  ens::L_BFGS lbfgsL;
  lbfgsL.MinGradientNorm() = std::sqrt(arma::datum::eps);
  lbfgsL.Optimize(optL, vL);
  Q.update();
}

Rcpp::List SEs(Rcpp::NumericMatrix& Yr,
               Rcpp::CharacterVector& famr,
               ghQ& Q,
               arma::cube& b,
               arma::mat& pD,
               bool FLAGCORLV){
  
  int p = Yr.ncol(), n = Yr.nrow() ;
  arma::mat Yc(Yr.begin(),n,p,false);
  arma::mat Qp = Q.get_grid() ;
  arma::vec Qw = Q.get_weig() ;
  int q = Q.get_q() ;

  arma::mat nbs = count_b(b,p);
  int npars = arma::accu(nbs), ncor = 0 ;

  arma::mat SSt(npars,npars,arma::fill::zeros) ;
  arma::mat sc(npars,n);
  sc.fill(arma::datum::nan);
  arma::mat grad, scout, L; 
  arma::vec vL;
  int mem = 0;

  for(int i = 0; i < p; i++){
    std::string fam = Rcpp::as<std::string>(famr(i));
    Dist Y(fam) ;
    arma::cube d1Yi = Y.d1y(i,Yc,Qp,b);
    int K = d1Yi.n_slices;
    for(int k = 0; k < K; k++){
      arma::mat tmpb = b.slice(k);
      arma::mat tmpM = d1Yi.slice(k);
      arma::uvec idx1 = arma::find(tmpb.row(i) != 0);
      arma::mat Zk = Qp.cols(idx1);
      for(int m = 0; m < n; m++){
        if(std::isnan(tmpM(m,0))){
          continue;
        } else {
          arma::rowvec vec1 = tmpM.row(m) % pD.row(m) % Qw.t();
          arma::vec vec2 = (vec1 * Zk).t();
          sc(arma::span(mem, (mem + nbs(i,k) - 1)),m) = vec2 ;
        }
      }
      mem += nbs(i,k);
    }
  }
  
  if(!FLAGCORLV){
    scout = sc;
  }
  else{
    arma::mat Sig = Q.get_sigm();
    arma::mat Qpi = Q.get_poin();
    
    L = arma::chol(Sig,"lower");
    arma::uvec Li = arma::trimatl_ind(arma::size(L));
    vL = L(Li);
    arma::uvec idx = arma::regspace<arma::uvec>(1,Li.n_elem-1);
    Li = Li(idx);
    vL = vL(idx);
    
    arma::mat EZm = Ezm(q,Qpi,Qw,pD);
    arma::cube VZm = Vzm(q,Qpi,Qw,pD,EZm);
    
    ncor = vL.n_elem;
    
    arma::mat iS = arma::inv_sympd(Sig);
    grad.zeros(ncor,n) ;
    scout.zeros(npars + ncor,n);
    
    for(int i = 0; i < ncor; i++){
      arma::uword ii = Li(i); 
      arma::mat D(q,q,arma::fill::zeros);
      D(ii) = 1.0 ;
      arma::mat SiD = arma::solve(Sig,D);
      arma::mat G = SiD*L.t()*iS;
      for(int m = 0; m < n; m++){
        grad(i,m) = -n*arma::trace(L.t()*SiD) + arma::trace(G*VZm.slice(m)) + arma::as_scalar(EZm.row(m)*G*EZm.row(m).t());
      }
    }
    scout = arma::join_cols(sc,grad);
    SSt.set_size(npars+ncor,npars+ncor);
    SSt.zeros();
  }

  int nPar = scout.n_rows;
  for(int i = 0; i < nPar; i++){
    for(int j = 0; j < nPar; j++){
      for(int m = 0; m < n; m++){
        SSt(i,j) += (std::isnan(scout(i,m)) | std::isnan(scout(j,m))) ? 0.0 : arma::as_scalar(scout(i,m)) * arma::as_scalar(scout(j,m));
      }
    }
  }

  arma::mat iSSt = arma::inv(SSt);
  arma::uvec idxb = arma::regspace<arma::uvec>(0, npars-1);
  arma::mat iSStb = iSSt(idxb,idxb);

  arma::vec seb = arma::sqrt(arma::diagvec(iSStb));

  if(FLAGCORLV){
    arma::cube SEb = b;
    SEb = Vb2Cb(seb,SEb,q);
    
    int nPhi = q*(q-1)/2;
    arma::uvec idxL = arma::regspace<arma::uvec>(npars, npars+ncor-1);
    arma::mat iSStL = iSSt(idxL,idxL);
    
    Rcpp::NumericVector vLtemp = Rcpp::NumericVector(vL.begin(),vL.end());
    Rcpp::NumericMatrix dLdRtemp = Jchol(vLtemp,q);
    arma::mat dLdR(dLdRtemp.begin(),nPhi,ncor,false,false);
    
    arma::mat SEPhi = dLdR*iSStL*dLdR.t();
    arma::vec SEphi = arma::sqrt(arma::diagvec(SEPhi));
    
    return Rcpp::List::create(Rcpp::Named("seb") = SEb,
                              Rcpp::Named("seR") = SEphi) ;
  } else {
    arma::cube SEb = b;
    SEb = Vb2Cb(seb,SEb,q);
    return Rcpp::List::create(Rcpp::Named("seb") = SEb) ;
  }
}

void EM_step(Rcpp::NumericMatrix& Yr,
             Rcpp::CharacterVector& famr,
             arma::cube& b,
             ghQ& Q,
             fYZ& fyz,
             Rcpp::List& control){

  int iL = control["EM.iter.lim"], nqp = control["nQP"];
  bool flagFISHER = control["Fisher"], flagVERBO = control["verbose"];
  bool flagCORLV = control["corLV"] ;
  double tol = control["tolerance"] ;
  
  int q = Q.get_q();
  arma::vec llK(iL + 1, fill::zeros) ;
  llK(0) = fyz.get_ll();
  arma::mat pD = fyz.get_pD();
  arma::vec Vb = Cb2Vb(b);

  for(int i = 0; i < iL; i++){
    Vb -= solve_EM(Yr,famr,Q, b, pD, flagFISHER);
    b = Vb2Cb(Vb,b,q);
    
    if(flagCORLV) update_sigm(Q,fyz);

    fyz.update(Yr,famr,Q,b);
    llK(i+1) = fyz.get_ll();
    pD = fyz.get_pD();
    double eps0 = llK(i+1) - llK(i);

    if(flagVERBO) std::cout << "\r EM iteration: " << i+1 << ", (llk: " << llK(i+1) << ")";
    if(eps0 < 0){
      std::cout << "\n Warning! No increase in log-likelihood (Suggestion: try DM)." ;
      break;
    }
    if(eps0 < tol) break;
    if (i % 2 == 0) Rcpp::checkUserInterrupt();
  }
}

class mLLK {
  // Private class members
  private:
  Rcpp::NumericMatrix& Yr;
  Rcpp::CharacterVector& famr ;
  ghQ& Q ;
  arma::cube& b ;
  fYZ& fyz;
  double llkO;

  public:
  // Constructor
  mLLK(Rcpp::NumericMatrix& Yi,
       Rcpp::CharacterVector& fami,
       ghQ& Qi, fYZ& fyzi,
       arma::cube& bi) : Yr(Yi), famr(fami), Q(Qi), b(bi), fyz(fyzi), llkO(0.0) {
    fyz.update(Yr,famr,Q,b);
  }

  double EvaluateWithGradient(const arma::mat& Vb, arma::mat& grad){
    
      int q = Q.get_q();
      arma::vec tVb(Vb);
      b = Vb2Cb(tVb,b,q);
      fyz.update(Yr,famr,Q,b);
      arma::mat pD = fyz.get_pD();
      arma::vec tgr1 = -d1ll(Yr,famr,Q,b,pD);
      grad = vec2mat(tgr1,tgr1.n_elem,1);
      llkO = fyz.get_ll();
      return -1.0*llkO ;
  }

  double get_ll(){
    return llkO;
  }

};

void DM_step(Rcpp::NumericMatrix& Yr,
             Rcpp::CharacterVector& famr,
             arma::cube& b,
             ghQ& Q,
             fYZ& fyz,
             Rcpp::List& control){

  int iL = control["DM.iter.lim"] ;
  bool flagCORLV = control["corLV"] , flagVERBO = control["verbose"];
  double tol = control["tolerance"] ;

  arma::vec Vb = Cb2Vb(b);

  if(!flagCORLV){
    
    mLLK llk(Yr,famr,Q,fyz,b);
    ens::L_BFGS lbfgs;
    lbfgs.MaxIterations() = iL;
    lbfgs.MinGradientNorm() = std::sqrt(arma::datum::eps);
    
    if(flagVERBO) std::cout << "\n Direct MML estimation (L-BFGS) ... " ;
    lbfgs.Optimize(llk, Vb);
    if(flagVERBO) std::cout << "\r Direct MML estimation (L-BFGS) ... Converged! (llk: " << llk.get_ll() << ")" ;

  } else {
    
    int sigup = 0;
    double eps = 1, llk0 ;

    while(std::abs(eps) > tol & sigup < 30){
      
      llk0 = fyz.get_ll();
      mLLK llk(Yr,famr,Q,fyz,b);
      ens::L_BFGS lbfgs;
      lbfgs.MaxIterations() = iL;
      lbfgs.MinGradientNorm() = std::sqrt(arma::datum::eps);
      
      if(flagVERBO & (sigup == 0)) std::cout << "\n Direct MML estimation (L-BFGS, Phi updates: " << sigup << ") ... " ;
      lbfgs.Optimize(llk, Vb);
      update_sigm(Q,fyz);
      fyz.update(Yr,famr,Q,b);
      eps = fyz.get_ll() - llk0;
      if(flagVERBO) std::cout << "\r Direct MML estimation (L-BFGS, Phi updates: " << sigup << ") ... Converged! (llk: " << llk.get_ll() << ")" ;
      sigup ++;
      
    }
  }
}
