#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

NumericVector prb (NumericVector x) {
 const double eps = pow(2.220446e-16,0.5) ; 
 int p = x.size() ; 
 NumericVector res(p); 
 for(int i = 0; i < p; i ++){
  res[i] = R::plogis(x[i], 0.0, 1.0, true, false) ;
  if(res[i] == 1) { res[i] = 1 - eps ; } ;
  if(res[i] == 0) { res[i] = eps ; } ;
 }
 return(res);
}

Rcpp::NumericVector Clogit (NumericVector x) {
 int n = x.size();
 NumericVector result(n);
 for(int i = 0; i < n; ++i) {
  result[i] = log( x[i] / (1.0 - x[i]) );
 }
 return result;
}

NumericVector dZIPo (NumericVector Y,
                     NumericVector mu,
                     NumericVector sg,
                     bool rlog = true) {
  
  NumericVector u(Y.size(), 0.0);
  NumericVector lf(Y.size()); 
  
  for(int i = 0; i < Y.size(); i++){
    if(Y[i] == 0) u[i] = 1; 
  }
  for(int j = 0; j < Y.size(); j++){
    lf[j] = u[j]*std::log(sg[j] + (1-sg[j])*exp(-mu[j])) + (1-u[j])*(log(1-sg[j]) - mu[j] + Y[j]*log(mu[j]) - std::lgamma(Y[j] + 1)) ;
  }
  if(rlog == true) return lf ;
  else  {
   return exp(lf);
  } 

}

arma::mat Zreg(DataFrame df,
               Formula formula) {
 Rcpp::Environment stats_env("package:stats");
 Rcpp::Function model_matrix = stats_env["model.matrix"];
 arma::mat df_new = as<arma::mat>(model_matrix(_["object"] = formula, _["data"] = df));
 return(df_new);
}

NumericMatrix make_mat(List input_list){
   
   unsigned int n = input_list.length();
   if(n == 0) { 
      Rcpp::stop("Must supply a list with more than 1 element.");
   }
   Rcpp::NumericVector testvals = input_list[0];
   unsigned int elems = testvals.length();
   Rcpp::NumericMatrix result_mat = Rcpp::no_init(elems,n);
   // fill by column
   for(unsigned int i = 0; i < n; i++) {
    Rcpp::NumericVector row_val = input_list[i];
    if(elems != row_val.length()) {
     Rcpp::stop("Length of row does not match matrix requirements"); 
    }
   result_mat(_,i) = row_val;
   }
 return result_mat;
}

NumericMatrix row_erase(NumericMatrix &x,
                        IntegerVector &rowID) {
   rowID = rowID.sort();
   NumericMatrix x2(Dimension(x.nrow()- rowID.size(), x.ncol()));
   int iter = 0; 
   int del = 1; // to count deleted elements
   for (int i = 0; i < x.nrow(); i++) {
      if (i != rowID[del - 1]) {
         x2.row(iter) = x.row(i);
         iter++;
      } else {
         del++;
      }
   }
   return x2;
}

arma::vec rowS(const arma::mat &X){
 int nRows = X.n_rows;
 arma::vec out(nRows);
 for(int i = 0; i < nRows; i++){
  out(i) = sum(X.row(i));
 }
return(out);
}

arma::rowvec colS(const arma::mat &X) {
 int nCols = X.n_cols;
 arma::rowvec out(nCols);
 for(int i = 0; i < nCols; i++){
  out(i) = sum(X.col(i));
 }
 return(out);
}

LogicalVector tP (CharacterVector x,
                  CharacterVector par){
 Rcpp::LogicalVector r(x.size()) ;
 for(int i = 0; i < x.size(); i++){
  r[i] = (x[i] == par[0]) ;
 }
 return any(r); 
}

// [[Rcpp::export]]
Rcpp::List mvGHq (int dm,
                  int qp,
                  Nullable<NumericVector> mu = R_NilValue,
                  Nullable<NumericMatrix> Sg = R_NilValue,
                  Nullable<double> prune = R_NilValue) {
  
 Environment GHq("package:fastGHQuad") ;
 Function ghq = GHq["gaussHermiteData"] ;
 Function expGrid("expand.grid") ;
 Function p0("paste0") ;
 const double pi = M_PI ; 
 
 List uGHq = ghq(qp) ;
 NumericVector x = uGHq["x"] ;
 NumericVector w = uGHq["w"] ;
 w = w * 1/pow(2*pi,0.5) * exp(pow(x,2)/2) ;
 CharacterVector nam = p0("Z", seq_len(dm)) ;
 int qlim = pow(qp, dm) ;
   
 Rcpp::List tmp1 = expGrid(rep(List::create(seq_len(qp)-1), dm)) ;
 NumericMatrix poi(make_mat(tmp1)) ;
 NumericMatrix ghp(qlim,dm) ;
 NumericVector ghw(qlim) ;
 NumericMatrix tmp2(qlim,dm) ; 

 for(int it1 = 0 ; it1 < ghp.nrow(); it1 ++) {
  for(int it2 = 0; it2 < ghp.ncol(); it2 ++) {
   ghp(it1,it2) = x[poi(it1,it2)] ;
   tmp2(it1,it2) = w[poi(it1,it2)] ; 
  }
  ghw[it1] = std::accumulate(tmp2(it1,_).begin(), tmp2(it1,_).end(), 1.0, std::multiplies<double>());
 }

 Environment stats("package:stats") ;
 Function quantile = stats["quantile"] ;
 Function Rwhich("which");
 Function Rrev("rev");
 Function Rmatrix("matrix"); 
 NumericVector lim(ghw.size()) ;
 IntegerVector tt(ghw.size()), tta(ghw.size()) ; 
 NumericMatrix out(ghw.size(), ghp.ncol()) ;
 arma::mat eiVe(dm,dm) ;
 arma::vec eiVa(dm) ;
 arma::mat rot(dm,dm) ;
   
 if(prune.isNotNull()){
  lim = rep(as<NumericVector>(quantile(ghw,prune)), ghw.size()); 
  tt = as<NumericVector>(Rwhich(ghw > lim))-1 ;
  tta = as<NumericVector>(Rwhich(ghw <= lim))-1 ;
  ghw = ghw[tt] ;
  for(int it3 = 0; it3 < tt.size(); it3 ++) {
   out(tt[it3], _) = ghp(tt[it3], _) ;
  }
  ghp = row_erase(out, tta); 
 }

 if(Sg.isNotNull()){
  arma::eig_sym(eiVa,eiVe,as<arma::mat>(Sg)) ;
  eiVa = as<arma::vec>(Rrev(wrap(eiVa))); 
  rot = eiVe.cols(as<arma::uvec>(Rrev(seq_len(eiVe.n_cols)-1))) * diagmat(pow(eiVa,0.5)) ;
 if(mu.isNotNull()){
  ghp = wrap(trans(rot * trans(as<arma::mat>(ghp)) + as<arma::mat>(Rmatrix(mu,dm,ghp.nrow(),false)))) ;
 }
 else {
  ghp = wrap(trans(rot * trans(as<arma::mat>(ghp)))) ;
 }
}
 colnames(ghp) = nam ;
 return(List::create(_["x"] = as<DataFrame>(ghp), _["w"] = ghw)) ;
}

// [[Rcpp::export]]
arma::mat cdF (arma::mat Y,
               Rcpp::List b,
               Rcpp::List form,
               Rcpp::CharacterVector fam,
               Rcpp::List mGQobj,
               int evqp) {
                     
 int n = Y.n_rows ;
 int p = Y.n_cols ;
 arma::mat R(n,p, arma::fill::zeros) ;
 arma::mat Zmu ;
 arma::mat Zsg ;
 arma::mat Zta ;
 arma::mat Znu ;
 arma::mat bmu ;
 arma::mat bsg ;
 arma::mat bta ;
 arma::mat bnu ;
 NumericVector mu ;
 NumericVector sg ;
 NumericVector tau ;
 NumericVector nu ;
 
 Zmu = Zreg(as<DataFrame>(mGQobj["x"]), as<Formula>(form["mu"])).row(evqp) ;
 bmu = as<arma::mat>(make_mat(b["mu"])) ;
 bmu.reshape(p, Zmu.n_cols) ;

// Algorithm:

 for(int i = 0; i < p; i++){

 if(fam[i] == "normal"){
  Zsg = Zreg(as<DataFrame>(mGQobj["x"]), as<Formula>(form["sigma"])).row(evqp) ;
  bsg = as<arma::mat>(make_mat(b["sigma"])) ;
  bsg.reshape(p, Zsg.n_cols) ;
  mu = wrap(Zmu * bmu.row(i).t()) ;
  sg = wrap(Zsg * bsg.row(i).t()) ;
  for(int j = 0; j < n; j ++) {
   R(j,i) = R::dnorm(Y(j,i), as<double>(mu), exp(as<double>(sg)), true) ;
  }
 }
  
 if(fam[i] == "poisson"){
  mu = wrap(Zmu * bmu.row(i).t()) ;
  for(int j = 0; j < n; j ++) {
   R(j,i) = R::dpois(Y(j,i), exp(as<double>(mu)), true) ;
  }
 }
 
 if(fam[i] == "gamma"){
  Zsg = Zreg(as<DataFrame>(mGQobj["x"]), as<Formula>(form["sigma"])).row(evqp) ;
  bsg = as<arma::mat>(make_mat(b["sigma"])) ;
  bsg.reshape(p, Zsg.n_cols) ;
  mu = wrap(Zmu * bmu.row(i).t()) ;
  sg = wrap(Zsg * bsg.row(i).t()) ;
  for(int j = 0; j < n; j ++) {
   R(j,i) = R::dgamma(Y(j,i), exp(as<double>(mu)), exp(as<double>(sg)), true) ;
  }
 }
 
 if(fam[i] == "binom"){
  mu = wrap(Zmu * bmu.row(i).t()) ;
  for(int j = 0; j < n; j ++) {
   R(j,i) = R::dbinom(Y(j,i), 1, as<double>(prb(mu)), true) ;
  }
 }
 
 if(fam[i] == "ZIpoisson"){
  Zsg = Zreg(as<DataFrame>(mGQobj["x"]), as<Formula>(form["sigma"])).row(evqp) ;
  bsg = as<arma::mat>(make_mat(b["sigma"])) ;
  bsg.reshape(p, Zsg.n_cols) ; 
  mu = wrap(Zmu * bmu.row(i).t()) ;
  sg = wrap(Zsg * bsg.row(i).t()) ;
  for(int j = 0; j < n; j ++) {
   R(j,i) = as<double>(dZIPo(Y(j,i), exp(mu), prb(sg), true)) ;
  }
 }
 
 }
 return(R) ;
}

// [[Rcpp::export]]
arma::vec mdF (arma::mat Y,
               Rcpp::List b,
               Rcpp::List form,
               Rcpp::CharacterVector fam,
               Rcpp::List mGQobj) {

 int n = Y.n_rows ;
 int p = Y.n_cols ;
 NumericVector ws = as<NumericVector>(mGQobj["w"]) ;
 int e = ws.size();
 
 arma::cube cfYz(n,p,e) ;
 arma::mat tM(n,1) ;
 arma::cube Res0(n,1,e, arma::fill::zeros);
 arma::vec Res1(n, arma::fill::zeros); 
 
 for(int it1 = 0; it1 < e; it1 ++){
  tM.fill(ws[it1]);
  cfYz.slice(it1) = cdF(Y,b,form,fam,mGQobj,it1) ; 
  Res0.slice(it1) = exp(rowS(cfYz.slice(it1))) % tM ;
  Res1 += Res0.slice(it1) ;
 }

 return(Res1); 
}

// [[Rcpp::export]]
Rcpp::List dvl (arma::mat Y,
                Rcpp::List b,
                Rcpp::List form,
                Rcpp::CharacterVector fam,
                Rcpp::List mGQobj) {

int n = Y.n_rows ;
int p = Y.n_cols ;
NumericVector ws = as<NumericVector>(mGQobj["w"]) ;
int e = ws.size();

// Data containers -----

arma::cube d1mu(n,p,e, arma::fill::zeros) ;
arma::cube d1sg(n,p,e, arma::fill::zeros) ;
arma::cube d1ta(n,p,e, arma::fill::zeros) ;
arma::cube d1nu(n,p,e, arma::fill::zeros) ;
arma::cube d2mu(n,p,e, arma::fill::zeros) ;
arma::cube d2sg(n,p,e, arma::fill::zeros) ;
arma::cube d2ta(n,p,e, arma::fill::zeros) ;
arma::cube d2nu(n,p,e, arma::fill::zeros) ;

arma::mat Zmu ;
arma::mat Zsg ;
arma::mat Zta ;
arma::mat Znu ;
arma::mat bmu ;
arma::mat bsg ;
arma::mat bta ;
arma::mat bnu ;
NumericVector mu ;
NumericVector sg ;
NumericVector tau ;
NumericVector nu ;

// ---------------------

bmu = as<arma::mat>(make_mat(b["mu"])) ;

for(int evqp = 0; evqp < e; evqp ++) {
 Zmu = Zreg(as<DataFrame>(mGQobj["x"]), as<Formula>(form["mu"])).row(evqp) ;
 bmu.reshape(p, Zmu.n_cols) ;
 
 for(int i = 0; i < p; i ++){
  
  if(fam[i] == "normal"){
   Zsg = Zreg(as<DataFrame>(mGQobj["x"]), as<Formula>(form["sigma"])).row(evqp) ;
   bsg = as<arma::mat>(make_mat(b["sigma"])) ;
   bsg.reshape(p, Zsg.n_cols) ;
   mu = wrap(Zmu * bmu.row(i).t()) ;
   sg = exp(Zsg * bsg.row(i).t()) ;
   for(int j = 0; j < n; j ++) {
    d1mu(j,i,evqp) = (Y(j,i) - as<double>(mu)) * pow(as<double>(sg),-2) ;
    d2mu(j,i,evqp) = -pow(as<double>(sg),-2) ;
    d1sg(j,i,evqp) = (pow((Y(j,i) - as<double>(mu)),2) - pow(as<double>(sg),2)) * pow(as<double>(sg),-3) * as<double>(sg) ;
    d2sg(j,i,evqp) = -2*pow(as<double>(sg),-2) * pow(as<double>(sg),2) ; 
   }
  }
  
  if(fam[i] == "poisson"){
   mu = exp(Zmu * bmu.row(i).t()) ;
   for(int j = 0; j < n; j ++) {
    d1mu(j,i,evqp) = ((Y(j,i) - as<double>(mu))/as<double>(mu)) * as<double>(mu) ;
    d2mu(j,i,evqp) = -pow(as<double>(mu),-1) * pow(as<double>(mu),2) ;
   }
  }
  
  if(fam[i] == "gamma"){
   Zsg = Zreg(as<DataFrame>(mGQobj["x"]), as<Formula>(form["sigma"])).row(evqp) ;
   bsg = as<arma::mat>(make_mat(b["sigma"])) ;
   bsg.reshape(p, Zsg.n_cols) ;
   mu = exp(Zmu * bmu.row(i).t()) ;
   sg = exp(Zsg * bsg.row(i).t()) ;
   for(int j = 0; j < n; j ++) {
    d1mu(j,i,evqp) = (log(Y(j,i)) - R::digamma(as<double>(mu)) - log(as<double>(sg))) * as<double>(mu) ;
    d2mu(j,i,evqp) = - R::trigamma(as<double>(mu)) * pow(as<double>(mu),2) ;
    d1sg(j,i,evqp) = (Y(j,i)*pow(as<double>(sg),-2) -as<double>(mu)*pow(as<double>(sg),-1)) * as<double>(sg) ;
    d2sg(j,i,evqp) = -as<double>(mu)*pow(as<double>(sg),-1)*pow(as<double>(sg),2) ;
   }
  }
  
  if(fam[i] == "binom"){
   mu = prb(wrap(Zmu * bmu.row(i).t())) ;
   for(int j = 0; j < n; j ++) {
    d1mu(j,i,evqp) = (Y(j,i) - as<double>(mu))*pow(as<double>(mu)*(1 - as<double>(mu)),-1)*as<double>(mu)*(1 - as<double>(mu)) ;
    d2mu(j,i,evqp) = -pow(as<double>(mu)*(1 - as<double>(mu)),-1) * pow(as<double>(mu)*(1 - as<double>(mu)),2) ;
   }
  }
  
  if(fam[i] == "ZIpoisson"){
   Zsg = Zreg(as<DataFrame>(mGQobj["x"]), as<Formula>(form["sigma"])).row(evqp) ;
   bsg = as<arma::mat>(make_mat(b["sigma"])) ;
   bsg.reshape(p, Zsg.n_cols) ;
   arma::vec u(n) ;  
   u.fill(0.0) ;
   mu = prb(wrap(Zmu * bmu.row(i).t())) ;
   sg = exp(Zsg * bsg.row(i).t()) ;
   for(int j = 0; j < n; j ++) {
    if(Y(j,i) == 0){ u[j] = pow(1 + exp(-as<double>(Clogit(sg)) -as<double>(mu)),-1) ; }
     d1mu(j,i,evqp) = (1-u[j])*(Y(j,i)/as<double>(mu)-1) * as<double>(mu) ;
     d2mu(j,i,evqp) = (1-u[j])* -pow(as<double>(mu),-1) * pow(as<double>(mu),2) ; // expected value is -1/mu; original is (-Y/mu^2)
     d1sg(j,i,evqp) = (u[j]-as<double>(sg)) * pow(as<double>(sg)*(1-as<double>(sg)),-1) * as<double>(sg)*(1-as<double>(sg)) ;
     d2sg(j,i,evqp) = -(u[j]*pow(as<double>(sg),-2) + (1-u[j])*pow(1-as<double>(sg),-2)) * pow(as<double>(sg)*(1-as<double>(sg)),2) ; 
   }
  }
  
 }
 
}

List dmu = Rcpp::List::create(_["d1"] = d1mu, _["d2"] = d2mu); 
List dsg = Rcpp::List::create(_["d1"] = d1sg, _["d2"] = d2sg); 
List dta = Rcpp::List::create(_["d1"] = d1ta, _["d2"] = d2ta);
List dnu = Rcpp::List::create(_["d1"] = d1nu, _["d2"] = d2nu);

// if(form.size() == 1) { return(Rcpp::List::create(_["mu"] = dmu)) ; };
// if(form.size() == 2) { return(Rcpp::List::create(_["mu"] = dmu, _["sigma"] = dsg)) ; };
// if(form.size() == 3) { return(Rcpp::List::create(_["mu"] = dmu, _["sigma"] = dsg, _["tau"] = dta)) ; };
// if(form.size() == 4) {
  return(Rcpp::List::create(_["mu"] = dmu, _["sigma"] = dsg, _["tau"] = dta, _["nu"] = dnu)) ;
  // };

}

// [[Rcpp::export]]
Rcpp::List NRb (arma::mat &Y,
                Rcpp::List form,
                Rcpp::CharacterVector fam,
                int dm,
                int qp,
                int maxit = 100,
                double tol = 1e-7){
                
 int n = Y.n_rows ;
 int p = Y.n_cols ;
 List qOb = mvGHq(dm,qp) ;
 NumericVector ws = as<NumericVector>(qOb["w"]) ;
 int e = ws.size() ;
 
 // Data containers
 
 arma::mat  Zsg, Zta, Znu, bsg, bta , bnu, hsgt, htat, hnut ;
 arma::vec  bsgt, btat, bnut, vLL0, vLL1 ;
 arma::cube bSsg, bSta , bSnu, bHsg, bHta, bHnu ;
 arma::cube cDs(n,p,e,arma::fill::zeros) ;
 arma::cube _mu1(1,p,e,arma::fill::zeros) ;
 arma::cube _mu2(1,p,e,arma::fill::zeros) ;
 arma::cube _sg1(1,p,e,arma::fill::zeros) ;
 arma::cube _sg2(1,p,e,arma::fill::zeros) ;
 arma::cube _ta1(1,p,e,arma::fill::zeros) ;
 arma::cube _ta2(1,p,e,arma::fill::zeros) ;
 arma::cube _nu1(1,p,e,arma::fill::zeros) ;
 arma::cube _nu2(1,p,e,arma::fill::zeros) ;
 arma::mat tM(n,p) ;
 List b1, D ; 
 double llk0, llk1 ;
 double eps = 2*tol ;
 int iter = 0; 
 
 arma::mat  Zmu = Zreg(as<DataFrame>(qOb["x"]), as<Formula>(form["mu"])) ;
 arma::mat  bmu(p, Zmu.n_cols, arma::fill::ones) ;
 arma::cube bSmu(Zmu.n_cols, 1, e, arma::fill::zeros) ;
 arma::vec  bmut(Zmu.n_cols, arma::fill::zeros) ;
 arma::cube bHmu(Zmu.n_cols, Zmu.n_cols, e, arma::fill::zeros) ;
 arma::mat  hmut(Zmu.n_cols, Zmu.n_cols, arma::fill::zeros) ;
 List b0 = List::create(_["mu"] = bmu) ;

 if(form.size() == 2) {
  Zsg = Zreg(as<DataFrame>(qOb["x"]), as<Formula>(form["sigma"])) ;
  bsg.ones(p, Zsg.n_cols) ;
  b0.push_back(bsg, "sigma") ;
  bSsg.zeros(Zsg.n_cols, 1, e) ;
  bHsg.zeros(Zsg.n_cols, Zsg.n_cols, e) ;
  bsgt.zeros(Zsg.n_cols);
  hsgt.zeros(Zsg.n_cols, Zsg.n_cols);
 };
 if(form.size() == 3) {
  Zsg = Zreg(as<DataFrame>(qOb["x"]), as<Formula>(form["sigma"])) ;
  bsg.zeros(p, Zsg.n_cols) ;
  b0.push_back(bsg, "sigma") ;
  bSsg.zeros(Zsg.n_cols, 1, e) ;
  bHsg.zeros(Zsg.n_cols, Zsg.n_cols, e) ;
  Zta = Zreg(as<DataFrame>(qOb["x"]), as<Formula>(form["tau"])) ;
  bta.zeros(p, Zta.n_cols) ;
  b0.push_back(bta, "tau") ;
  bSta.zeros(Zta.n_cols, 1, e) ;
  bHta.zeros(Zta.n_cols, Zta.n_cols, e) ;
  bsgt.zeros(Zsg.n_cols);
  hsgt.zeros(Zsg.n_cols);
  btat.zeros(Zta.n_cols);
  htat.zeros(Zta.n_cols);
 };
 if(form.size() == 4) {
  Zsg = Zreg(as<DataFrame>(qOb["x"]), as<Formula>(form["sigma"])) ;
  bsg.zeros(p, Zsg.n_cols) ;
  b0.push_back(bsg, "sigma") ;
  bSsg.zeros(Zsg.n_cols, 1, e) ;
  bHsg.zeros(Zsg.n_cols, Zsg.n_cols, e) ;
  Zta = Zreg(as<DataFrame>(qOb["x"]), as<Formula>(form["tau"])) ;
  bta.zeros(p, Zta.n_cols) ;
  b0.push_back(bta, "tau") ;
  bSta.zeros(Zta.n_cols, 1, e) ;
  bHta.zeros(Zta.n_cols, Zta.n_cols, e) ;
  Znu = Zreg(as<DataFrame>(qOb["x"]), as<Formula>(form["nu"])) ;
  bnu.zeros(p, Znu.n_cols) ;
  b0.push_back(bnu, "nu") ;
  bSnu.zeros(Znu.n_cols, 1, e) ;
  bHnu.zeros(Znu.n_cols, Znu.n_cols, e) ;
  bsgt.zeros(Zsg.n_cols);
  hsgt.zeros(Zsg.n_cols);
  btat.zeros(Zta.n_cols);
  htat.zeros(Zta.n_cols);
  bnut.zeros(Znu.n_cols);
  hnut.zeros(Znu.n_cols);
 };
 for(int i0 = 0; i0 < e; i0++) {
 bSmu.slice(i0) = Zmu.row(i0).t() ;
 bHmu.slice(i0) = Zmu.row(i0).t() * Zmu.row(i0) ; 
 if(form.size() == 2) {
  bSsg.slice(i0) = Zsg.row(i0).t() ;
  bHsg.slice(i0) = Zsg.row(i0).t() * Zsg.row(i0) ;
 } ;
 if(form.size() == 3) {
  bSsg.slice(i0) = Zsg.row(i0).t() ;
  bHsg.slice(i0) = Zsg.row(i0).t() * Zsg.row(i0) ;
  bSta.slice(i0) = Zta.row(i0).t() ;
  bHta.slice(i0) = Zta.row(i0).t() * Zta.row(i0) ;
 } ;
 if(form.size() == 4) {
  bSsg.slice(i0) = Zsg.row(i0).t() ;
  bHsg.slice(i0) = Zsg.row(i0).t() * Zsg.row(i0) ;
  bSta.slice(i0) = Zta.row(i0).t() ;
  bHta.slice(i0) = Zta.row(i0).t() * Zta.row(i0) ;
  bSnu.slice(i0) = Znu.row(i0).t() ;
  bHnu.slice(i0) = Znu.row(i0).t() * Znu.row(i0) ;
 } ;
 }

 while(iter < maxit && eps > tol) {

 // Preparation
 
 if(iter % 10 == 0) Rcpp::checkUserInterrupt();
 vLL0 = mdF(Y,b0,form,fam,qOb) ;
 llk0 = as<double>(wrap(sum(log(vLL0)))) ; 
 D = dvl(Y,b0,form,fam,qOb) ;
 List Dmu = D["mu"] ;
 
 for(int i0 = 0; i0 < e; i0++) {
 tM.fill(ws[i0]); 
 for(int i00 = 0; i00 < p; i00++) {
  cDs.slice(i0).col(i00) = exp(rowS(cdF(Y,b0,form,fam,qOb,i0))) % (1/vLL0) % tM.col(i00) ;
 }
 _mu1.slice(i0) = colS(as<arma::cube>(Dmu["d1"]).slice(i0) % cDs.slice(i0)) ;
 _mu2.slice(i0) = colS(as<arma::cube>(Dmu["d2"]).slice(i0) % cDs.slice(i0)) ;
 if(form.size() == 2) {
  List Dsg = D["sigma"] ;
  _sg1.slice(i0) = colS(as<arma::cube>(Dsg["d1"]).slice(i0) % cDs.slice(i0)) ;
  _sg2.slice(i0) = colS(as<arma::cube>(Dsg["d2"]).slice(i0) % cDs.slice(i0)) ;
 } ;
 if(form.size() == 3) {
  List Dsg = D["sigma"] ;
  _sg1.slice(i0) = colS(as<arma::cube>(Dsg["d1"]).slice(i0) % cDs.slice(i0)) ;
  _sg2.slice(i0) = colS(as<arma::cube>(Dsg["d2"]).slice(i0) % cDs.slice(i0)) ;
  List Dta = D["tau"] ;
  _ta1.slice(i0) = colS(as<arma::cube>(Dta["d1"]).slice(i0) % cDs.slice(i0)) ;
  _ta2.slice(i0) = colS(as<arma::cube>(Dta["d2"]).slice(i0) % cDs.slice(i0)) ;
 } ;
 if(form.size() == 4) {
  List Dsg = D["sigma"] ;
  _sg1.slice(i0) = colS(as<arma::cube>(Dsg["d1"]).slice(i0) % cDs.slice(i0)) ;
  _sg2.slice(i0) = colS(as<arma::cube>(Dsg["d2"]).slice(i0) % cDs.slice(i0)) ;
  List Dta = D["tau"] ;
  _ta1.slice(i0) = colS(as<arma::cube>(Dta["d1"]).slice(i0) % cDs.slice(i0)) ;
  _ta2.slice(i0) = colS(as<arma::cube>(Dta["d2"]).slice(i0) % cDs.slice(i0)) ;
  List Dnu = D["nu"] ;
  _nu1.slice(i0) = colS(as<arma::cube>(Dnu["d1"]).slice(i0) % cDs.slice(i0)) ;
  _nu2.slice(i0) = colS(as<arma::cube>(Dnu["d2"]).slice(i0) % cDs.slice(i0)) ;
 } ;
 } ;

 for(int i1 = 0; i1 < p; i1++){

 // E-Step

 for(int i2 = 0; i2 < e; i2++){
  bmut += bSmu.slice(i2) * as<double>(wrap(_mu1.slice(i2).col(i1))) ; 
  hmut += bHmu.slice(i2) * as<double>(wrap(_mu2.slice(i2).col(i1))) ; 
 if(form.size() == 2) {
  bsgt += bSsg.slice(i2) * as<double>(wrap(_sg1.slice(i2).col(i1))) ;
  hsgt += bHsg.slice(i2) * as<double>(wrap(_sg2.slice(i2).col(i1))) ;
 } ;
 if(form.size() == 3) {
  bsgt += bSsg.slice(i2) * as<double>(wrap(_sg1.slice(i2).col(i1))) ;
  hsgt += bHsg.slice(i2) * as<double>(wrap(_sg2.slice(i2).col(i1))) ;
  btat += bSta.slice(i2) * as<double>(wrap(_ta1.slice(i2).col(i1))) ;
  htat += bHta.slice(i2) * as<double>(wrap(_ta2.slice(i2).col(i1))) ;
 } ;
 if(form.size() == 4) {
  bsgt += bSsg.slice(i2) * as<double>(wrap(_sg1.slice(i2).col(i1))) ;
  hsgt += bHsg.slice(i2) * as<double>(wrap(_sg2.slice(i2).col(i1))) ;
  btat += bSta.slice(i2) * as<double>(wrap(_ta1.slice(i2).col(i1))) ;
  htat += bHta.slice(i2) * as<double>(wrap(_ta2.slice(i2).col(i1))) ;
  bnut += bSnu.slice(i2) * as<double>(wrap(_nu1.slice(i2).col(i1))) ;
  hnut += bHnu.slice(i2) * as<double>(wrap(_nu2.slice(i2).col(i1))) ;
  } ;
 }

 // M-step

 bmu.row(i1) -= (inv_sympd(hmut)*bmut).t() ;
 if(form.size() == 2) {
  bsg.row(i1) -= (inv_sympd(hsgt)*bsgt).t() ;
 } ;
 if(form.size() == 3) {
  bsg.row(i1) -= (inv_sympd(hsgt)*bsgt).t() ;
  bta.row(i1) -= (inv_sympd(htat)*btat).t() ;
 } ;
 if(form.size() == 4) {
  bsg.row(i1) -= (inv_sympd(hsgt)*bsgt).t() ;
  bta.row(i1) -= (inv_sympd(htat)*btat).t() ;
  bnu.row(i1) -= (inv_sympd(hnut)*bnut).t() ;
 } ;
 
 bmut.fill(0.0); hmut.fill(0.0);
 bsgt.fill(0.0); hsgt.fill(0.0);
 btat.fill(0.0); htat.fill(0.0);
 bnut.fill(0.0); hnut.fill(0.0);

 }
 
 b1 = List::create(_["mu"] = bmu) ;
 if(form.size() == 2) {
  b1.push_back(bsg, "sigma") ;
 } ;
 if(form.size() == 3) {
   b1.push_back(bsg, "sigma") ;
   b1.push_back(bta, "tau") ;
 } ;
 if(form.size() == 4) {
   b1.push_back(bsg, "sigma") ;
   b1.push_back(bta, "tau") ;
   b1.push_back(bnu, "nu") ;
 } ;
 
 vLL1 = mdF(Y,b1,form,fam,qOb) ;
 llk1 = as<double>(wrap(sum(log(vLL1)))) ;
 eps = abs(llk1 - llk0) ;
 iter++ ;
 b0 = b1 ;
 Rprintf("\r EM - Interation: %i", iter) ;
 
 }
 
 return(List::create(_["b"] = b0, _["e"] = eps, _["llk"] = llk1)) ;
  
}



/*** R
#test <- NRb(as.matrix(Y),form,fam,1,5,maxit = 100)
*/

