// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <miscfun.hpp>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::cube DY(DataFrame& Y, List& ghQ, List& b, StringVector& fam){
 arma::mat Yt = DFtoNM(Y);
 List Z = ghQ["out"];
 arma::mat gP = DFtoNM(ghQ["points"]);
 arma::cube fyz(Yt.n_rows, gP.n_rows, Yt.n_cols, fill::value(0.0));
 NumericVector mu ;
 NumericVector sigma ;
 
 for(unsigned int z = 0; z < gP.n_rows; z++){
  for(unsigned int i = 0; i < Yt.n_cols; i++){
   NumericVector Yi = Rcpp::wrap(Yt.col(i));
   if(fam(i) == "normal"){
    arma::mat Zmu = DFtoNM(Z["mu"]) ;
    arma::mat bmu = DFtoNM(b["mu"]) ;
    mu = Zmu.row(z) * bmu.row(i).t() ;
    arma::mat Zsg = DFtoNM(Z["sigma"]) ;
    arma::mat bsg = DFtoNM(b["sigma"]) ;
    sigma = exp(Zsg.row(z) * bsg.row(i).t()) ;
    for(unsigned int j = 0; j < Yi.size(); j++){
     fyz.slice(i).row(j).col(z) = R::dnorm(Yi(j), as<double>(mu), as<double>(sigma), true); }
    continue;
   }
   if(fam(i) == "binomial"){
    arma::mat Zmu = DFtoNM(Z["mu"]) ;
    arma::mat bmu = DFtoNM(b["mu"]) ;
    mu = cprobs(Zmu.row(z) * bmu.row(i).t()) ;
    for(unsigned int j = 0; j < Yi.size(); j++){
     fyz.slice(i).row(j).col(z) = R::dbinom(Yi(j), 1, as<double>(mu), true); }
    continue;
   }
   if(fam(i) == "lognormal"){
    arma::mat Zmu = DFtoNM(Z["mu"]) ;
    arma::mat bmu = DFtoNM(b["mu"]) ;
    mu = Zmu.row(z) * bmu.row(i).t() ;
    arma::mat Zsg = DFtoNM(Z["sigma"]) ;
    arma::mat bsg = DFtoNM(b["sigma"]) ;
    sigma = exp(Zsg.row(z) * bsg.row(i).t()) ;
    for(unsigned int j = 0; j < Yi.size(); j++){
     fyz.slice(i).row(j).col(z) = R::dlnorm(Yi(j), as<double>(mu), as<double>(sigma), true); }
    continue;
   }
   if(fam(i) == "poisson"){
    arma::mat Zmu = DFtoNM(Z["mu"]) ;
    arma::mat bmu = DFtoNM(b["mu"]) ;
    mu = exp(Zmu.row(z) * bmu.row(i).t()) ;
    for(unsigned int j = 0; j < Yi.size(); j++){
     fyz.slice(i).row(j).col(z) = R::dpois(Yi(j), as<double>(mu), true); }
    continue;
   }
   if(fam(i) == "gamma"){
    arma::mat Zmu = DFtoNM(Z["mu"]) ;
    arma::mat bmu = DFtoNM(b["mu"]) ;
    mu = exp(Zmu.row(z) * bmu.row(i).t()) ;
    arma::mat Zsg = DFtoNM(Z["sigma"]) ;
    arma::mat bsg = DFtoNM(b["sigma"]) ;
    sigma = exp(Zsg.row(z) * bsg.row(i).t()) ;
    for(unsigned int j = 0; j < Yi.size(); j++){
     fyz.slice(i).row(j).col(z) = R::dgamma(Yi(j), as<double>(mu), as<double>(sigma), true); }
    continue;
   }
   if(fam(i) == "ZIpoisson"){
    arma::mat Zmu = DFtoNM(Z["mu"]) ;
    arma::mat bmu = DFtoNM(b["mu"]) ;
    mu = exp(Zmu.row(z) * bmu.row(i).t());
    arma::mat Zsg = DFtoNM(Z["sigma"]) ;
    arma::mat bsg = DFtoNM(b["sigma"]) ;
    sigma = cprobs(Zsg.row(z) * bsg.row(i).t()) ;
    for(unsigned int j = 0; j < Yi.size(); j++){
     fyz.slice(i).row(j).col(z) = dZIPo(Yi(j), as<double>(mu), as<double>(sigma), true); }
    continue;
   }
  }
 }
 return fyz ;
}


// Rcpp::List cDVL(DataFrame& Y, List& ghQ, List& b, StringVector& fam, StringVector& info){
//   Rcpp::List dvY = Rcpp::List::create( dvY ) ;
//   return(dvY);
// }

 
// //[[Rcpp::export]]
// List EFALasso(arma::mat loadingm, arma::mat psim, arma::mat scov, int indicatorn, int factorn, int maxstep, arma::vec penaltyvec, int penaltyvecn, double tol, arma::mat loadpenweight) {
// 
//   int iteration;
//   mat modelcov;
//   mat matrixm;
//   mat matrixam;
//   mat bmi;
//   mat s1;
//   double sumtemp = 0.0;
//   double thetatilde = 0.0;
//   double evasign = 0.0;
//   //double signcoef = 0.0;
//   double signtilde = 0.0;
//   mat lambdaupdate = zeros<mat>(indicatorn,factorn);
//   mat psiupdate = zeros<mat>(indicatorn,indicatorn);
//   mat lambdasub;
//   mat mattemp;
//   double mattodouble;
//   mat loadingdif;
//   mat psidif;
//   double loadingmaxdif = 1.0;
//   double psimaxdif = 1.0;
//   int tolproceed = 1;
//   double penal;
//   double penalstar;
//   mat lambdahatold;
//   mat psihatold;
//   mat loadinghatnow;
//   mat psihatnow;
//   mat lambdaout = zeros<mat>(indicatorn*factorn,penaltyvecn);
//   mat psiout = zeros<mat>(indicatorn,penaltyvecn);
//   mat penaltyvalue = zeros<mat>(1,penaltyvecn);
//   vec itercheck = zeros<vec>(penaltyvecn);
//   mat iterdiff = zeros<mat>(2,penaltyvecn);
//   
//   for(int ipenal=0; ipenal<penaltyvecn; ipenal++){
//   	    penal = penaltyvec(ipenal);
//     
//         if( ipenal==0 ){
//             loadinghatnow = loadingm;
//             psihatnow = psim;
//         } 
//         
//         iteration = 1;
//         tolproceed = 1;
//   	
//   	    while( iteration <= maxstep && tolproceed == 1){
//   	    	lambdahatold = loadinghatnow;
//   	    	psihatold = psihatnow;
//   	    	
//             modelcov = loadinghatnow*trans(loadinghatnow) + psihatnow;
//             matrixm = trans(loadinghatnow)*inv(psihatnow)*loadinghatnow + eye( factorn, factorn ) ;
//             matrixam = inv(matrixm) + inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov*inv(psihatnow)*loadinghatnow*inv(matrixm);
//     
//             lambdaupdate = zeros<mat>(indicatorn,factorn);
//             psiupdate = zeros<mat>(indicatorn,indicatorn);
//     
//             for (int i=0; i<indicatorn; i++) {
//                 bmi = inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov.cols(i,i);
//       
//                 for (int j=0; j<factorn; j++) {
//                     sumtemp = 0.0;
//                     for (int row=0; row<factorn; row++) {
//                         if(row!=j) sumtemp = sumtemp + matrixam(row,j)*loadinghatnow(i,row);
//                     }
//                     thetatilde = ( bmi(j,0) - sumtemp )/matrixam(j,j);
//                     penalstar = loadpenweight(i,j)*psihatnow(i,i)*penal/matrixam(j,j);
//                     evasign = std::abs(thetatilde) - penalstar;
//                     signtilde = 0.0;
//                     if( thetatilde > 0.0 ) signtilde = 1.0;
//                     if( thetatilde < 0.0 ) signtilde = -1.0;
//                     if( evasign > 0.0 ) {
//                        lambdaupdate(i,j) = signtilde*evasign;
//                     }
//                     else {
//                         lambdaupdate(i,j) = 0.0;
//                     }
//         
//                 }
//                 lambdasub = lambdaupdate.rows(i,i);
//                 mattemp = -2.0*lambdasub*bmi+lambdasub*matrixam*trans(lambdasub);
//                 mattodouble = mattemp(0,0);
//                 //mattodouble = as_scalar(mattemp);
//                 psiupdate(i,i) =  scov(i,i) + mattodouble;
//             }
//     
//             loadingdif = abs(loadinghatnow - lambdaupdate);
//             psidif = abs(psiupdate - psihatnow);
//             loadingmaxdif = loadingdif.max();
//             psimaxdif = psidif.max();
//     
//             if( loadingmaxdif < tol && psimaxdif < tol) tolproceed = 0L;
//             
// 			loadinghatnow = lambdaupdate;
//             psihatnow = psiupdate;
//     
//             iteration = iteration + 1;
//         }
//         
//         lambdaout.col(ipenal) = vectorise(loadinghatnow);
//         psiout.col(ipenal) = diagvec(psihatnow);
//         penaltyvalue(0,ipenal) = penal;
//         itercheck(ipenal) = iteration-1;
//         iterdiff(0,ipenal) =  loadingmaxdif;
//         iterdiff(1,ipenal) =  psimaxdif;
//     }
//   
//   
//   return Rcpp::List::create(
//     Rcpp::Named("loading") = lambdaout,
//     Rcpp::Named("unique") = psiout,
//     Rcpp::Named("lam") = penaltyvalue,
//     Rcpp::Named("steps") = itercheck,
//     Rcpp::Named("diff") = iterdiff
//   );
//   
// }
// 
