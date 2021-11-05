// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace arma;



//[[Rcpp::export]]
List EFALasso(arma::mat loadingm, arma::mat psim, arma::mat scov, int indicatorn, int factorn, int maxstep, arma::vec penaltyvec, int penaltyvecn, double tol, arma::mat loadpenweight) {

  int iteration;
  mat modelcov;
  mat matrixm;
  mat matrixam;
  mat bmi;
  mat s1;
  double sumtemp = 0.0;
  double thetatilde = 0.0;
  double evasign = 0.0;
  //double signcoef = 0.0;
  double signtilde = 0.0;
  mat lambdaupdate = zeros<mat>(indicatorn,factorn);
  mat psiupdate = zeros<mat>(indicatorn,indicatorn);
  mat lambdasub;
  mat mattemp;
  double mattodouble;
  mat loadingdif;
  mat psidif;
  double loadingmaxdif = 1.0;
  double psimaxdif = 1.0;
  int tolproceed = 1;
  double penal;
  double penalstar;
  mat lambdahatold;
  mat psihatold;
  mat loadinghatnow;
  mat psihatnow;
  mat lambdaout = zeros<mat>(indicatorn*factorn,penaltyvecn);
  mat psiout = zeros<mat>(indicatorn,penaltyvecn);
  mat penaltyvalue = zeros<mat>(1,penaltyvecn);
  vec itercheck = zeros<vec>(penaltyvecn);
  mat iterdiff = zeros<mat>(2,penaltyvecn);
  
  for(int ipenal=0; ipenal<penaltyvecn; ipenal++){
  	    penal = penaltyvec(ipenal);
    
        if( ipenal==0 ){
            loadinghatnow = loadingm;
            psihatnow = psim;
        } 
        
        iteration = 1;
        tolproceed = 1;
  	
  	    while( iteration <= maxstep && tolproceed == 1){
  	    	lambdahatold = loadinghatnow;
  	    	psihatold = psihatnow;
  	    	
            modelcov = loadinghatnow*trans(loadinghatnow) + psihatnow;
            matrixm = trans(loadinghatnow)*inv(psihatnow)*loadinghatnow + eye( factorn, factorn ) ;
            matrixam = inv(matrixm) + inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov*inv(psihatnow)*loadinghatnow*inv(matrixm);
    
            lambdaupdate = zeros<mat>(indicatorn,factorn);
            psiupdate = zeros<mat>(indicatorn,indicatorn);
    
            for (int i=0; i<indicatorn; i++) {
                bmi = inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov.cols(i,i);
      
                for (int j=0; j<factorn; j++) {
                    sumtemp = 0.0;
                    for (int row=0; row<factorn; row++) {
                        if(row!=j) sumtemp = sumtemp + matrixam(row,j)*loadinghatnow(i,row);
                    }
                    thetatilde = ( bmi(j,0) - sumtemp )/matrixam(j,j);
                    penalstar = loadpenweight(i,j)*psihatnow(i,i)*penal/matrixam(j,j);
                    evasign = std::abs(thetatilde) - penalstar;
                    signtilde = 0.0;
                    if( thetatilde > 0.0 ) signtilde = 1.0;
                    if( thetatilde < 0.0 ) signtilde = -1.0;
                    if( evasign > 0.0 ) {
                       lambdaupdate(i,j) = signtilde*evasign;
                    }
                    else {
                        lambdaupdate(i,j) = 0.0;
                    }
        
                }
                lambdasub = lambdaupdate.rows(i,i);
                mattemp = -2.0*lambdasub*bmi+lambdasub*matrixam*trans(lambdasub);
                mattodouble = mattemp(0,0);
                //mattodouble = as_scalar(mattemp);
                psiupdate(i,i) =  scov(i,i) + mattodouble;
            }
    
            loadingdif = abs(loadinghatnow - lambdaupdate);
            psidif = abs(psiupdate - psihatnow);
            loadingmaxdif = loadingdif.max();
            psimaxdif = psidif.max();
    
            if( loadingmaxdif < tol && psimaxdif < tol) tolproceed = 0L;
            
			loadinghatnow = lambdaupdate;
            psihatnow = psiupdate;
    
            iteration = iteration + 1;
        }
        
        lambdaout.col(ipenal) = vectorise(loadinghatnow);
        psiout.col(ipenal) = diagvec(psihatnow);
        penaltyvalue(0,ipenal) = penal;
        itercheck(ipenal) = iteration-1;
        iterdiff(0,ipenal) =  loadingmaxdif;
        iterdiff(1,ipenal) =  psimaxdif;
    }
  
  
  return Rcpp::List::create(
    Rcpp::Named("loading") = lambdaout,
    Rcpp::Named("unique") = psiout,
    Rcpp::Named("lam") = penaltyvalue,
    Rcpp::Named("steps") = itercheck,
    Rcpp::Named("diff") = iterdiff
  );
  
}


//[[Rcpp::export]]
List EFALassoCorr(arma::mat loadingm, arma::mat psim, arma::mat scov, int indicatorn, int factorn, int maxstep, arma::vec penaltyvec, int penaltyvecn, double tol, arma::mat loadpenweight) {

  int iteration;
  mat modelcov;
  mat matrixm;
  mat matrixam;
  mat bmi;
  mat s1;
  double sumtemp = 0.0;
  double thetatilde = 0.0;
  double evasign = 0.0;
  //double signcoef = 0.0;
  double signtilde = 0.0;
  mat lambdaupdate = zeros<mat>(indicatorn,factorn);
  mat psiupdate = zeros<mat>(indicatorn,indicatorn);
  mat lambdasub;
  mat mattemp;
  double mattodouble;
  mat loadingdif;
  mat psidif;
  double loadingmaxdif = 1.0;
  double psimaxdif = 1.0;
  int tolproceed = 1;
  double penal;
  double penalstar;
  int nonsingular;
  mat lambdahatold;
  mat psihatold;
  mat loadinghatnow;
  mat psihatnow;
  mat lambdaout = zeros<mat>(indicatorn*factorn,penaltyvecn);
  mat psiout = zeros<mat>(indicatorn,penaltyvecn);
  mat penaltyvalue = zeros<mat>(1,penaltyvecn);
  vec itercheck = zeros<vec>(penaltyvecn);
  mat iterdiff = zeros<mat>(2,penaltyvecn);
  
  for(int ipenal=0; ipenal<penaltyvecn; ipenal++){
  	    penal = penaltyvec(ipenal);
    
        if( ipenal==0 ){
            loadinghatnow = loadingm;
            psihatnow = psim;
        } 
        
        iteration = 1;
        tolproceed = 1;
  	
  	    while( iteration <= maxstep && tolproceed == 1){
  	    	lambdahatold = loadinghatnow;
  	    	psihatold = psihatnow;
  	    	
  	    	nonsingular = 0L;
            for(int ni=0; ni<indicatorn; ni++) {
                if( psihatnow(ni,ni)>0.000001 ) {
                	nonsingular = nonsingular + 1;
				}
 			}
 				
 			if( nonsingular != indicatorn){
 				loadinghatnow = loadingm;
 				psihatnow = psim;
 				iteration = maxstep + 2;
 				loadingmaxdif = 1.0;
 				psimaxdif = 1.0;
 				break;
			}
				 
            modelcov = loadinghatnow*trans(loadinghatnow) + psihatnow;
            matrixm = trans(loadinghatnow)*inv(psihatnow)*loadinghatnow + eye( factorn, factorn ) ;
            matrixam = inv(matrixm) + inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov*inv(psihatnow)*loadinghatnow*inv(matrixm);
    
            lambdaupdate = zeros<mat>(indicatorn,factorn);
            psiupdate = zeros<mat>(indicatorn,indicatorn);
    
            for (int i=0; i<indicatorn; i++) {
                bmi = inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov.cols(i,i);
      
                for (int j=0; j<factorn; j++) {
                    sumtemp = 0.0;
                    for (int row=0; row<factorn; row++) {
                        if(row!=j) sumtemp = sumtemp + matrixam(row,j)*loadinghatnow(i,row);
                    }
                    thetatilde = ( bmi(j,0) - sumtemp )/matrixam(j,j);
                    penalstar = loadpenweight(i,j)*psihatnow(i,i)*penal/matrixam(j,j);
                    evasign = std::abs(thetatilde) - penalstar;
                    signtilde = 0.0;
                    if( thetatilde > 0.0 ) signtilde = 1.0;
                    if( thetatilde < 0.0 ) signtilde = -1.0;
                    if( evasign > 0.0 ) {
                       lambdaupdate(i,j) = signtilde*evasign;
                    }
                    else {
                        lambdaupdate(i,j) = 0.0;
                    }
        
                }
                lambdasub = lambdaupdate.rows(i,i);
                mattemp = lambdasub*trans(lambdasub);
                //mattodouble = mattemp(0,0);
                mattodouble = as_scalar(mattemp);
                //psiupdate(i,i) =  scov(i,i) + mattodouble;
                psiupdate(i,i) =  1.0 - mattodouble;
            }
    
            loadingdif = abs(loadinghatnow - lambdaupdate);
            psidif = abs(psiupdate - psihatnow);
            loadingmaxdif = loadingdif.max();
            psimaxdif = psidif.max();
    
            if( loadingmaxdif < tol && psimaxdif < tol) tolproceed = 0L;
            
			loadinghatnow = lambdaupdate;
            psihatnow = psiupdate;
    
            iteration = iteration + 1;
        }
        
        lambdaout.col(ipenal) = vectorise(loadinghatnow);
        psiout.col(ipenal) = diagvec(psihatnow);
        penaltyvalue(0,ipenal) = penal;
        itercheck(ipenal) = iteration-1;
        iterdiff(0,ipenal) =  loadingmaxdif;
        iterdiff(1,ipenal) =  psimaxdif;
    }
  
  
  return Rcpp::List::create(
    Rcpp::Named("loading") = lambdaout,
    Rcpp::Named("unique") = psiout,
    Rcpp::Named("lam") = penaltyvalue,
    Rcpp::Named("steps") = itercheck,
    Rcpp::Named("diff") = iterdiff
  );
  
}



//[[Rcpp::export]]
List LassoNoInt(arma::mat y, arma::mat x, int n, int p, arma::vec xnorm, arma::mat betainitial, arma::vec lamvec, 
int lamn, int maxstep, double tol, arma::vec penweight) {

mat betachange = betainitial;
mat ytilde;
double betatilde = 0.0;
double sumtemp;
double lamstar = 0.0;
double absbeta = 0.0;
double signbeta = 0.0;
int iteration;
mat betaold;
mat betadif;
double maxdif;
double lam;
mat betaall;
mat cbbetaall;
vec itercheck = zeros(lamn);
mat lamgam = zeros(lamn,1);
int count = 0;

for (int lamloop=0; lamloop<lamn; lamloop++) {
    lam = lamvec(lamloop);

        mat takeout = zeros(n,1);
        iteration = 1;
        maxdif = 1.0;
        ytilde = y - x*betachange;

        while( iteration <= maxstep && maxdif > tol){
            betaold = betachange;
            for (int k=0; k<p; k++) {
                sumtemp = 0.0;
                for (int m=0; m<n; m++) {
                    ytilde(m,0) = ytilde(m,0) + x(m,k)*betachange(k,0) - takeout(m,0);
                    sumtemp = sumtemp + ytilde(m,0)*x(m,k);
                }
                betatilde = sumtemp/xnorm(k);
		
		        if( penweight(k)==0.0 ){
		            betachange(k,0) = betatilde;
		        }
		        else {
		        	lamstar = penweight(k)*lam/xnorm(k);
		            absbeta = std::abs(betatilde);
  
		            if( absbeta <= lamstar ){
                        betachange(k,0) = 0.0;
		            } 
		            else {
		                signbeta = 0.0;
			            if( betatilde > 0.0 ){
			                signbeta = 1.0;
			            } 
			            else {
			                signbeta = -1.0;
			            }
			            betachange(k,0) = signbeta*( absbeta-lamstar );
		            } 
		            
		        }
		
		        for (int m=0; m<n; m++) {
		            takeout(m,0) = x(m,k)*betachange(k,0);
		        }
	        }  
	
	        betadif = abs(betaold - betachange);
	        maxdif = betadif.max();
	        iteration = iteration + 1;
        }

        if( lamloop==0 ){
            betaall = betachange;
        } else {
            cbbetaall = join_rows(betaall, betachange);
            betaall = cbbetaall;
        }

        itercheck(count) = iteration;
        lamgam(count,0) = lam;
        count = count + 1;
}

return Rcpp::List::create(
        Rcpp::Named("parall") = betaall,
        Rcpp::Named("beta1beta2") = lamgam,
        Rcpp::Named("steps") = itercheck
);

}


//[[Rcpp::export]]
List EFASCAD(arma::mat loadingm, arma::mat psim, arma::mat scov, int indicatorn, int factorn, int maxstep, arma::vec lamvec, int lamvecn, arma::vec gamvec, int gamvecn, double tol, arma::mat loadpenweight, int LLA) {

  int iteration;
  mat modelcov;
  mat matrixm;
  mat matrixam;
  mat bmi;
  mat s1;
  double sumtemp = 0.0;
  double thetatilde = 0.0;
  //double evasign = 0.0;
  //double signcoef = 0.0;
  double signtilde = 0.0;
  mat lambdaupdate = zeros<mat>(indicatorn,factorn);
  mat psiupdate = zeros<mat>(indicatorn,indicatorn);
  mat lambdasub;
  mat mattemp;
  double mattodouble;
  mat loadingdif;
  mat psidif;
  double loadingmaxdif = 1.0;
  double psimaxdif = 1.0;
  int tolproceed = 1;
  double penalvaluelam;
  double penalvaluegam;
  //double penalvaluelamgam;
  double penallamstar;
  //double penalgamstar;
  double posiloading;
  double scadweight;
  double temp;
  double posithetatilde;
  mat lambdahatold;
  mat psihatold;
  mat loadinghatnow;
  mat psihatnow;
  mat loadinghatwarm;
  mat psihatwarm;
  mat lambdaout = zeros<mat>(indicatorn*factorn,lamvecn*gamvecn);
  mat psiout = zeros<mat>(indicatorn,lamvecn*gamvecn);
  mat penaltyvalue = zeros<mat>(2,lamvecn*gamvecn);
  vec itercheck = zeros<vec>(lamvecn*gamvecn);
  mat iterdiff = zeros<mat>(2,lamvecn*gamvecn);
  
  int count = 0;
  for(int ipenallam=0; ipenallam<lamvecn; ipenallam++){
  	    penalvaluelam = lamvec(ipenallam);
    
        if( ipenallam==0 ){
            loadinghatnow = loadingm;
            psihatnow = psim;
        } 
		else {
			loadinghatnow = loadinghatwarm;
            psihatnow = psihatwarm;
		}
        
        for(int ipenalgam=0; ipenalgam<gamvecn; ipenalgam++) {
        	penalvaluegam = gamvec(ipenalgam);
        	
        	iteration = 1;
            tolproceed = 1;
  	
  	        while( iteration <= maxstep && tolproceed == 1){
  	    	    lambdahatold = loadinghatnow;
  	    	    psihatold = psihatnow;
  	    	
                modelcov = loadinghatnow*trans(loadinghatnow) + psihatnow;
                matrixm = trans(loadinghatnow)*inv(psihatnow)*loadinghatnow + eye( factorn, factorn ) ;
                matrixam = inv(matrixm) + inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov*inv(psihatnow)*loadinghatnow*inv(matrixm);
    
                lambdaupdate = zeros<mat>(indicatorn,factorn);
                psiupdate = zeros<mat>(indicatorn,indicatorn);
    
                for (int i=0; i<indicatorn; i++) {
                    bmi = inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov.cols(i,i);
      
                    for (int j=0; j<factorn; j++) {
                        sumtemp = 0.0;
                        for (int row=0; row<factorn; row++) {
                            if(row!=j) sumtemp = sumtemp + matrixam(row,j)*loadinghatnow(i,row);
                        }
                        thetatilde = ( bmi(j,0) - sumtemp )/matrixam(j,j);
                        
                        if( LLA==1L ){
						
                        	posiloading = std::abs( loadingm(i,j) );
                        	if( posiloading<=penalvaluelam ){
                        		scadweight = 1.0;
							} 
							else {
								temp = penalvaluegam*penalvaluelam - posiloading;
								if( temp<=0.0 ){
									scadweight = 0.0;
								}
								else {
									scadweight = temp/( (penalvaluegam-1)*penalvaluelam );
								}
							
							}
                        
                        	penallamstar = scadweight*loadpenweight(i,j)*psihatnow(i,i)*penalvaluelam/matrixam(j,j);
                    	   	 //penalgamstar = 2*penalvaluegam*psihatnow(i,i)*loadpenweight(i,j)/matrixam(j,j);
                   		     //thetatildetilde = thetatilde/( 1+penalgamstar );
							//penalvaluelamgam = penallamstar/( 1+penalgamstar );
                        	posithetatilde = std::abs(thetatilde) - penallamstar;
                        
                        	if( posithetatilde <= 0.0 ){
                        		lambdaupdate(i,j) = 0.0;
							}
							else {
                        		signtilde = 0.0;
                        		if( thetatilde > 0.0 ) signtilde = 1.0;
                        		if( thetatilde < 0.0 ) signtilde = -1.0;
                        	
                        		lambdaupdate(i,j) = signtilde*posithetatilde;
							}
						}
						
                    }
                    lambdasub = lambdaupdate.rows(i,i);
                    mattemp = -2.0*lambdasub*bmi+lambdasub*matrixam*trans(lambdasub);
                    //mattodouble = mattemp(0,0);
                    mattodouble = as_scalar(mattemp);
                    psiupdate(i,i) =  scov(i,i) + mattodouble;
                }
    
                loadingdif = abs(loadinghatnow - lambdaupdate);
                psidif = abs(psiupdate - psihatnow);
                loadingmaxdif = loadingdif.max();
                psimaxdif = psidif.max();
    
                if( loadingmaxdif < tol && psimaxdif < tol) tolproceed = 0L;
            
			    loadinghatnow = lambdaupdate;
                psihatnow = psiupdate;
    
                iteration = iteration + 1;
            }
            
            if( ipenalgam==0 ){
            	loadinghatwarm = loadinghatnow ;
                psihatwarm = psihatnow;
			}
        
            lambdaout.col(count) = vectorise(loadinghatnow);
            psiout.col(count) = diagvec(psihatnow);
            penaltyvalue(0,count) = penalvaluelam;
            penaltyvalue(1,count) = penalvaluegam;
            itercheck(count) = iteration-1;
            iterdiff(0,count) =  loadingmaxdif;
            iterdiff(1,count) =  psimaxdif;
            
            count = count + 1;
		}
        
   }
  
  
  return Rcpp::List::create(
    Rcpp::Named("loading") = lambdaout,
    Rcpp::Named("unique") = psiout,
    Rcpp::Named("lamgam") = penaltyvalue,
    Rcpp::Named("steps") = itercheck,
    Rcpp::Named("diff") = iterdiff
  );
  
}



//[[Rcpp::export]]
List EFASCADCorr(arma::mat loadingm, arma::mat psim, arma::mat scov, int indicatorn, int factorn, int maxstep, arma::vec lamvec, int lamvecn, arma::vec gamvec, int gamvecn, double tol, arma::mat loadpenweight, int LLA) {

  int iteration;
  mat modelcov;
  mat matrixm;
  mat matrixam;
  mat bmi;
  mat s1;
  double sumtemp = 0.0;
  double thetatilde = 0.0;
  //double evasign = 0.0;
  //double signcoef = 0.0;
  double signtilde = 0.0;
  mat lambdaupdate = zeros<mat>(indicatorn,factorn);
  mat psiupdate = zeros<mat>(indicatorn,indicatorn);
  mat lambdasub;
  mat mattemp;
  double mattodouble;
  mat loadingdif;
  mat psidif;
  double loadingmaxdif = 1.0;
  double psimaxdif = 1.0;
  int tolproceed = 1;
  double penalvaluelam;
  double penalvaluegam;
  //double penalvaluelamgam;
  double penallamstar;
  //double penalgamstar;
  double posiloading;
  double scadweight;
  double temp;
  double posithetatilde;
  int nonsingular;
  mat lambdahatold;
  mat psihatold;
  mat loadinghatnow;
  mat psihatnow;
  mat loadinghatwarm;
  mat psihatwarm;
  mat lambdaout = zeros<mat>(indicatorn*factorn,lamvecn*gamvecn);
  mat psiout = zeros<mat>(indicatorn,lamvecn*gamvecn);
  mat penaltyvalue = zeros<mat>(2,lamvecn*gamvecn);
  vec itercheck = zeros<vec>(lamvecn*gamvecn);
  mat iterdiff = zeros<mat>(2,lamvecn*gamvecn);
  
  int count = 0;
  for(int ipenallam=0; ipenallam<lamvecn; ipenallam++){
  	    penalvaluelam = lamvec(ipenallam);
    
        if( ipenallam==0 ){
            loadinghatnow = loadingm;
            psihatnow = psim;
        } 
		else {
			loadinghatnow = loadinghatwarm;
            psihatnow = psihatwarm;
		}
        
        for(int ipenalgam=0; ipenalgam<gamvecn; ipenalgam++) {
        	penalvaluegam = gamvec(ipenalgam);
        	
        	iteration = 1;
            tolproceed = 1;
  	
  	        while( iteration <= maxstep && tolproceed == 1){
  	    	    lambdahatold = loadinghatnow;
  	    	    psihatold = psihatnow;
  	    	
                modelcov = loadinghatnow*trans(loadinghatnow) + psihatnow;
                nonsingular = 0L;
                for(int ni=0; ni<indicatorn; ni++) {
                	if( psihatnow(ni,ni)>0.000001 ) {
                		nonsingular = nonsingular + 1;
					}
 				}
 				
 				if( nonsingular != indicatorn){
 					if( ipenalgam !=0 ) {
 						loadinghatnow = loadinghatwarm;
 					    psihatnow = psihatwarm;
					 } else {
					 	loadinghatnow = loadingm;
 				    	psihatnow = psim;
					 }
 					iteration = maxstep + 2;
 					loadingmaxdif = 1.0;
 					psimaxdif = 1.0;
 					break;
				 }
 				
                matrixm = trans(loadinghatnow)*inv(psihatnow)*loadinghatnow + eye( factorn, factorn ) ;
                matrixam = inv(matrixm) + inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov*inv(psihatnow)*loadinghatnow*inv(matrixm);
    
                lambdaupdate = zeros<mat>(indicatorn,factorn);
                psiupdate = zeros<mat>(indicatorn,indicatorn);
    
                for (int i=0; i<indicatorn; i++) {
                    bmi = inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov.cols(i,i);
      
                    for (int j=0; j<factorn; j++) {
                        sumtemp = 0.0;
                        for (int row=0; row<factorn; row++) {
                            if(row!=j) sumtemp = sumtemp + matrixam(row,j)*loadinghatnow(i,row);
                        }
                        thetatilde = ( bmi(j,0) - sumtemp )/matrixam(j,j);
                        
                        if( LLA==1L ){
						
                        	posiloading = std::abs( loadingm(i,j) );
                        	if( posiloading<=penalvaluelam ){
                        		scadweight = 1.0;
							} 
							else {
								temp = penalvaluegam*penalvaluelam - posiloading;
								if( temp<=0.0 ){
									scadweight = 0.0;
								}
								else {
									scadweight = temp/( (penalvaluegam-1)*penalvaluelam );
								}
							
							}
                        
                        	penallamstar = scadweight*loadpenweight(i,j)*psihatnow(i,i)*penalvaluelam/matrixam(j,j);
                    	   	 //penalgamstar = 2*penalvaluegam*psihatnow(i,i)*loadpenweight(i,j)/matrixam(j,j);
                   		     //thetatildetilde = thetatilde/( 1+penalgamstar );
							//penalvaluelamgam = penallamstar/( 1+penalgamstar );
                        	posithetatilde = std::abs(thetatilde) - penallamstar;
                        
                        	if( posithetatilde <= 0.0 ){
                        		lambdaupdate(i,j) = 0.0;
							}
							else {
                        		signtilde = 0.0;
                        		if( thetatilde > 0.0 ) signtilde = 1.0;
                        		if( thetatilde < 0.0 ) signtilde = -1.0;
                        	
                        		lambdaupdate(i,j) = signtilde*posithetatilde;
							}
						}
						
                    }
                    lambdasub = lambdaupdate.rows(i,i);
                    mattemp = lambdasub*trans(lambdasub);
                    //mattodouble = mattemp(0,0);
                    mattodouble = as_scalar(mattemp);
                    //psiupdate(i,i) =  scov(i,i) + mattodouble;
                    psiupdate(i,i) =  1.0 - mattodouble;
                }
    
                loadingdif = abs(loadinghatnow - lambdaupdate);
                psidif = abs(psiupdate - psihatnow);
                loadingmaxdif = loadingdif.max();
                psimaxdif = psidif.max();
    
                if( loadingmaxdif < tol && psimaxdif < tol) tolproceed = 0L;
            
			    loadinghatnow = lambdaupdate;
                psihatnow = psiupdate;
    
                iteration = iteration + 1;
            }
            
            if( ipenalgam==0 ){
            	loadinghatwarm = loadinghatnow ;
                psihatwarm = psihatnow;
			}
        
            lambdaout.col(count) = vectorise(loadinghatnow);
            psiout.col(count) = diagvec(psihatnow);
            penaltyvalue(0,count) = penalvaluelam;
            penaltyvalue(1,count) = penalvaluegam;
            itercheck(count) = iteration-1;
            iterdiff(0,count) =  loadingmaxdif;
            iterdiff(1,count) =  psimaxdif;
            
            count = count + 1;
		}
        
   }
  
  
  return Rcpp::List::create(
    Rcpp::Named("loading") = lambdaout,
    Rcpp::Named("unique") = psiout,
    Rcpp::Named("lamgam") = penaltyvalue,
    Rcpp::Named("steps") = itercheck,
    Rcpp::Named("diff") = iterdiff
  );
  
}



//[[Rcpp::export]]
List EFAMCP(arma::mat loadingm, arma::mat psim, arma::mat scov, int indicatorn, int factorn, int maxstep, arma::vec lamvec, int lamvecn, arma::vec gamvec, int gamvecn, double tol, arma::mat loadpenweight) {

  int iteration;
  mat modelcov;
  mat matrixm;
  mat matrixam;
  mat bmi;
  mat s1;
  double sumtemp = 0.0;
  double thetatilde = 0.0;
  //double evasign = 0.0;
  //double signcoef = 0.0;
  double signtilde = 0.0;
  mat lambdaupdate = zeros<mat>(indicatorn,factorn);
  mat psiupdate = zeros<mat>(indicatorn,indicatorn);
  mat lambdasub;
  mat mattemp;
  double mattodouble;
  mat loadingdif;
  mat psidif;
  double loadingmaxdif = 1.0;
  double psimaxdif = 1.0;
  int tolproceed = 1;
  double penalvaluelam;
  double penalvaluegam;
  double penalvaluelamgam;
  double penallamstar;
  double penalgamstar;
  double posithetatilde;
  mat lambdahatold;
  mat psihatold;
  mat loadinghatnow;
  mat psihatnow;
  mat loadinghatwarm;
  mat psihatwarm;
  mat lambdaout = zeros<mat>(indicatorn*factorn,lamvecn*gamvecn);
  mat psiout = zeros<mat>(indicatorn,lamvecn*gamvecn);
  mat penaltyvalue = zeros<mat>(2,lamvecn*gamvecn);
  vec itercheck = zeros<vec>(lamvecn*gamvecn);
  mat iterdiff = zeros<mat>(2,lamvecn*gamvecn);
  
  int count = 0;
  for(int ipenallam=0; ipenallam<lamvecn; ipenallam++){
  	    penalvaluelam = lamvec(ipenallam);
    
        if( ipenallam==0 ){
            loadinghatnow = loadingm;
            psihatnow = psim;
        } 
		else {
			loadinghatnow = loadinghatwarm;
            psihatnow = psihatwarm;
		}
        
        for(int ipenalgam=0; ipenalgam<gamvecn; ipenalgam++) {
        	penalvaluegam = gamvec(ipenalgam);
        	penalvaluelamgam = penalvaluelam*penalvaluegam;
        	
        	iteration = 1;
            tolproceed = 1;
  	
  	        while( iteration <= maxstep && tolproceed == 1){
  	    	    lambdahatold = loadinghatnow;
  	    	    psihatold = psihatnow;
  	    	
                modelcov = loadinghatnow*trans(loadinghatnow) + psihatnow;
                matrixm = trans(loadinghatnow)*inv(psihatnow)*loadinghatnow + eye( factorn, factorn ) ;
                matrixam = inv(matrixm) + inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov*inv(psihatnow)*loadinghatnow*inv(matrixm);
    
                lambdaupdate = zeros<mat>(indicatorn,factorn);
                psiupdate = zeros<mat>(indicatorn,indicatorn);
    
                for (int i=0; i<indicatorn; i++) {
                    bmi = inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov.cols(i,i);
      
                    for (int j=0; j<factorn; j++) {
                        sumtemp = 0.0;
                        for (int row=0; row<factorn; row++) {
                            if(row!=j) sumtemp = sumtemp + matrixam(row,j)*loadinghatnow(i,row);
                        }
                        thetatilde = ( bmi(j,0) - sumtemp )/matrixam(j,j);
                        penallamstar = loadpenweight(i,j)*psihatnow(i,i)*penalvaluelam/matrixam(j,j);
                        penalgamstar = matrixam(j,j)*penalvaluegam/( loadpenweight(i,j)*psihatnow(i,i) );
                        posithetatilde = std::abs(thetatilde);
                        
                        if( posithetatilde <= penallamstar){
                        	lambdaupdate(i,j) = 0.0;
						}
						else if( posithetatilde <= penalvaluelamgam ){
                        	signtilde = 0.0;
                        	if( thetatilde > 0.0 ) signtilde = 1.0;
                        	if( thetatilde < 0.0 ) signtilde = -1.0;
                        	
                        	lambdaupdate(i,j) = signtilde*( posithetatilde-penallamstar )/(1-1/penalgamstar);
						}
						else {
							lambdaupdate(i,j) = thetatilde;
						}

                    }
                    lambdasub = lambdaupdate.rows(i,i);
                    mattemp = -2.0*lambdasub*bmi+lambdasub*matrixam*trans(lambdasub);
                    //mattodouble = mattemp(0,0);
                    mattodouble = as_scalar(mattemp);
                    psiupdate(i,i) =  scov(i,i) + mattodouble;
                }
    
                loadingdif = abs(loadinghatnow - lambdaupdate);
                psidif = abs(psiupdate - psihatnow);
                loadingmaxdif = loadingdif.max();
                psimaxdif = psidif.max();
    
                if( loadingmaxdif < tol && psimaxdif < tol) tolproceed = 0L;
            
			    loadinghatnow = lambdaupdate;
                psihatnow = psiupdate;
    
                iteration = iteration + 1;
            }
            
            if( ipenalgam==0 ){
            	loadinghatwarm = loadinghatnow ;
                psihatwarm = psihatnow;
			}
        
            lambdaout.col(count) = vectorise(loadinghatnow);
            psiout.col(count) = diagvec(psihatnow);
            penaltyvalue(0,count) = penalvaluelam;
            penaltyvalue(1,count) = penalvaluegam;
            itercheck(count) = iteration-1;
            iterdiff(0,count) =  loadingmaxdif;
            iterdiff(1,count) =  psimaxdif;
            
            count = count + 1;
		}
        
   }
  
  
  return Rcpp::List::create(
    Rcpp::Named("loading") = lambdaout,
    Rcpp::Named("unique") = psiout,
    Rcpp::Named("lamgam") = penaltyvalue,
    Rcpp::Named("steps") = itercheck,
    Rcpp::Named("diff") = iterdiff
  );
  
}


//[[Rcpp::export]]
List EFAMCPCorr(arma::mat loadingm, arma::mat psim, arma::mat scov, int indicatorn, int factorn, int maxstep, arma::vec lamvec, int lamvecn, arma::vec gamvec, int gamvecn, double tol, arma::mat loadpenweight) {

  int iteration;
  mat modelcov;
  mat matrixm;
  mat matrixam;
  mat bmi;
  mat s1;
  double sumtemp = 0.0;
  double thetatilde = 0.0;
  //double evasign = 0.0;
  //double signcoef = 0.0;
  double signtilde = 0.0;
  mat lambdaupdate = zeros<mat>(indicatorn,factorn);
  mat psiupdate = zeros<mat>(indicatorn,indicatorn);
  mat lambdasub;
  mat mattemp;
  double mattodouble;
  mat loadingdif;
  mat psidif;
  double loadingmaxdif = 1.0;
  double psimaxdif = 1.0;
  int tolproceed = 1;
  double penalvaluelam;
  double penalvaluegam;
  double penalvaluelamgam;
  double penallamstar;
  double penalgamstar;
  double posithetatilde;
  int nonsingular;
  mat lambdahatold;
  mat psihatold;
  mat loadinghatnow;
  mat psihatnow;
  mat loadinghatwarm;
  mat psihatwarm;
  mat lambdaout = zeros<mat>(indicatorn*factorn,lamvecn*gamvecn);
  mat psiout = zeros<mat>(indicatorn,lamvecn*gamvecn);
  mat penaltyvalue = zeros<mat>(2,lamvecn*gamvecn);
  vec itercheck = zeros<vec>(lamvecn*gamvecn);
  mat iterdiff = zeros<mat>(2,lamvecn*gamvecn);
  
  int count = 0;
  for(int ipenallam=0; ipenallam<lamvecn; ipenallam++){
  	    penalvaluelam = lamvec(ipenallam);
    
        if( ipenallam==0 ){
            loadinghatnow = loadingm;
            psihatnow = psim;
        } 
		else {
			loadinghatnow = loadinghatwarm;
            psihatnow = psihatwarm;
		}
        
        for(int ipenalgam=0; ipenalgam<gamvecn; ipenalgam++) {
        	penalvaluegam = gamvec(ipenalgam);
        	penalvaluelamgam = penalvaluelam*penalvaluegam;
        	
        	iteration = 1;
            tolproceed = 1;
  	
  	        while( iteration <= maxstep && tolproceed == 1){
  	    	    lambdahatold = loadinghatnow;
  	    	    psihatold = psihatnow;

                nonsingular = 0L;
                for(int ni=0; ni<indicatorn; ni++) {
                	if( psihatnow(ni,ni)>0.000001 ) {
                		nonsingular = nonsingular + 1;
					}
 				}
 				
 				if( nonsingular != indicatorn){
 					if( ipenalgam !=0 ) {
 						loadinghatnow = loadinghatwarm;
 					    psihatnow = psihatwarm;
					 } else {
					 	loadinghatnow = loadingm;
 				    	psihatnow = psim;
					 }
 					
 					iteration = maxstep + 2;
 					loadingmaxdif = 1.0;
 					psimaxdif = 1.0;
 					break;
				 }
				 
				modelcov = loadinghatnow*trans(loadinghatnow) + psihatnow;
                matrixm = trans(loadinghatnow)*inv(psihatnow)*loadinghatnow + eye( factorn, factorn ) ;
                matrixam = inv(matrixm) + inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov*inv(psihatnow)*loadinghatnow*inv(matrixm);
    
                lambdaupdate = zeros<mat>(indicatorn,factorn);
                psiupdate = zeros<mat>(indicatorn,indicatorn);
    
                for (int i=0; i<indicatorn; i++) {
                    bmi = inv(matrixm)*trans(loadinghatnow)*inv(psihatnow)*scov.cols(i,i);
      
                    for (int j=0; j<factorn; j++) {
                        sumtemp = 0.0;
                        for (int row=0; row<factorn; row++) {
                            if(row!=j) sumtemp = sumtemp + matrixam(row,j)*loadinghatnow(i,row);
                        }
                        thetatilde = ( bmi(j,0) - sumtemp )/matrixam(j,j);
                        penallamstar = loadpenweight(i,j)*psihatnow(i,i)*penalvaluelam/matrixam(j,j);
                        penalgamstar = matrixam(j,j)*penalvaluegam/( loadpenweight(i,j)*psihatnow(i,i) );
                        posithetatilde = std::abs(thetatilde);
                        
                        if( posithetatilde <= penallamstar){
                        	lambdaupdate(i,j) = 0.0;
						}
						else if( posithetatilde <= penalvaluelamgam ){
                        	signtilde = 0.0;
                        	if( thetatilde > 0.0 ) signtilde = 1.0;
                        	if( thetatilde < 0.0 ) signtilde = -1.0;
                        	
                        	lambdaupdate(i,j) = signtilde*( posithetatilde-penallamstar )/(1-1/penalgamstar);
						}
						else {
							lambdaupdate(i,j) = thetatilde;
						}

                    }
                    lambdasub = lambdaupdate.rows(i,i);
                    mattemp = lambdasub*trans(lambdasub);
                    //mattodouble = mattemp(0,0);
                    mattodouble = as_scalar(mattemp);
                    //psiupdate(i,i) =  scov(i,i) + mattodouble;
                    psiupdate(i,i) =  1.0 - mattodouble;
                }
    
                loadingdif = abs(loadinghatnow - lambdaupdate);
                psidif = abs(psiupdate - psihatnow);
                loadingmaxdif = loadingdif.max();
                psimaxdif = psidif.max();
    
                if( loadingmaxdif < tol && psimaxdif < tol) tolproceed = 0L;
            
			    loadinghatnow = lambdaupdate;
                psihatnow = psiupdate;
    
                iteration = iteration + 1;
            }
            
            if( ipenalgam==0 ){
            	loadinghatwarm = loadinghatnow ;
                psihatwarm = psihatnow;
			}
        
            lambdaout.col(count) = vectorise(loadinghatnow);
            psiout.col(count) = diagvec(psihatnow);
            penaltyvalue(0,count) = penalvaluelam;
            penaltyvalue(1,count) = penalvaluegam;
            itercheck(count) = iteration-1;
            iterdiff(0,count) =  loadingmaxdif;
            iterdiff(1,count) =  psimaxdif;
            
            count = count + 1;
		}
        
   }
  
  
  return Rcpp::List::create(
    Rcpp::Named("loading") = lambdaout,
    Rcpp::Named("unique") = psiout,
    Rcpp::Named("lamgam") = penaltyvalue,
    Rcpp::Named("steps") = itercheck,
    Rcpp::Named("diff") = iterdiff
  );
  
}



//[[Rcpp::export]]
List MCPNoInt(arma::mat y, arma::mat x, int n, int p, arma::vec xnorm, arma::mat betainitial, arma::vec lamvec, 
int lamn, arma::vec gamvec, int gamn, int maxstep, double tol, arma::vec penweight) {

mat betachange = betainitial;
mat betawarm = betachange;
mat ytilde;
double betatilde = 0.0;
double sumtemp;
double lamstar = 0.0;
double gamstar = 0.0;
double absbeta = 0.0;
double lamstargam = 0.0;
double signbeta = 0.0;
int iteration;
mat betaold;
mat betadif;
double maxdif;
double lam;
double gam;
mat betaall;
mat cbbetaall;
vec itercheck = zeros(lamn*gamn);
mat lamgam = zeros(lamn*gamn,2);
int count = 0;

for (int lamloop=0; lamloop<lamn; lamloop++) {
    lam = lamvec(lamloop);
    betachange = betawarm;
    for (int gamloop=0; gamloop<gamn; gamloop++) {
        gam = gamvec(gamloop);

        mat takeout = zeros(n,1);
        iteration = 1;
        maxdif = 1.0;
        ytilde = y - x*betachange;

        while( iteration <= maxstep && maxdif > tol){
            betaold = betachange;
            for (int k=0; k<p; k++) {
                sumtemp = 0.0;
                for (int m=0; m<n; m++) {
                    ytilde(m,0) = ytilde(m,0) + x(m,k)*betachange(k,0) - takeout(m,0);
                    sumtemp = sumtemp + ytilde(m,0)*x(m,k);
                }
                betatilde = sumtemp/xnorm(k);
		
		        if( penweight(k)==0.0 ){
		            betachange(k,0) = betatilde;
		        }
		        else {
		        	lamstar = penweight(k)*lam/xnorm(k);
                    gamstar = gam*xnorm(k)/penweight(k);
		            absbeta = std::abs(betatilde);
                    lamstargam = lamstar*gamstar;
		            if( absbeta <= lamstar ){
                        betachange(k,0) = 0.0;
		            } 
		            else if( absbeta <= lamstargam ){
		                signbeta = 0.0;
			            if( betatilde > 0.0 ){
			                signbeta = 1.0;
			            } 
			            else {
			                signbeta = -1.0;
			            }
			            betachange(k,0) = signbeta*( absbeta-lamstar )/( 1.0 - 1.0/gamstar );
		            } 
		            else {
		                betachange(k,0) = betatilde;
		            }
		        }
		
		        for (int m=0; m<n; m++) {
		            takeout(m,0) = x(m,k)*betachange(k,0);
		        }
	        }  
	
	        betadif = abs(betaold - betachange);
	        maxdif = betadif.max();
	        iteration = iteration + 1;
        }

        if( lamloop==0 && gamloop==0 ){
            betaall = betachange;
        } else {
            cbbetaall = join_rows(betaall, betachange);
                betaall = cbbetaall;
        }

        if( gamloop==0 ){
            betawarm = betachange;
        }
        itercheck(count) = iteration;
        lamgam(count,0) = lam;
        lamgam(count,1) = gam;
        count = count + 1;
    }
}

return Rcpp::List::create(
        Rcpp::Named("parall") = betaall,
        Rcpp::Named("beta1beta2") = lamgam,
        Rcpp::Named("steps") = itercheck
);

}
