#ifndef fYZClass_H
#define fYZClass_H

class fYZ {
  // Private class members
  arma::mat pD ;
  double ll ;

public:

  // Constructors

  fYZ(){}

  fYZ(Rcpp::NumericMatrix& Yr,
      Rcpp::CharacterVector& famr,
      ghQ& Q,
      arma::cube& b) {

    int p = Yr.ncol(), n = Yr.nrow() ;
    arma::mat Yc(Yr.begin(),n,p,false);
    arma::mat Qp = Q.get_grid() ;
    arma::vec Qw = Q.get_weig() ;
    int tqp = Qp.n_rows ;

    pD.set_size(n,tqp);
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

    arma::vec Cm(n), llv(n);
    for(int jj = 0; jj < n; jj++){
      Cm(jj) = arma::max(pD.row(jj));
      llv(jj) = arma::as_scalar(Cm(jj)) + arma::as_scalar(arma::log(arma::exp(pD.row(jj) - arma::as_scalar(Cm(jj))) * Qw));
    }

    ll = arma::accu(llv);

    pD.each_col() -= Cm ;
    pD = arma::exp(pD) ;
    pD.each_col() /= arma::exp(llv - Cm);

  }

  fYZ(fYZ& fYZobj){
    pD = fYZobj.get_pD();
    ll = fYZobj.get_ll();
  }

  // Member functions
  void update(Rcpp::NumericMatrix& Yr,
              Rcpp::CharacterVector& famr,
              ghQ& Q,
              arma::cube& b);
  arma::mat get_pD();
  double get_ll();
};

#endif
