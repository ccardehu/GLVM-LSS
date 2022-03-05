// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp ;

/*
******************************************************************************
Offset and Array classes based on code by Romain Francois copied from
http://comments.gmane.org/gmane.comp.lang.r.rcpp/5932 on 2014-01-07.
******************************************************************************
*/

class Offset{
private:
 IntegerVector dim ;

public:
 Offset( IntegerVector dim ) : dim(dim) {}
 int operator()( IntegerVector ind ){
  int ret = ind[0] ;
  int offset = 1 ;
   for(int d=1; d < dim.size(); d++) {
    offset = offset * dim[d-1] ; 
    ret = ret + ind[d] * offset ;
   }
  return ret ;
 } ;
 IntegerVector getDims() const {
  return(dim) ;
 };
};

class Array : public NumericVector {
private:
 // NumericVector value;
 Offset dims ;

public:
 //Rcpp:as
 Array(SEXP x) : NumericVector(x), dims( (IntegerVector)((RObject)x).attr("dim") ) {}
 Array( NumericVector x,  Offset d ): NumericVector(x), dims(d) {}
 Array( Dimension d ): NumericVector( d ), dims( d ) {}
 IntegerVector getDims() const { return(dims.getDims()); };
 NumericVector getValue()  const { return(*((NumericVector*)(this))); };
 inline double& operator()( IntegerVector ind) {
  int vecind = dims(ind);
  NumericVector value = this->getValue();  
  return value(vecind);
 } ;
 // change dims without changing order of elements (!= aperm)
 void resize(IntegerVector newdim) {
 int n = std::accumulate((this->getDims()).begin(), (this->getDims()).end(), 1,  std::multiplies<int>());
 int nnew = std::accumulate(newdim.begin(), newdim.end(), 1, std::multiplies<int>());
 if(n != nnew)  stop("old and new old dimensions don't match.");
 this->dims = Offset(newdim);
 } ;
} ;

namespace Rcpp {
 // wrap(): converter from Array to an R array
 template <> SEXP wrap(const Array& A) {
 IntegerVector dims = A.getDims();
 //Dimension dims = A.getDims();
 Vector<REALSXP> x = A;
 x.attr( "dim" ) = wrap(dims);
 return x; 
 }
}
