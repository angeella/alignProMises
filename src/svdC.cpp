#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 

List svdC(arma::mat X) {
// Output matrices
arma::mat U;
arma::vec S;
arma::mat V;

// Perform SVD
arma::svd(U, S, V, X);
List out = List::create(Named("U")=U, Named("S") = S, Named("V") = V);
return out;
}

/*** R
#X <- matrix(rnorm(5*5),5,5)
#out<- svdC(X)
#str(out)
*/