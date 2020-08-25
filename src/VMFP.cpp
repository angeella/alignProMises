#include <RcppArmadillo.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
// [[Rcpp::export]] 

List vMFP(arma::mat X, float k, arma::mat Q, arma::mat ref_ds, bool scaling, bool reflection) {
  int nc = X.n_cols;

  arma::mat U;
  arma::vec s;
  arma::mat V;
//Put transposes to save memory.
  arma::mat trasp = (ref_ds.t() * X + k * Q.t()).t();
  arma::svd(U, s, V, trasp);
  arma::mat Vt = V.t();
  if(!reflection & scaling){
    arma::mat s_new;
    s_new.eye(nc,nc);
    double  d = det(U * Vt);
    if (d > 0){ 
      s_new(nc,nc) =  1;
      }
    if (d < 0) {
      s_new(nc,nc) =  -1;
      }
    arma::mat R = (U * s_new) * Vt;
    double scale = arma::sum(s_new * s);
    arma::mat Xest = X * R.t() * scale;
    List out = List::create(Named("Xest") = Xest, Named("R") = R);
    return out; 
  }
  if(!reflection & !scaling){
    arma::mat s_new;
    s_new.eye(nc,nc);
    double  d = det(U * Vt);
    if (d > 0){ 
      s_new(nc,nc) =  1;
    }
    if (d < 0) {
      s_new(nc,nc) =  -1;
    }
    arma::mat R = (U * s_new) * Vt;
    arma::mat Xest = X * R.t();
    List out = List::create(Named("Xest") = Xest, Named("R") = R);
    return out;  
  }
  if(reflection & scaling){
    arma::mat R =  U * V;
    double scale =  arma::sum(s);
    arma::mat Xest = X * R.t() * scale;
    List out = List::create(Named("Xest") = Xest, Named("R") = R);
    return out;  
  }

  if(reflection & !scaling){
    arma::mat R =  U * V;
    arma::mat Xest = X * R.t();
    List out = List::create(Named("Xest") = Xest, Named("R") = R);
    return out; 
  }

}

/*** R
X <- matrix(rnorm(5*5),5,5)
out<- vMFP(X, k = 1, Q = matrix(rnorm(5*5),5,5), ref_ds = matrix(rnorm(5*5),5,5), scaling = TRUE, reflection = TRUE)
str(out)
*/
