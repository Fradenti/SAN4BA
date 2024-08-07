#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::cube loops_post_dens(arma::mat mu,
                           arma::mat sig,
                           arma::colvec x_grid,
                           arma::colvec est_L,
                           arma::mat S,
                           arma::cube omega ){
  int nsim   = mu.n_rows;
  int n_grid = x_grid.n_elem;
  arma::cube dens(n_grid,S.n_cols,nsim);
  double dens_dummy = 0.0;
  
  for(int i = 0; i < nsim; i++){
    
    for(int k = 0; k < S.n_cols; k++){
      
      for(int q = 0; q < n_grid; q++){
        
        for(int l = 0; l < est_L[i]; l++){
          
          dens_dummy += omega(l,S(i,k)-1,i) * R::dnorm4( x_grid(q), mu(i,l), sig(i,l), 0 ) ;
          
        }
        dens(q,k,i) = dens_dummy;          
        dens_dummy = 0;
      }
    }
  }
  return(dens);                     
}

