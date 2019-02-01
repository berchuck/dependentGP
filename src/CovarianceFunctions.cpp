#include <RcppArmadillo.h>

//Function to calculate temporal correlation structure-------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat H(double Phi, arma::mat const& Dist, int KernInd, int T) {

  //Create ouput object
  arma::mat out(T, T);

  //exponential
  if (KernInd == 0) out = arma::exp(-Phi * Dist);

  //squared exponential
  if (KernInd == 1) {
    out = arma::exp(-(Dist / (Phi * Phi)));
    arma::mat EyeT(T, T, arma::fill::eye);
    out += (0.000001 * EyeT);
  }

  //ar(1) continuous
  if (KernInd == 2) {
    out = arma::eye(T, T);
    for (int j = 0; j < T; j++) {
      for (int k = 0; k < j; k++) {
        out(j, k) = std::pow(Phi, Dist(j, k));
      }
    }
    out = arma::symmatl(out);
  }
  return out;
}
