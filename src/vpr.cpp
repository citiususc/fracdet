#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <R.h>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends("RcppArmadillo")]]


// Implements the R operation convolve(x, y, type="o")
// Since the usual definition of convolution is implemented in R as
// convolve(x, rev(y), type="o"), and this is equivalent to arma::conv(x,y), we
// implement the required operation as arma::conv(x,rev(y))
arma::vec r_convolve(const arma::vec &x,const arma::vec &y) {
  arma::vec y_reverse(y.n_elem);
  for (int i = 0; i < y.n_elem; i++) {
    y_reverse[i] = y[y.n_elem - i - 1];
  }
  return arma::conv(x, y_reverse);
}

// Downsamples the symmetric signal x[n] whose first sample starts in n=first_index
// The function ensures that the downsampled vector is also symmetric by
// using the data from the first half of data (from n=first_index to n=0)
arma::vec downsampleSymmetricSignal(const arma::vec &x, int first_index) {
  // start at 0 if first_index is even
  int start_at = std::abs(first_index % 2);
  if (start_at == 1) {
    first_index++;
  }
  // final vector will have -first_index + 1 values
  arma::vec out(-first_index + 1);
  // create a symmetric signal using the first half values from x
  // x is iterated in steps of two to implement downsampling
  int n_half_x = (int)(-first_index / 2);
  int i;
  for (i = 0; i < n_half_x; i++){
    out[i] = x[start_at + i * 2];
    out[-first_index - i] = x[start_at + i * 2];
  }
  // sample at the middle
  out[i] = x[start_at + i * 2];
  return out;
}

// [[Rcpp::export]]
NumericVector theoreticalVpr(NumericVector& gr, NumericVector& hr,
                      double H, double sigma2, int nlevels) {
  arma::vec g(gr.begin(), gr.size(), false);
  arma::vec h(hr.begin(), hr.size(), false);

  arma::vec g_acf = r_convolve(g, g);
  // we will use only the reversed version of h_acf vector
  arma::vec rev_h_acf = r_convolve(h, h);
  std::reverse(rev_h_acf.begin(), rev_h_acf.end());
  // vector for storing the variance per resolution level (vpr)
  NumericVector vpr(nlevels);
  int max_D_index = (pow(2, nlevels) - 1) * (g.n_elem - 1) + 1;
  arma::vec D(2 * max_D_index + 1);
  for (int i = 0; i < D.n_elem; i++) {
    D[i] = pow(std::abs(-max_D_index + i), 2 * H);
  }
  // Iterate each resolution level
  int central_position;
  for (int i = 0; i < nlevels; i ++) {
    central_position = (int)((D.n_elem - 1) / 2);
    vpr[i] = -sigma2 /  2.0 * arma::dot(D.subvec(central_position - g.n_elem + 1,
                                                 central_position + g.n_elem -1),
                                        g_acf);
    // update D
    D = r_convolve(D, rev_h_acf);
    int first_index = (int) -((D.n_elem - 1) / 2);
    D = downsampleSymmetricSignal(D, first_index);
  }
  return  rev(vpr);
}

