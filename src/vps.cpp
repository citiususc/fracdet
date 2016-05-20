#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <R.h>

using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends("RcppArmadillo")]]


// implements the R operation convolve(x, y, "o")
arma::vec arma_convolve(arma::vec x, arma::vec y) {
  // extend x
  int x_len = x.n_elem;
  arma::vec v_zeros(y.n_elem - 1);
  v_zeros.fill(0.0);
  x = join_cols(v_zeros, x);
  // extend y
  v_zeros.set_size(x_len - 1);
  v_zeros.fill(0.0);
  y = join_cols(y, v_zeros);
  // compute the convolution using the FFT for efficiency
  return arma::real(arma::ifft( arma::fft(x) % arma::conj(arma::fft(y))) );
}

// Downsamples the symmetric signal x[n] whose first sample starts in n=first_index
// The function ensures that the downsampled vector is also symmetric by
// using the data from the first half of data (from n=first_index to n=0)
arma::vec downsampleSymmetricSignal(arma::vec x, int first_index) {
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
NumericVector fastVps(NumericVector gr, NumericVector hr, double H, double sigma2, double Ts, int nlevels) {
  arma::vec g(gr.begin(), gr.size(), false);
  arma::vec h(hr.begin(), hr.size(), false);

  arma::vec g_acf = arma_convolve(g, g);
  // we will use only the reversed version of h_acf vector
  arma::vec rev_h_acf = arma_convolve(h, h);
  std::reverse(rev_h_acf.begin(), rev_h_acf.end());

  NumericVector vps(nlevels);
  int max_D_index = (pow(2, nlevels) - 1) * (g.n_elem - 1) + 1;
  arma::vec D(2 * max_D_index + 1);
  for (int i = 0; i < D.n_elem; i++) {
    D[i] = pow(std::abs(-max_D_index + i), 2 * H);
  }
  // Iterate each resolution level
  int central_position;
  double vps_factor = -sigma2 * pow(Ts, 2 * H) /  2.0;
  for (int i = 0; i < nlevels; i ++) {
    central_position = (int)((D.n_elem - 1) / 2);
    vps[i] = vps_factor * arma::dot(D.subvec(central_position - g.n_elem + 1, central_position + g.n_elem -1),
                                    g_acf);
    // update D
    D = arma_convolve(D, rev_h_acf);
    int first_index = (int) -((D.n_elem - 1) / 2);
    D = downsampleSymmetricSignal(D, first_index);
  }
  return  vps;
}
/*** R
getWaveletFilters <- function(filter.number = 1, family = "DaubExPhase"){
  h <- filter.select(filter.number = filter.number,
                     family=family)$H
  len <- length(h)
  # Definition of g[n]:
  # g[n] = (-1)^(1-n)h[1-n]
  index <- 1 + - (len - 1):0
  g <- rev(h) * (-1)^(1 - index)
  list(h = h, g = g)
}


filter.number = 8
family = "DaubLeAsymm"
H = 0.3
sigma2 = 1
Ts =  1
nlevels = 7
plot(a<-expectedVps(H,sigma2,Ts,filter.number = filter.number,family = family,nlevels))
filters = getWaveletFilters(filter.number,family)
points(b <-fastVps(filters$g, filters$h,H, sigma2,Ts ,nlevels), pch = 3, col = 3)

print(a)
print(b)

library(microbenchmark)
microbenchmark(R={a<-expectedVps(H,sigma2,Ts,filter.number = filter.number,family = family,nlevels)},
               C={filters = getWaveletFilters(filter.number,family)
               b <-fastVps(filters$g, filters$h,H, sigma2,Ts ,nlevels)})

*/


