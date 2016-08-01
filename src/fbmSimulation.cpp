#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

NumericVector computeFgnCov(int n, double H){
  NumericVector fgn_cov(n);
  for (int i = 0; i < n; i++){
    fgn_cov[i] =  (std::pow(std::abs(i - 1), 2 * H) - 2 * std::pow(std::abs(i), 2 * H) +
      std::pow(std::abs(i + 1), 2 * H)) / 2;
  }
  return fgn_cov;
}

// Functor to sum two vectors a and b after multiplying them by some factors
// factor_a * a + factor_b * b
class MultiplyVectorsAndSum {
private:
  double mfactor_a;
  double mfactor_b;

public:
  MultiplyVectorsAndSum(double factor_a, double factor_b) : mfactor_a(factor_a), mfactor_b(factor_b) { }

  double operator() (const double& a, const double& b) const {
    return this->mfactor_a * a + this->mfactor_b * b;
  }
};


// [[Rcpp::export]]
NumericVector simulateFbmCpp(int n, double H) {
  // Covariances of fGn:
  NumericVector v1 = computeFgnCov(n, H);
 // Initialization of algorithm:
  NumericVector y = rnorm(n);
  NumericVector fgn(n, 0.0);
  NumericVector v2(n + 1);
  std::copy(v1.begin() + 1, v1.end(), v2.begin() + 1);
  double k =  -v2[1];
  double sigma = std::sqrt(v1[0]);
  // Levinson's algorithm:
  for(int j = 1; j < n; j++) {
    sigma = sigma * sqrt(1 - k * k);
    NumericVector v = v1.size() - j;
    std::transform(v1.begin() + j - 1, v1.end() - 1,
                   v2.begin() + j,
                   v.begin(),
                   MultiplyVectorsAndSum(1,k));
    std::transform(v2.begin() + j, v2.end() - 1,
                   v1.begin() + j - 1,
                   v2.begin() + j,
                   MultiplyVectorsAndSum(1,k)
    );
    std::copy(v.begin(), v.end(), v1.begin() + j);
    std::transform(fgn.begin() + j, fgn.end(),
                   v1.begin() + j,
                   fgn.begin() + j,
                   MultiplyVectorsAndSum(1, y[j] / sigma));
    // update k
    k =  - v2[j + 1] / (sigma * sigma);
  }
  NumericVector fbm = cumsum(fgn);
  return fbm;
}
