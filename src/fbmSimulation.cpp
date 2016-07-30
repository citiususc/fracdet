#include <Rcpp.h>
//#include <cstdlib>
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


class MultiplyVectorsAndSum {
private:
  double factor_a;
  double factor_b;
  
public:
  MultiplyVectorsAndSum(double x, double y) : factor_a(x), factor_b(y) { }
  
  double operator () (const double& a, const double& b) const {
    return this->factor_a * a + this->factor_b * b;
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
  double aa = std::sqrt(v1[0]);
  double bb;
  // Levinson's algorithm:
  for(int j = 1; j < n; j++) {
    aa = aa * sqrt(1 - k * k);
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
    bb = y[j] / aa;
    std::transform(fgn.begin() + j, fgn.end(),
                   v1.begin() + j,
                   fgn.begin() + j,
                   MultiplyVectorsAndSum(1,bb));
    k =  - v2[j + 1] / (aa * aa);
  }
  NumericVector fbm = cumsum(fgn);
  return fbm;
}
