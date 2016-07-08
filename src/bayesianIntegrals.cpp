// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>

bool is_integral_error(int error) {
  int integral_errors[] = {GSL_EMAXITER, GSL_EROUND, GSL_ESING, GSL_EDIVERGE};
  for (int i = 0; i < 4; i++) {
    if (error == integral_errors[i]) {
      return true;
    }
  }
  return false;
}
// prototype of the function to integrate
double integral_I(double x, void *p);


// In addition to the standard error codes for invalid arguments the functions can return the following values
// GSL_EMAXITER  the maximum number of subdivisions was exceeded.
// GSL_EROUND cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table.
// GSL_ESING a non-integrable singularity or other bad integrand behavior was found in the integration interval.
// GSL_EDIVERGE the integral is divergent, or too slowly convergent to be integrated numerically.

// [[Rcpp::export]]
Rcpp::List qawf_integral(double alpha, double lower, double upper,
                         double abs_tol = 0.0, double rel_tol = 1e-7,
                         unsigned int size = 1e5)  {
  // close the GSL error handler
  gsl_set_error_handler_off();
  // define integration workspace
  gsl_integration_workspace *work = gsl_integration_workspace_alloc(size);

  // declare output varaibles
  double result = 0.0;
  double abserr = 0.0;

  gsl_function F;
  F.function = &integral_I;
  F.params = &alpha;


  int res = gsl_integration_qags(&F, lower, upper, abs_tol, rel_tol, size,
                                 work, &result, &abserr);
  gsl_integration_workspace_free(work);

  if (res) {
    if (is_integral_error(res)){
      Rcpp::warning(gsl_strerror(res));
    } else {
      Rcpp::stop(gsl_strerror(res));
    }
  }

  return Rcpp::List::create(Rcpp::Named("value") = result,
                            Rcpp::Named("abs_error") = abserr,
                            Rcpp::Named("error_code") = res);            // return vector
}

double integral_I(double x, void* p){
  double alpha = * (double *) p;
  return exp(-alpha * x);
}


// /*** R
// alpha = 2
// lower = 0
// upper = 1
// # demonstrate equivalence of R and C results
// C.result = qawf_integral(alpha, lower = lower, upper = upper, rel_tol = 1e-10)$value
// R.result = integrate(function(x) exp(-alpha * x), lower = lower, upper = upper)
//
// print(C.result)
// print(R.result)
// */
