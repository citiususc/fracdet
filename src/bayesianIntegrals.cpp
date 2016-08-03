// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>
#include "bayesianIntegrals.h"
#include "statisticalFunctions.h"
using namespace Rcpp;


/* The product of two gaussian kernels N(mu1, sd1) * N(mu2, sd2) can be
 * approximated by just one gaussian kernel: N(mu, sd). Since all the integrands
 * involved in the computation of the bayes rule are similar to gaussian kernels,
 * we may compute a finite interval for all the integrals by using these approach
 */

/* Force the computation of a finite sd of a non-standardized Student's t (the
 * variance is not well-defined if df < 2) to use in the computations of the
 * finite limits of the integrals
 */
double finiteSdStudent(double df, double sigma) {
  int  kEXPAND = 20;
  double st_sd = sdStudent(df, sigma);
  // check if the variance is defined
  if (!R_FINITE(st_sd) || R_IsNA(st_sd)) {
    // use large value for st_sd
    st_sd = sigma * kEXPAND;
  }
  return st_sd;
}

/* The product of two gaussian kernels N(mu1, sd1) * N(mu2, sd2) can be
 * approximated by just one gaussian kernel: N(mu, sd). This function computes
 * mu and sd for the case mu2 = 0 involved in the integrals.
 */
NumericVector approxGaussianProduct(double mu1, double sd1, double sd2){
  double v1 = std::pow(sd1, 2);
  double v2 = std::pow(sd2, 2);
  return NumericVector::create(
    _["mean"] = mu1 * v2 / (v1 + v2),
    _["sd"] = std::sqrt( v1 * v2 / (v1 + v2))
  );
}

/* Compute a finite interval for the integrals consisting on the product of
 * two gaussian kernels:
 */
NumericVector getIntegralLimits(double mu1, double sd1, double sd2) {
  int kSIGMATIMES = 6;
  NumericVector gaussian_pars = approxGaussianProduct(mu1, sd1, sd2);
  return NumericVector::create(
    _["lower"] = gaussian_pars["mean"] - kSIGMATIMES * gaussian_pars["sd"],
    _["upper"] = gaussian_pars["mean"] + kSIGMATIMES * gaussian_pars["sd"]);
}


/*
 * Get finite limits for the integrals of type 0: those involving one Student's t
 * and a gaussian density functions (noted with subscript 0 in the paper).
 */
NumericVector getIntegralLimitsT0(const gsl_function& F) {
  parsT0 &pars = * (parsT0 *) F.params;
  double st_sd = finiteSdStudent(2 * pars.alpha, std::sqrt(pars.beta / pars.alpha));
  return getIntegralLimits(pars.d, st_sd, pars.tau);
}

/*
 * Get finite limits for the integrals of type 1: those involving two Student's t
 * density functions (noted with subscript 1 in the paper).
 */
NumericVector getIntegralLimitsT1(const gsl_function& F) {
  parsT1 &pars = * (parsT1 *) F.params;
  double st_sd1 = finiteSdStudent(2 * pars.alpha, std::sqrt(pars.beta / pars.alpha));
  double st_sd2 = finiteSdStudent(pars.n, pars.mu);
  return getIntegralLimits(pars.d, st_sd1, st_sd2);
}

/*
 * Definitions of the integrals involved in the computation of the bayes
 * decision rule
 */
double integrandPi0(double x, void* p) {
  parsT0 &pars = * (parsT0 *) p;
  return dnst(pars.d, 2 * pars.alpha, x, std::sqrt(pars.beta / pars.alpha)) *
    R::dnorm(x, 0.0, pars.tau, 0);
}

double integrandLambda0(double x, void* p) {
  parsT0 &pars = * (parsT0 *) p;
  return x * dnst(pars.d, 2 * pars.alpha, x, std::sqrt(pars.beta / pars.alpha)) *
    R::dnorm(x, 0.0, pars.tau, 0);
}

double integrandPi1(double x, void* p) {
  parsT1 &pars = * (parsT1 *) p;
  return dnst(pars.d, 2 * pars.alpha, x, std::sqrt(pars.beta / pars.alpha)) *
    dnst(x, pars.n, 0.0, pars.mu);
}

double integrandLambda1(double x, void* p) {
  parsT1 &pars = * (parsT1 *) p;
  return x * dnst(pars.d, 2 * pars.alpha, x, std::sqrt(pars.beta / pars.alpha)) *
    dnst(x, pars.n, 0.0, pars.mu);
}



/*
Check if the status code returned the gsl integration function corresponds
to an error related with the integration process:
* GSL_EMAXITER  the maximum number of subdivisions was exceeded.
* GSL_EROUND cannot reach tolerance because of roundoff error, or roundoff error
was detected in the extrapolation table.
* GSL_ESING a non-integrable singularity or other bad integrand behavior was
found in the integration interval.
* GSL_EDIVERGE the integral is divergent, or too slowly convergent to be
integrated numerically.
*/
bool isIntegrationError(int error) {
  int integral_errors[] = {GSL_EMAXITER, GSL_EROUND, GSL_ESING, GSL_EDIVERGE};
  for (int i = 0; i < 4; i++) {
    if (error == integral_errors[i]) {
      return true;
    }
  }
  return false;
}

/*
 * Handle possible integration errors
 */
void handleIntegrationError(int error_code, gsl_integration_workspace * workspace) {
  if (error_code) {
    if (isIntegrationError(error_code)){
      Rcpp::warning(gsl_strerror(error_code));
    } else {
      gsl_integration_workspace_free(workspace);
      Rcpp::stop(gsl_strerror(error_code));
    }
  }
}

/*
 * Calculate integrals of type 0
 */
int calculateIntT0(double d, List& priors, double (&f)(double x, void* p),
                  double abs_tol, double rel_tol, int size, int key,
                  gsl_integration_workspace * workspace, double& result, double& abserr) {
  parsT0 params;
  params.alpha = priors["alpha"];
  params.beta = priors["beta"];
  params.d = d;
  params.tau = priors["tau"];

  gsl_function F;
  F.function = &f;
  F.params = (void *) &params;

  NumericVector limits = getIntegralLimitsT0(F);
  int error_code =  gsl_integration_qag(&F, limits["lower"], limits["upper"],
                                        abs_tol, rel_tol, size, key,
                                        workspace, &result, &abserr);
  handleIntegrationError(error_code, workspace);
  return error_code;

}

/*
 * Calculate integrals of type 1
 */
int calculateIntT1(double d, List& priors, double (&f)(double x, void* p),
                  double abs_tol, double rel_tol, int size, int key,
                  gsl_integration_workspace * workspace, double& result, double& abserr) {
  parsT1 params;
  params.alpha = priors["alpha"];
  params.beta = priors["beta"];
  params.d = d;
  params.mu = priors["mu"];
  params.n = priors["n"];

  gsl_function F;
  F.function = &f;
  F.params = (void *) &params;

  NumericVector limits = getIntegralLimitsT1(F);
  int error_code =   gsl_integration_qag(&F, limits["lower"], limits["upper"],
                                     abs_tol, rel_tol, size, key,
                                     workspace, &result, &abserr);
  handleIntegrationError(error_code, workspace);
  return error_code;

}


// [[Rcpp::export]]
NumericVector bayesRuleCpp(NumericVector input, List priors, int key = 6,
                    double abs_tol = 0.0, double rel_tol = 1e-7,
                    unsigned int size = 1e5)  {
  // close the GSL error handler
  gsl_set_error_handler_off();
  // define integration workspace
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(size);

  // functions for storing the results of each of the integrals
  double pi0, pi0_abserr, pi1, pi1_abserr,
  lambda0, lambda0_abserr, lambda1, lambda1_abserr;
  int error_code;

  double p = priors["p"];
  double q = 1 - p;

  NumericVector output(input.size());
  for(int i=0; i < input.size(); i++) {
    error_code = calculateIntT0(input[i], priors, integrandPi0,
                               abs_tol,rel_tol, size, key,
                               workspace,
                               pi0, pi0_abserr);

    error_code = calculateIntT0(input[i], priors, integrandLambda0,
                               abs_tol,rel_tol, size, key,
                               workspace,
                               lambda0, lambda0_abserr);
    // type 1 integrals
    error_code = calculateIntT1(input[i], priors, integrandPi1,
                               abs_tol,rel_tol, size, key,
                               workspace,
                               pi1, pi1_abserr);

    error_code = calculateIntT1(input[i], priors, integrandLambda1,
                               abs_tol,rel_tol, size, key,
                               workspace,
                               lambda1, lambda1_abserr);

    output[i] = (p * lambda1  + q * lambda0) / (p * pi1 + q * pi0);
  }
  gsl_integration_workspace_free(workspace);

  return output;

}
