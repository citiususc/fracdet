// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>

using namespace Rcpp;


double dnst(double x, double df, double mean, double sd) {
  return dt(NumericVector::create((x - mean) / sd),
            df)[0] / sd;
}

// Calculate the Standard Deviation (sd) of a non-standardized Student's t X
// X = mu + sigma * T (where T is a standardized Student's t)
double sdStudent(double df, double sigma) {
  double st_df = NA_REAL;
  if (df > 2) {
    st_df = sigma * std::sqrt(df / (df - 2));
  } else if (df > 1 && df <=2) {
    st_df = R_PosInf;
  }
  return st_df;
}

// force the computation of a finite sd of a non-standardized Student's t to
// use in the computations of the limits of the integrals

double finiteSdStudent(double df, double sigma) {
  int  kEXPAND = 20.0;
  double st_sd = sdStudent(df, sigma);
  //traits::is_infinite<REALSXP>
  if (!R_FINITE(st_sd) || R_IsNA(st_sd)) {
    // use large value for st_sd
    st_sd = sigma * kEXPAND;
  }
  return st_sd;
}

// dnorm wrapper
double dnormw(double x, double mean, double sd) {
  return dnorm(NumericVector::create(x), mean, sd)[0];
}



// In addition to the standard error codes for invalid arguments the functions can return the following values
// GSL_EMAXITER  the maximum number of subdivisions was exceeded.
// GSL_EROUND cannot reach tolerance because of roundoff error, or roundoff error was detected in the extrapolation table.
// GSL_ESING a non-integrable singularity or other bad integrand behavior was found in the integration interval.
// GSL_EDIVERGE the integral is divergent, or too slowly convergent to be integrated numerically.
bool is_integral_error(int error) {
  int integral_errors[] = {GSL_EMAXITER, GSL_EROUND, GSL_ESING, GSL_EDIVERGE};
  for (int i = 0; i < 4; i++) {
    if (error == integral_errors[i]) {
      return true;
    }
  }
  return false;
}

Rcpp::List calculate_integral(gsl_function F, double lower, double upper,
                              double abs_tol = 0.0, double rel_tol = 1e-7,
                              unsigned int size = 1e5)  {
  // close the GSL error handler
  gsl_set_error_handler_off();
  // define integration workspace
  gsl_integration_workspace *work = gsl_integration_workspace_alloc(size);

  // declare output variables
  double result = 0.0;
  double abserr = 0.0;


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
                            Rcpp::Named("error_code") = res);
}


// the product of two gaussian kernels N(mu1, sd1) * N(mu2, sd2) can be
//  approximated by just one gaussian kernel: N(mu, sd). This function computes
//  mu and sd for the case mu2 = 0
NumericVector approximateGaussianProduct(double mu1, double sd1, double sd2){
  double v1 = std::pow(sd1, 2);
  double v2 = std::pow(sd2, 2);
  return NumericVector::create(
    _["mean"] = mu1 * v2 / (v1 + v2),
    _["sd"] = std::sqrt( v1 * v2 / (v1 + v2))
  );
}

NumericVector getLimits(double mu1, double sd1, double sd2) {
  int kSIGMATIMES = 6;
  NumericVector gaussian_pars = approximateGaussianProduct(mu1, sd1, sd2);
  return NumericVector::create(
    _["lower"] = gaussian_pars["mean"] - kSIGMATIMES * gaussian_pars["sd"],
    _["upper"] = gaussian_pars["mean"] + kSIGMATIMES * gaussian_pars["sd"]);
}


// structures defining the parameters required for computing the integrals
struct pars0 {
  double d;
  double alpha;
  double beta;
  double tau;
};

struct parsI {
  double d;
  double alpha;
  double beta;
  double n;
  double mu;
};

NumericVector getLimitsInt0(const gsl_function& F) {
  pars0 &pars = * (pars0 *) F.params;
  // calculate the sd of the non-standardized students' t and compute proper
  // limits for the integration based on the parameters of the distributions
  double st_sd = finiteSdStudent(2 * pars.alpha, std::sqrt(pars.beta / pars.alpha));
  return getLimits(pars.d, st_sd, pars.tau);
}

NumericVector getLimitsInt1(const gsl_function& F) {
  parsI &pars = * (parsI *) F.params;

  // calculate the sd of the non-standardized students' ts and get proper limits
  // for the integral based on the parameters of the distributions
  double st_sd1 = finiteSdStudent(2 * pars.alpha, std::sqrt(pars.beta / pars.alpha));
  double st_sd2 = finiteSdStudent(pars.n, pars.mu);
  return getLimits(pars.d, st_sd1, st_sd2);
}



double integralPi0(double x, void* p) {
  pars0 &pars = * (pars0 *) p;
  return dnst(pars.d, 2 * pars.alpha, x, std::sqrt(pars.beta / pars.alpha)) *
    dnormw(x, 0.0, pars.tau);
}

double integralTheta0(double x, void* p) {
  pars0 &pars = * (pars0 *) p;
  return x * dnst(pars.d, 2 * pars.alpha, x, std::sqrt(pars.beta / pars.alpha)) *
    dnormw(x, 0.0, pars.tau);
}

double integralPi1(double x, void* p) {
  parsI &pars = * (parsI *) p;
  return dnst(pars.d, 2 * pars.alpha, x, std::sqrt(pars.beta / pars.alpha)) *
    dnst(x, pars.n, 0.0, pars.mu);
}

double integralTheta1(double x, void* p) {
  parsI &pars = * (parsI *) p;
  return x * dnst(pars.d, 2 * pars.alpha, x, std::sqrt(pars.beta / pars.alpha)) *
    dnst(x, pars.n, 0.0, pars.mu);
}

void createPars0(double d, List& priors, pars0& params) {
  params.alpha = priors["alpha"];
  params.beta = priors["beta"];
  params.d = d;
  params.tau = priors["tau"];
}

void createPars1(double d, List& priors, parsI& params) {
  params.alpha = priors["alpha"];
  params.beta = priors["beta"];
  params.d = d;
  params.mu = priors["mu"];
  params.n = priors["n"];
}


int calculateInt0(double d, List& priors, double (&f)(double x, void* p),
                  double abs_tol, double rel_tol, int size, int key,
                  gsl_integration_workspace * workspace, double& result, double& abserr) {
  pars0 params;
  createPars0(d, priors, params);
  gsl_function F;
  F.function = &f;
  F.params = (void *) &params;

  NumericVector limits = getLimitsInt0(F);
  return  gsl_integration_qag(&F, limits["lower"], limits["upper"],
                              abs_tol, rel_tol, size, key,
                              workspace, &result, &abserr);

}

int calculateInt1(double d, List& priors, double (&f)(double x, void* p),
                  double abs_tol, double rel_tol, int size, int key,
                  gsl_integration_workspace * workspace, double& result, double& abserr) {
  parsI params;
  createPars1(d, priors, params);
  gsl_function F;
  F.function = &f;
  F.params = (void *) &params;

  NumericVector limits = getLimitsInt1(F);
  return  gsl_integration_qag(&F, limits["lower"], limits["upper"],
                              abs_tol, rel_tol, size, key,
                              workspace, &result, &abserr);

}

void handleIntegrationError(int error_code) {
  if (error_code) {
    if (is_integral_error(error_code)){
      Rcpp::warning(gsl_strerror(error_code));
    } else {
      Rcpp::stop(gsl_strerror(error_code));
    }
  }
}

// [[Rcpp::export]]
NumericVector bayesRuleCpp(NumericVector input, List priors,int key = 6,
                    double abs_tol = 0.0, double rel_tol = 1e-7,
                    unsigned int size = 1e5)  {
  // close the GSL error handler
  gsl_set_error_handler_off();
  // define integration workspace
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(size);

  double pi0_result, pi0_abserr, pi1_result, pi1_abserr,
  theta0_result, theta0_abserr, theta1_result, theta1_abserr;
  double p = priors["p"];
  double q = 1 - p;
  int error_code;
  NumericVector output(input.size());
  for(int i=0; i < input.size(); i++) {
    error_code = calculateInt0(input[i], priors, integralPi0,
                               abs_tol,rel_tol, size, key,
                               workspace,
                               pi0_result, pi0_abserr);
    handleIntegrationError(error_code);
    error_code = calculateInt0(input[i], priors, integralTheta0,
                               abs_tol,rel_tol, size, key,
                               workspace,
                               theta0_result, theta0_abserr);
    handleIntegrationError(error_code);
    // type 1 integrals
    error_code = calculateInt1(input[i], priors, integralPi1,
                               abs_tol,rel_tol, size, key,
                               workspace,
                               pi1_result, pi1_abserr);
    handleIntegrationError(error_code);
    error_code = calculateInt1(input[i], priors, integralTheta1,
                               abs_tol,rel_tol, size, key,
                               workspace,
                               theta1_result, theta1_abserr);
    handleIntegrationError(error_code);

    output[i] = (p * theta1_result  + q * theta0_result) / (p * pi1_result + q * pi0_result);
  }
  gsl_integration_workspace_free(workspace);

  return output;

}

