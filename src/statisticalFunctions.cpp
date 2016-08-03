# include <Rcpp.h>

/* density function of a non-standardized Student's t X
 X = mu + sigma * T
 */
double dnst(double x, double df, double mean, double sd) {
  return R::dt((x - mean) / sd, df, 0) / sd;
}

// Calculate the Standard Deviation (sd) of a non-standardized Student's t X
// X = mu + sigma * T (where T is a standardized Student's t)
double sdStudent(double df, double sigma) {
  double st_df = NA_REAL;
  if (df > 2) {
    st_df = sigma * std::sqrt(df / (df - 2));
  } else if (df > 1 && df <= 2) {
    st_df = R_PosInf;
  }
  return st_df;
}

