#ifndef CITIUS_FRACDET_H
#define CITIUS_FRACDET_H

/*
 * parameter structures required for the defining gsl_functions involved in
 * the integration of the bayes decision rule.
 * We distinguish two types of integrands:
 * -type 0: those that involve a Student's t and a gaussian function (noted with
 *  the subscript 0 in the article).
 * -type 1: those that involve two Student's ts distributions (noted with the
 *  subscript 1 in the article).
 *
 *  In the structures: d represents the input value to the bayes rule
 */
struct parsT0 {
  double d;
  double alpha;
  double beta;
  double tau;
};

struct parsT1 {
  double d;
  double alpha;
  double beta;
  double n;
  double mu;
};

#endif