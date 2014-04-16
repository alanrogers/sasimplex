#ifndef __GSL_SASIMPLEX_H__
#define __GSL_SASIMPLEX_H__

#include <gsl/gsl_types.h>
#include <gsl/gsl_multimin.h>

GSL_VAR const gsl_multimin_fminimizer_type *gsl_multimin_fminimizer_sasimplex;

void        sasimplex_random_seed(gsl_multimin_fminimizer * minimizer,
                                  unsigned seed);
void        sasimplex_set_temp(gsl_multimin_fminimizer * minimizer,
                               double temperature);
int         sasimplex_set_bounds(gsl_multimin_fminimizer * minimizer,
                                 const gsl_vector *lbound,
                                 const gsl_vector *ubound);
int         sasimplex_randomize_state(gsl_multimin_fminimizer * minimizer,
                                      int rotate, gsl_vector * lo,
                                      gsl_vector * hi,
                                      const gsl_vector * step_size);
double      sasimplex_vertical_scale(gsl_multimin_fminimizer * minimizer);
int         sasimplex_n_iterations(gsl_multimin_fminimizer * minimizer,
                                   double *size,
                                   double tol_fval,
                                   double tol_size,
                                   int nItr, double temperature, int verbose);
int         sasimplex_converged(gsl_multimin_fminimizer * minimizer,
                                double tol_fval, double tol_size);
void        sasimplex_print(gsl_multimin_fminimizer * minimizer);

#endif                       /* __GSL_SASIMPLEX_H__ */
