#ifndef __GSL_SASIMPLEX_H__
#define __GSL_SASIMPLEX_H__

#include <gsl/gsl_types.h>
#include <gsl/gsl_multimin.h>

GSL_VAR const gsl_multimin_fminimizer_type *gsl_multimin_fminimizer_sasimplex;

void        sasimplex_random_seed(gsl_multimin_fminimizer * minimizer,
                                  unsigned seed);
void        sasimplex_set_temp(gsl_multimin_fminimizer * minimizer,
                               double temperature);
int         sasimplex_randomize_state(gsl_multimin_fminimizer *minimizer,
                          int rotate, gsl_vector * lo,
                          gsl_vector * hi,
                          const gsl_vector * step_size);
#endif                       /* __GSL_SASIMPLEX_H__ */
