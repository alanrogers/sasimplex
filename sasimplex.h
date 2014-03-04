#ifndef __GSL_SASIMPLEX_H__
#define __GSL_SASIMPLEX_H__

#include <gsl/gsl_types.h>
#include <gsl/gsl_multimin.h>

GSL_VAR const gsl_multimin_fminimizer_type *
    gsl_multimin_fminimizer_sasimplex;
GSL_VAR const gsl_multimin_fminimizer_type *
    gsl_multimin_fminimizer_sasimplexrand;

void sasimplex_seed_rng(gsl_multimin_fminimizer *minimizer, unsigned long seed);
void sasimplex_set_temp(gsl_multimin_fminimizer *minimizer, double temperature);

#endif /* __GSL_SASIMPLEX_H__ */

