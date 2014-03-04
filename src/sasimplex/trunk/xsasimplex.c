#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include "sasimplex.h"

/*
 * This code was copied from the gsl manual section on
 * multidimensional minimization.
 *
 * Alan Rogers
 */

double my_f (const gsl_vector *v, void *params);

/*
 * Multiple minima wherever x and y are integers. Global minimum
 * at (x,y)=(0,0)
 */
double
my_f (const gsl_vector *v, void *params)
{
    double x, y;
    double fx, fy; /* fractional parts of x and y */
    double rval;
  
    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);
    fx = x - floor(x + 0.5);
    fy = y - floor(y + 0.5);

    rval = 1.0 + fabs(fx) + fabs(fy);
    rval *= 1.0 + fabs(x) + fabs(y);
    
    return rval; 
}

int 
main(void){
  double par[5] = {1.0, 2.0, 10.0, 20.0, 30.0};
  const size_t stateDim = 2; /**< dimension of state space */
  const double tol = 1e-4;
  const int maxItr = 1000;
  double temperature = 10.0;
  double initStepSize=1.0;

  /* initial coordinates */
  double initVal[stateDim];
  initVal[0] = 5.5;
  initVal[1] = 7.5;

#if 0
  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_sasimplex;
#else
  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_sasimplexrand;
#endif
  gsl_vector *x = gsl_vector_alloc (stateDim);
  gsl_vector *ss = gsl_vector_alloc (stateDim);

  unsigned i;
  int status, itr = 0;
  double size; 

  /* Starting point */
  for(i=0; i<stateDim; ++i)
      gsl_vector_set (x, i, initVal[i]);

  /* Set initial step sizes to initStepSize */
  gsl_vector_set_all (ss, initStepSize);

  /* Initialize method and iterate */
  gsl_multimin_function minex_func;
  minex_func.n = stateDim;
  minex_func.f = my_f;
  minex_func.params = par;

  gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, stateDim);
  if(s==NULL) {
      fprintf(stderr, "%s:%d: bad allocation\n", __FILE__,__LINE__);
      exit(1);
  }
  sasimplex_seed_rng(s, 0); /* 0 ==> use clock */
  sasimplex_set_temp(s, temperature);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  printf("Using minimizer %s.\n", gsl_multimin_fminimizer_name(s));
  printf("%s:%d: temperature=%lf\n", __FILE__,__LINE__,temperature);
  printf ("%5s %10s %10s %7s %8s\n", "itr",
          "x", "y", "f", "size");
  do{
      itr++;
      status = gsl_multimin_fminimizer_iterate(s);
      if(status) {
          printf("%s:%d:%s: rtn val %d from gsl_multimin_fminimizer_iterate\n",
                 __FILE__,__LINE__,__func__,status);
          break;
      }

      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, tol);

      printf ("%5d %10.3e %10.3e %7.3f %8.3f\n", 
              itr,
              gsl_vector_get(s->x, 0), 
              gsl_vector_get(s->x, 1), 
              s->fval, size);
      if(temperature > 0.0 && itr % 10 == 0) {
          temperature -= 1.0;
          sasimplex_set_temp(s, temperature);
          printf("%s:%d: temperature=%lf\n", __FILE__,__LINE__,temperature);
      }
  }while (status == GSL_CONTINUE && itr < maxItr);
  switch(status) {
  case GSL_SUCCESS:
      printf("converged to minimum\n");
      break;
  default:
      printf("no convergence: status=%d itr=%d\n", status, itr);
  }
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return status;
}

