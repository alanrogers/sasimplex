#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_multimin.h>
#include "sasimplex.h"
#include "annealsched.h"

double      my_f(const gsl_vector * v, void *params);

/*
 * Multiple minima wherever x and y are integers. Global minimum
 * at (x,y)=(par[0],par[1]).
 */
double my_f(const gsl_vector * v, void *params) {
    double      x, y;
    double      fx, fy;         /* fractional parts of x and y */
    double      rval;
    double     *par = (double *) params;

    x = gsl_vector_get(v, 0) - par[0];
    y = gsl_vector_get(v, 1) - par[1];
    fx = x - floor(x + 0.5);
    fy = y - floor(y + 0.5);

    rval = 1.0 + fabs(fx) + fabs(fy);
    rval *= 1.0 + fabs(x) + fabs(y);

    return rval;
}

#define STATEDIM 2

int main(void) {
    double      par[STATEDIM] = { 3.0, 2.0 };  /* (x,y) at minumum */
    const double tol = 1e-4;
    const int   maxItr = 1000;
    double      initStepSize = 0.05;
    const int   randomInitialPoint = 1; /* start at a random point? */
    const int   rotate = 1;             /* random rotation of init simplex?*/
    unsigned long seed = time(NULL);    /* for random numbers */

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_sasimplex;
    unsigned    i;
    int         status, itr = 0;
    double      size;

    /* Set initial step sizes to initStepSize */
    gsl_vector *ss = gsl_vector_alloc(STATEDIM);
    gsl_vector_set_all(ss, initStepSize);

    /* Set up annealing schedule */
    AnnealSched *sched = AnnealSched_alloc(5,   /* number of temperatures */
                                           100, /* iterations per temperature */
                                           4.0, /* initial temperature */
                                           0.5  /* ratio of adjacent temps */
        );

    /* Starting point */
    gsl_vector *x = gsl_vector_alloc(STATEDIM);
    double      initVal[STATEDIM] = {5.5, 7.5};
    for(i = 0; i < STATEDIM; ++i)
        gsl_vector_set(x, i, initVal[i]);

    /* Initialize method and iterate */
    gsl_multimin_function minex_func;

    minex_func.n = STATEDIM;
    minex_func.f = my_f;
    minex_func.params = par;

    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, STATEDIM);
    if(s == NULL) {
        fprintf(stderr, "%s:%d: bad allocation\n", __FILE__, __LINE__);
        exit(1);
    }

    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
    sasimplex_random_seed(s, seed);
    if(randomInitialPoint) {
        gsl_vector *loInit = gsl_vector_alloc(STATEDIM);
        gsl_vector *hiInit = gsl_vector_alloc(STATEDIM);
        gsl_vector_set_all(loInit, -4.5);
        gsl_vector_set_all(hiInit, 4.5);

        sasimplex_randomize_state(s, rotate, loInit, hiInit, ss);
    }

    printf("Using minimizer %s.\n", gsl_multimin_fminimizer_name(s));
    printf("%5s %10s %10s %7s %8s %8s %8s\n", "itr",
           "x", "y", "f", "size", "AbsErr", "temp");
    do {
        double temperature = AnnealSched_next(sched);
        sasimplex_set_temp(s, temperature);
        itr++;
        status = gsl_multimin_fminimizer_iterate(s);
        if(status) {
            printf("%s:%d:%s: rtn val %d from %s\n",
                   __FILE__, __LINE__, __func__,
                   status, "gsl_multimin_fminimizer_iterate");
            break;
        }

        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, tol);

        /* absErr is summed absolute error */
        double      errx = gsl_vector_get(s->x, 0) - par[0];
        double      erry = gsl_vector_get(s->x, 1) - par[1];
        double      absErr = fabs(errx) + fabs(erry);

        printf("%5d %10.3e %10.3e %7.3f %8.3f %8.4f %8.4f\n",
               itr,
               gsl_vector_get(s->x, 0),
               gsl_vector_get(s->x, 1), s->fval, size, absErr, temperature);
    } while(status == GSL_CONTINUE && itr < maxItr);
    switch (status) {
    case GSL_SUCCESS:
        printf("converged to minimum\n");
        break;
    default:
        printf("no convergence: status=%d itr=%d\n", status, itr);
    }
    printf("True minimum: (%lf, %lf)\n", par[0], par[1]);

    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    return status;
}
