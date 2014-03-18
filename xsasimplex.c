#include "sasimplex.h"
#include "annealsched.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_multimin.h>

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
    double      initStepSize = 0.05;
    const int   rotate = 1;             /* random rotation of init simplex?*/
	const int   verbose = 1;	
    unsigned long seed = time(NULL);    /* for random numbers */
	unsigned    nTries = 10;            /* number of random starts */
    int         nT  = 10;               /* number of temperatures */
    int         nPerT = 100;            /* iterations per temperature */
    double      initT = 3.0;            /* initial temperature */
    double      decay = 0.7;

    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_sasimplex;
    unsigned    i, try;
    int         status;
    double      size, absErr, temperature;

    /* Set initial step sizes to initStepSize */
    gsl_vector *ss = gsl_vector_alloc(STATEDIM);
    gsl_vector_set_all(ss, initStepSize);

    /* Set up annealing schedule */
    AnnealSched *sched = AnnealSched_alloc(nT, initT, decay);

    /* Initial state vector */
    double      initVal[STATEDIM] = {5.5, 7.5};
    gsl_vector *x = gsl_vector_alloc(STATEDIM);
    for(i = 0; i < STATEDIM; ++i)
        gsl_vector_set(x, i, initVal[i]);

	/* for random restarts */
	gsl_vector *loInit = gsl_vector_alloc(STATEDIM);
	gsl_vector *hiInit = gsl_vector_alloc(STATEDIM);
	gsl_vector_set_all(loInit, -4.5);
	gsl_vector_set_all(hiInit, 4.5);

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

    printf("Using minimizer %s.\n", gsl_multimin_fminimizer_name(s));
    if(!verbose) {
        printf(" %5s %7s %8s %8s %8s\n",
               "itr", "fval", "size", "vscale", "tmptr");
    }
	for(try=0; try < nTries; ++try) {
        int itr=0;   /* iteration within current temperature */
        int iT = 0;  /* index of current temperature */
		AnnealSched_reset(sched);
		if(try > 0)
			sasimplex_randomize_state(s, rotate, loInit, hiInit, ss);
		printf("%3d %5d %10.3e %10.3e\n",
			   try, itr,
			   gsl_vector_get(s->x, 0),
			   gsl_vector_get(s->x, 1));
        for(iT=0; iT<nT; ++iT) {  /* iterate over temperatures */
            AnnealSched_print(sched, stdout);
            temperature = AnnealSched_next(sched);
            status = sasimplex_n_iterations(s,
                                            &size,
                                            tol,
                                            nPerT,
                                            temperature,
                                            1);

            /* absErr is summed absolute error */
            double      errx = gsl_vector_get(s->x, 0) - par[0];
            double      erry = gsl_vector_get(s->x, 1) - par[1];
            absErr = fabs(errx) + fabs(erry);

            if(verbose) {
                printf("try=%d x=%.4f y=%.4f abserr=%.4e\n",
                       try, gsl_vector_get(s->x, 0),
                       gsl_vector_get(s->x, 1), absErr);
            }
            if(status != GSL_CONTINUE)
                break;
        }
        if(!verbose) {
            printf(" %5d %7.3f %8.3f %8.4f %8.4f\n",
                   itr, s->fval, size,
                   sasimplex_vertical_scale(s),
                   temperature);
            printf("try=%d x=%.4f y=%.4f abserr=%.4e\n",
                   try, gsl_vector_get(s->x, 0),
                   gsl_vector_get(s->x, 1), absErr);
        }
		switch (status) {
		case GSL_SUCCESS:
			printf("converged\n");
			break;
		default:
			printf("no convergence: status=%d\n", status);
		}
	}

	printf("True minimum: (%lf, %lf)\n", par[0], par[1]);
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    return status;
    }
