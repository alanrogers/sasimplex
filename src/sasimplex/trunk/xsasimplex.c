#include "sasimplex.h"
#include "annealsched.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_multimin.h>

double      my_f(const gsl_vector * v, void *params);
int         pr_vector(FILE *fp, const char *fmt, const gsl_vector * v);

int pr_vector(FILE *fp, const char *fmt, const gsl_vector * v) {
    int rval=0;
    size_t i, n = v->size;
    putc('[', fp);
    if(n>0) {
        rval = fprintf(fp, fmt, gsl_vector_get(v, 0));
        if(rval <= 0)
            return rval;
        
    }
    for(i=1; i < n; ++i) {
        fputs(", ", fp);
        rval = fprintf(fp,fmt,  gsl_vector_get(v, i));
        if(rval <= 0)
            return rval;
    }
    putc(']', fp);
    return rval;
}

/*
 * Local minima wherever x-par[0] and y-par[1] are integers. Global
 * minimum at (x,y)=(par[0],par[1]).
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
	const int   verbose = 0;
    unsigned long seed = time(NULL);    /* for random numbers */
	unsigned    nTries = 20;            /* number of random starts */
    int         nT  = 10;               /* number of temperatures */
    int         nPerT = 200;            /* iterations per temperature */
    double      initT = 3.0;            /* initial temperature */
    double      decay = 0.7;

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

    const gsl_multimin_fminimizer_type *fmType 
        = gsl_multimin_fminimizer_sasimplex;
    gsl_multimin_fminimizer *s 
        = gsl_multimin_fminimizer_alloc(fmType, STATEDIM);
    if(s == NULL) {
        fprintf(stderr, "%s:%d: bad allocation\n", __FILE__, __LINE__);
        exit(1);
    }

    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
    sasimplex_random_seed(s, seed);

    printf("Using minimizer %s.\n", gsl_multimin_fminimizer_name(s));
	for(try=0; try < nTries; ++try) {
        int iT = 0;  /* index of current temperature */
		AnnealSched_reset(sched);
		if(try > 0)
			sasimplex_randomize_state(s, rotate, loInit, hiInit, ss);
        for(iT=0; iT<nT; ++iT) {  /* iterate over temperatures */
            temperature = AnnealSched_next(sched);
            status = sasimplex_n_iterations(s,
                                            &size,
                                            tol,
                                            nPerT,
                                            temperature,
                                            verbose);

            /* absErr is summed absolute error */
            double      errx = gsl_vector_get(s->x, 0) - par[0];
            double      erry = gsl_vector_get(s->x, 1) - par[1];
            absErr = fabs(errx) + fabs(erry);

            if(verbose) {
                printf("try=%d x=%.4lf y=%.4f abserr=%.4le\n",
                       try, gsl_vector_get(s->x, 0),
                       gsl_vector_get(s->x, 1), absErr);
            }
            if(status != GSL_CONTINUE)
                break;
        }
        if(!verbose) {
            printf("%2d: size=%.4le aberr=%.4le vscale=%.4lf x=",
                   try, size, absErr,
                   sasimplex_vertical_scale(s));
            pr_vector(stdout, "%.4lf", s->x);
        }
		switch (status) {
		case GSL_SUCCESS:
			printf(" converged\n");
			break;
        case GSL_CONTINUE:
			printf(" no convergence in %d iterations at tmptr %lf\n",
                   nPerT, temperature);
            sasimplex_print(s);
			break;
		default:
			printf(" unknown status: %d\n", status);
		}
	}

	printf("True minimum: (%lf, %lf)\n", par[0], par[1]);
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);

    return status;
    }
