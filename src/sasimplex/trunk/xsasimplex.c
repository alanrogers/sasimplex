#include <stdio.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include "sasimplex.h"

typedef struct AnnealSched AnnealSched;

/**
 * This structure represents a schedule of annealing temperatures.
 */
struct AnnealSched {
	int nTemps;   /* number of temperature values */
	int nPerTemp; /* iterations per temperature value */
	int currTemp, currIteration;
	double *temp; /* array of temperature values */
};

AnnealSched *AnnealSched_alloc(int nTemps, int nPerTemp,
											   double initTemp,
											   double deflationFactor);
double AnnealSched_next(AnnealSched *s);

double my_f (const gsl_vector *v, void *params);

/**
 * Allocate and initialize annealing schedule for sasimplex.
 *
 * nTemps  : the number of temperatures
 * nPerTemp: the number of iterations at each temperature
 * initTemp: initial temperature
 * deflationFactor: the ratio of temperature i+1 to temperature i
 */
AnnealSched *
AnnealSched_alloc(int nTemps, int nPerTemp, double initTemp,
						 double deflationFactor) {
	int i;
	AnnealSched *s = malloc(sizeof(AnnealSched));
	if(s == NULL) {
#if 0
        GSL_ERROR("failed to allocate AnnealSched", GSL_ENOMEM);
#endif
		fprintf(stderr,"bad malloc\n");
		exit(1);
	}
	s->temp = malloc(nTemps * sizeof(s->temp[0]));
	if(s->temp == NULL) {
		free(s);
#if 0
        GSL_ERROR("failed to allocate annealing array", GSL_ENOMEM);
#endif
		fprintf(stderr,"bad malloc\n");
		exit(1);
	}

	s->currIteration = s->currTemp = 0;
	s->nTemps = nTemps;
	s->nPerTemp = nPerTemp;
	s->temp[0] = initTemp;
	for(i=1; i < nTemps-1; ++i)
		s->temp[i] = deflationFactor * s->temp[i-1];
	s->temp[nTemps-1] = 0.0;
	return s;
}

/** Get next temperature */
double AnnealSched_next(AnnealSched *s) {
	double currTemp;
	currTemp = s->temp[s->currTemp];

	++s->currIteration;
	if(s->currIteration == s->nPerTemp) {
		s->currIteration = 0;
		++s->currTemp;
		if(s->currTemp == s->nTemps)
			s->currTemp = s->nTemps - 1;
	}
	return currTemp;
}

/*
 * Multiple minima wherever x and y are integers. Global minimum
 * at (x,y)=(par[0],par[1]).
 */
double
my_f (const gsl_vector *v, void *params)
{
    double x, y;
    double fx, fy; /* fractional parts of x and y */
    double rval;
    double *par = (double *) params;
  
    x = gsl_vector_get(v, 0) - par[0];
    y = gsl_vector_get(v, 1) - par[1];
    fx = x - floor(x + 0.5);
    fy = y - floor(y + 0.5);

    rval = 1.0 + fabs(fx) + fabs(fy);
    rval *= 1.0 + fabs(x) + fabs(y);
    
    return rval; 
}

int 
main(void){
    double par[2] = {3.0, 2.0}; /* (x,y) at minumum */
    const size_t stateDim = 2; /**< dimension of state space */
    const double tol = 1e-4;
    const int maxItr = 1000;
    double initStepSize=0.05;

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

	/* Set up annealing schedule */
	AnnealSched *sched =
		AnnealSched_alloc(5,  /* number of temperatures */
								 100, /* iterations per temperature */
								 4.0,/* initial temperature */
								 0.5  /* ratio of adjacent temperatures */
			);

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
	unsigned long seed = time(NULL); /* get seed from clock */
    sasimplex_seed_rng(s, seed); 
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    printf("Using minimizer %s.\n", gsl_multimin_fminimizer_name(s));
    printf ("%5s %10s %10s %7s %8s %8s %8s\n", "itr",
            "x", "y", "f", "size", "MAE", "temp");
    do{
		double temperature;
		temperature = AnnealSched_next(sched);
		sasimplex_set_temp(s, temperature);
        itr++;
        status = gsl_multimin_fminimizer_iterate(s);
        if(status) {
            printf("%s:%d:%s: rtn val %d from %s\n",
                   __FILE__,__LINE__,__func__,
				   status,
				   "gsl_multimin_fminimizer_iterate");
            break;
        }

        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, tol);

		/* mae is mean absolute error */
		double errx = gsl_vector_get(s->x, 0) - par[0];
		double erry = gsl_vector_get(s->x, 1) - par[1];
		double mae = 0.5*(fabs(errx) + fabs(erry));
        printf ("%5d %10.3e %10.3e %7.3f %8.3f %8.4f %8.4f\n", 
                itr,
                gsl_vector_get(s->x, 0), 
                gsl_vector_get(s->x, 1), 
                s->fval, size, mae, temperature);
    }while (status == GSL_CONTINUE && itr < maxItr);
    switch(status) {
    case GSL_SUCCESS:
        printf("converged to minimum\n");
        break;
    default:
        printf("no convergence: status=%d itr=%d\n", status, itr);
    }
    printf("True minimum: (%lf, %lf)\n", par[0], par[1]);
 
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    return status;
}

