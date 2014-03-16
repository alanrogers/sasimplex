#include <stdio.h>
#include <stdlib.h>
#include "annealsched.h"

/**
 * This structure represents a schedule of annealing temperatures.
 */
struct AnnealSched {
    int         nTmptrs;         /* number of temperature values */
    int         nPerTmptr;       /* iterations per temperature value */
	double      initTmptr;       /* initial temperature */
	double      param;           /* controls rate of temperature decline */
    int         iTmptr, itr;
    double     *tmptr;           /* array of temperature values */
};

#define LINEAR_ANNEALING_SCHEDULE

#ifdef GEOMETRIC_ANNEALING_SCHEDULE
/**
 * Allocate and initialize annealing schedule for sasimplex.
 *
 * nTmptrs  : the number of temperatures
 * nPerTmptr: the number of iterations at each temperature
 * initTmptr: initial temperature
 * deflationFactor: the ratio of temperature i+1 to temperature i
 */
AnnealSched *AnnealSched_alloc(int nTmptrs, int nPerTmptr, double initTmptr,
                               double deflationFactor) {
	int i;
    AnnealSched *s = malloc(sizeof(AnnealSched));
    if(s == NULL) {
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }
    s->tmptr = malloc(nTmptrs * sizeof(s->tmptr[0]));
    if(s->tmptr == NULL) {
        free(s);
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }

    s->nTmptrs = nTmptrs;
    s->nPerTmptr = nPerTmptr;
	s->initTmptr = initTmptr;
	s->param = deflationFactor;
    s->tmptr[0] = s->initTmptr;
    for(i = 1; i < s->nTmptrs - 1; ++i)
        s->tmptr[i] = s->param * s->tmptr[i - 1];
    s->tmptr[s->nTmptrs - 1] = 0.0;
	AnnealSched_reset(s);
	return s;
}
#elif defined(LINEAR_ANNEALING_SCHEDULE)
/**
 * Allocate and initialize annealing schedule for sasimplex.
 *
 * nTmptrs  : the number of temperatures
 * nPerTmptr: the number of iterations at each temperature
 * initTmptr: initial temperature
 * deflationFactor: the ratio of temperature i+1 to temperature i
 */
AnnealSched *AnnealSched_alloc(int nTmptrs, int nPerTmptr, double initTmptr,
                               double deflationFactor) {
	int i;
    AnnealSched *s = malloc(sizeof(AnnealSched));
    if(s == NULL) {
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }
    s->tmptr = malloc(nTmptrs * sizeof(s->tmptr[0]));
    if(s->tmptr == NULL) {
        free(s);
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }

    double  step = initTmptr/ (nTmptrs - 1);
    s->nTmptrs = nTmptrs;
    s->nPerTmptr = nPerTmptr;
	s->initTmptr = initTmptr;
	s->param = deflationFactor;
    s->tmptr[0] = s->initTmptr;
    for(i = 0; i < s->nTmptrs - 1; ++i)
        s->tmptr[i] = initTmptr - i*step;
    s->tmptr[s->nTmptrs - 1] = 0.0;
	AnnealSched_reset(s);
	return s;
}
#else
#  error "ANNEALING_SCHEDULE undefined"
#endif

void AnnealSched_reset(AnnealSched *s) {
    s->itr = s->iTmptr = 0;
}

/** Free memory allocated for annealing schedule */
void AnnealSched_free(AnnealSched * s) {
    free(s->tmptr);
    free(s);
}

/** Get next temperature */
double AnnealSched_next(AnnealSched * s) {
    double      currTmptr;
    currTmptr = s->tmptr[s->iTmptr];
    ++s->itr;
    if(s->itr == s->nPerTmptr) {
        s->itr = 0;
        ++s->iTmptr;
        if(s->iTmptr == s->nTmptrs)
            s->iTmptr = s->nTmptrs - 1;
    }
    return currTmptr;
}
