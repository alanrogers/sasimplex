#include <stdio.h>
#include <stdlib.h>
#include "annealsched.h"

/**
 * This structure represents a schedule of annealing temperatures.
 */
struct AnnealSched {
    int         nTmptrs;         /* number of temperature values */
    int         nPerTmptr;       /* iterations per temperature value */
    int         currTmptr, currIteration;
    double     *tmptr;           /* array of temperature values */
};

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
    int         i;
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

    s->currIteration = s->currTmptr = 0;
    s->nTmptrs = nTmptrs;
    s->nPerTmptr = nPerTmptr;
    s->tmptr[0] = initTmptr;
    for(i = 1; i < nTmptrs - 1; ++i)
        s->tmptr[i] = deflationFactor * s->tmptr[i - 1];
    s->tmptr[nTmptrs - 1] = 0.0;
    return s;
}

/** Free memory allocated for annealing schedule */
void AnnealSched_free(AnnealSched * s) {
    free(s->tmptr);
    free(s);
}

/** Get next temperature */
double AnnealSched_next(AnnealSched * s) {
    double      currTmptr;
    currTmptr = s->tmptr[s->currTmptr];
    ++s->currIteration;
    if(s->currIteration == s->nPerTmptr) {
        s->currIteration = 0;
        ++s->currTmptr;
        if(s->currTmptr == s->nTmptrs)
            s->currTmptr = s->nTmptrs - 1;
    }
    return currTmptr;
}
