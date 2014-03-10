#include <stdio.h>
#include <stdlib.h>
#include "annealsched.h"


/**
 * This structure represents a schedule of annealing temperatures.
 */
struct AnnealSched {
    int         nTemps;         /* number of temperature values */
    int         nPerTemp;       /* iterations per temperature value */
    int         currTemp, currIteration;
    double     *temp;           /* array of temperature values */
};

/**
 * Allocate and initialize annealing schedule for sasimplex.
 *
 * nTemps  : the number of temperatures
 * nPerTemp: the number of iterations at each temperature
 * initTemp: initial temperature
 * deflationFactor: the ratio of temperature i+1 to temperature i
 */
AnnealSched *AnnealSched_alloc(int nTemps, int nPerTemp, double initTemp,
                               double deflationFactor) {
    int         i;

    AnnealSched *s = malloc(sizeof(AnnealSched));

    if(s == NULL) {
#if 0
        GSL_ERROR("failed to allocate AnnealSched", GSL_ENOMEM);
#endif
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }
    s->temp = malloc(nTemps * sizeof(s->temp[0]));
    if(s->temp == NULL) {
        free(s);
#if 0
        GSL_ERROR("failed to allocate annealing array", GSL_ENOMEM);
#endif
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }

    s->currIteration = s->currTemp = 0;
    s->nTemps = nTemps;
    s->nPerTemp = nPerTemp;
    s->temp[0] = initTemp;
    for(i = 1; i < nTemps - 1; ++i)
        s->temp[i] = deflationFactor * s->temp[i - 1];
    s->temp[nTemps - 1] = 0.0;
    return s;
}

/** Free memory allocated for annealing schedule */
void AnnealSched_free(AnnealSched * s) {
    free(s->temp);
    free(s);
}

/** Get next temperature */
double AnnealSched_next(AnnealSched * s) {
    double      currTemp;

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
