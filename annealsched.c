#include <stdio.h>
#include <stdlib.h>
#include "annealsched.h"

/**
 * This structure represents a schedule of relative annealing
 * temperatures. The "relative temperature" is the ratio of absolute
 * temperature to some measure of vertical scale--of the variation in
 * function values.  The measure of scale is not part of AnnealSched,
 * but is provided as an argument to AnnealSched_next.
 */
struct AnnealSched {
    int         nT;                /* number of temperature values */
    int         nPerT;             /* iterations per temperature value */
    int         iT;                /* index of current tmptr */
    int         itr;               /* iteration w/i current tmptr */
	double      currRelT;          /* currRelT = tmptr/scale */
    double      initRelT;          /* initial value of currRelT */
    double      decay;             /* rate at which currRelT decays */
};

/**
 * Allocate and initialize annealing schedule for sasimplex.
 *
 * nT      : the number of temperatures
 * nPerT   : the number of iterations at each temperature
 * initRelT: initial relative temperature
 * decay   : the rate at which relT declines
 */
AnnealSched *AnnealSched_alloc(int nT, int nPerT, double initRelT,
                               double decay) {
    AnnealSched *s = malloc(sizeof(AnnealSched));
    if(s == NULL) {
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }
    s->initRelT = initRelT;
    s->nT = nT;
    s->nPerT = nPerT;
    s->decay = decay;
	AnnealSched_reset(s);
	return s;
}

void AnnealSched_reset(AnnealSched *s) {
    s->itr = 0;            /* range: 0..(nPerT-1) */
    s->iT = 0;             /* range: 0..(nT-1) */
    s->currRelT = s->initRelT;
}

/** Free memory allocated for annealing schedule */
void AnnealSched_free(AnnealSched * s) {
    free(s);
}

/** Get next temperature */
double AnnealSched_next(AnnealSched * s, double scale) {
    ++s->itr;
    if(s->itr == s->nPerT) {
        s->itr = 0;
        if(s->iT == s->nT - 1) {
            s->currRelT = 0.0;
        }else{
            ++s->iT;
            s->currRelT *= s->decay;
        }
    }
    return scale * s->currRelT;
}

void AnnealSched_print(AnnealSched *s, FILE *fp) {
    fprintf(fp, "itr=%d/%d iT=%d/%d\n", s->itr, s->nPerT,
            s->iT, s->nT);
    fprintf(fp, "relT: init=%lf curr=%lf decay=%lf\n",
            s->initRelT, s->currRelT, s->decay);
}
