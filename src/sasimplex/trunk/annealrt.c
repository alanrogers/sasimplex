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
    int         iT;                /* index of current tmptr */
	double      relT;              /* relT = tmptr/scale */
    double      initRelT;          /* initial value of currRelT */
};

/**
 * Allocate and initialize annealing schedule for sasimplex.
 *
 * nT      : the number of temperatures
 * initRelT: initial relative temperature
 */
AnnealSched *AnnealSched_alloc(int nT, double initRelT) {
    AnnealSched *s = malloc(sizeof(AnnealSched));
    if(s == NULL) {
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }
    s->initRelT = initRelT;
    s->nT = nT;
	AnnealSched_reset(s);
	return s;
}

void AnnealSched_reset(AnnealSched *s) {
    s->iT = 0;             /* range: 0..(nT-1) */
    s->relT = s->initRelT;
}

/** Free memory allocated for annealing schedule */
void AnnealSched_free(AnnealSched * s) {
    free(s);
}

/** Get next temperature */
double AnnealSched_next(AnnealSched * s, double scale) {
    if(++s->iT == s->nT) {
        s->iT = 0;
        s->relT = 0.0;
    }
    return scale * s->relT;
}

void AnnealSched_print(AnnealSched *s, FILE *fp) {
    fprintf(fp, "%s: iT=%d/%d\n",  __func__, s->iT, s->nT);
    fprintf(fp, "relT: init=%lf curr=%lf\n",
            s->initRelT, s->relT);
}
