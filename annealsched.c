#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "annealsched.h"

/**
 * This structure represents a geometric schedule of annealing
 * temperatures.
 */
struct AnnealSched {
    int         nT;                /* number of temperature values */
    int         iT;                /* index of current tmptr */
	double      T;                 /* current temperature */
    double      initT;             /* initial temperature */
    double      decay;             /* rate at which temperature decays */
};

/**
 * Allocate and initialize annealing schedule for sasimplex.
 *
 * nT      : the number of temperatures
 * initT: initial relative temperature
 * decay   : the rate at which T declines
 */
AnnealSched *AnnealSched_alloc(int nT, double initT, double decay) {
    AnnealSched *s = malloc(sizeof(AnnealSched));
    if(s == NULL) {
        fprintf(stderr, "bad malloc\n");
        exit(1);
    }
    if(nT == 1)
        initT = 0.0;
    s->initT = initT;
    s->decay = decay;
    s->nT = nT;
	AnnealSched_reset(s);
	return s;
}

/** Return number of temperature values */
int AnnealSched_size(const AnnealSched *sched) {
    return sched->nT;
}

void AnnealSched_reset(AnnealSched *s) {
    s->iT = 0;             /* range: 0..(nT-1) */
    s->T = s->initT;
}

/** Free memory allocated for annealing schedule */
void AnnealSched_free(AnnealSched * s) {
    free(s);
}

AnnealSched * AnnealSched_copy(const AnnealSched *old) {
    AnnealSched *new = malloc(sizeof(*new));
    if(new == NULL) {
        fprintf(stderr,"%s%d%s: bad malloc\n",
                __FILE__,__LINE__,__func__);
        exit(ENOMEM);
    }
    assert(sizeof(*old) == sizeof(*new));
    memcpy(new, old, sizeof(*new));
    return new;
}

int AnnealSched_cmp(const AnnealSched *s1, const AnnealSched *s2) {
    return memcmp( (const void *) s1, (const void *) s2,
                   sizeof(AnnealSched));
}

/** Get next temperature */
double AnnealSched_next(AnnealSched * s) {
    double tmptr = s->T;

    ++s->iT;
    if(s->iT == s->nT-1) 
        s->T = 0.0;
    else
        s->T *= s->decay;

    return tmptr;
}

void AnnealSched_print(AnnealSched *s, FILE *fp) {
    fprintf(fp, "%s: iT=%d/%d",  __func__, s->iT, s->nT);
    fprintf(fp, " initT=%lf currT=%lf\n",
            s->initT, s->T);
}
