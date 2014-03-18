#ifndef __ARR_ANNEALSCHED__
#define __ARR_ANNEALSCHED__

typedef struct AnnealSched AnnealSched;

#include <stdio.h>

AnnealSched *AnnealSched_alloc(int nT, double initT, double decay);
void         AnnealSched_reset(AnnealSched *s);
double       AnnealSched_next(AnnealSched * s);
void         AnnealSched_free(AnnealSched * s);
void         AnnealSched_print(AnnealSched *s, FILE *fp);
#endif /* __ARR_ANNEALSCHED__ */
