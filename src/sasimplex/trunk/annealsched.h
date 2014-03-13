#ifndef __ARR_ANNEALSCHED__
#define __ARR_ANNEALSCHED__

typedef struct AnnealSched AnnealSched;

AnnealSched *AnnealSched_alloc(int nTmptrs, int nPerTmptr,
                               double initTmptr, double deflationFactor);
double      AnnealSched_next(AnnealSched * s);
void        AnnealSched_free(AnnealSched * s);

#endif /* __ARR_ANNEALSCHED__ */
