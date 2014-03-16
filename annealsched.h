#ifndef __ARR_ANNEALSCHED__
#define __ARR_ANNEALSCHED__

typedef struct AnnealSched AnnealSched;

AnnealSched *AnnealSched_alloc(int nTmptrs, int nPerTmptr,
                               double initTmptr, double deflationFactor);
void        AnnealSched_reset(AnnealSched *s);
double      AnnealSched_next(AnnealSched * s, double scale);
void        AnnealSched_free(AnnealSched * s);

#endif /* __ARR_ANNEALSCHED__ */
