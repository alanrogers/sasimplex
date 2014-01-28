static inline double ran_unif(unsigned long *seed);
static int sasimplex_iterate(void *vstate, gsl_multimin_function * f,
                             gsl_vector * x, double *size, double *fval,
                             unsigned long *seed);

static inline double 
ran_unif(unsigned long *seed) {
    unsigned long s = *seed;

    *seed = (s * 69069 + 1) & 0xffffffffUL;
    return (*seed) / 4294967296.0;
}

static int
sasimplex_iterate(void *vstate, gsl_multimin_function * f,
                  gsl_vector * x, double *size, double *fval,
                  unsigned long *seed) {

    /* Simplex iteration tries to minimize function f value */
    /* Includes corrections from Ivo Alxneit <ivo.alxneit@psi.ch> */

    nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

    /* xc and xc2 vectors store tried corner point coordinates */

    gsl_vector *xc = state->ws1;
    gsl_vector *xc2 = state->ws2;
    gsl_vector *y1 = state->y1;
    gsl_matrix *x1 = state->x1;

    const size_t n = y1->size;
    size_t      i;
    size_t      hi, s_hi, lo;
    double      dhi, ds_hi, dlo;
    int         status;
    double      val, val2;

    if(xc->size != x->size) {
        GSL_ERROR("incompatible size of x", GSL_EINVAL);
    }

    /* get index of highest, second highest and lowest point */

    dhi = dlo = gsl_vector_get(y1, 0) + temp*log(rand(&seed));
    hi = 0;
    lo = 0;

    ds_hi = gsl_vector_get(y1, 1) + temp*log(rand(&seed));
    s_hi = 1;

    for(i = 1; i < n; i++) {
        val = gsl_vector_get(y1, i) + temp*log(rand(&seed));
        if(val < dlo) {
            dlo = val;
            lo = i;
        } else if(val > dhi) {
            ds_hi = dhi;
            s_hi = hi;
            dhi = val;
            hi = i;
        } else if(val > ds_hi) {
            ds_hi = val;
            s_hi = i;
        }
    }

    /* try reflecting the highest value point */

    val = try_corner_move(-1.0, state, hi, xc, f);

    if(gsl_finite(val) && val < 
       gsl_vector_get(y1, lo) + temp*log(rand(&seed))) {
        /* reflected point is lowest, try expansion */

        val2 = try_corner_move(-2.0, state, hi, xc2, f);

        if(gsl_finite(val2) && val2 < 
           gsl_vector_get(y1, lo) + temp*log(rand(&seed))) {
            update_point(state, hi, xc2, val2);
        } else {
            update_point(state, hi, xc, val);
        }
    } else if(!gsl_finite(val) || val >
              gsl_vector_get(y1, s_hi) + temp*log(rand(&seed))) {
        /* reflection does not improve things enough, or we got a
         * non-finite function value */

        if(gsl_finite(val) && val <=
           gsl_vector_get(y1, hi) + temp*log(rand(&seed))) {
            /* if trial point is better than highest point, replace
             * highest point */

            update_point(state, hi, xc, val);
        }

        /* try one-dimensional contraction */

        val2 = try_corner_move(0.5, state, hi, xc2, f);

        if(gsl_finite(val2) && val2 <=
           gsl_vector_get(y1, hi) + temp*log(rand(&seed))) {
            update_point(state, hi, xc2, val2);
        } else {
            /* contract the whole simplex about the best point */

            status = contract_by_best(state, lo, xc, f);

            if(status != GSL_SUCCESS) {
                GSL_ERROR("contraction failed", GSL_EFAILED);
            }
        }
    } else {
        /* trial point is better than second highest point.  Replace
         * highest point by it */

        update_point(state, hi, xc, val);
    }

    /* return lowest point of simplex as x */

    lo = gsl_vector_min_index(y1);
    gsl_matrix_get_row(x, x1, lo);
    *fval = gsl_vector_get(y1, lo);

    /* Update simplex size */

    {
        double      S2 = state->S2;

        if(S2 > 0) {
            *size = sqrt(S2);
        } else {
            /* recompute if accumulated error has made size invalid */
            *size = compute_size(state, state->center);
        }
    }

    return GSL_SUCCESS;
}
