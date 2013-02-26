/* ========================================================================== */
/* === UMFPACK_global ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    Global variables.  UMFPACK uses these function pointers for several
    user-redefinable functions.   The amd_* functions are defined in
    AMD/Source/amd_global.c.

    Function pointer	    default	    for mexFunction
					    (see MATLAB/umfpackmex.c)
    ----------------	    -------	    ---------------
    amd_malloc		    malloc	    mxMalloc
    amd_free		    free	    mxFree
    amd_realloc		    realloc	    mxRealloc
    amd_calloc		    calloc	    mxCalloc
    amd_printf		    printf	    mexPrintf

    umfpack_hypot	    umf_hypot	    umf_hypot
    umfpack_divcomplex	    umf_divcomplex  umf_divcomplex

    This routine is compiled just once for all four versions of UMFPACK
    (int/UF_long, double/complex).
*/

#include "umf_internal.h"

double (*umfpack_hypot) (double, double) = umf_hypot ;
int (*umfpack_divcomplex) (double, double, double, double, double *, double *)
    = umf_divcomplex ;


/* ========================================================================== */
/* === umf_hypot ============================================================ */
/* ========================================================================== */

/* There is an equivalent routine called hypot in <math.h>, which conforms
 * to ANSI C99.  However, UMFPACK does not assume that ANSI C99 is available.
 * You can use the ANSI C99 hypot routine with:
 *
 *	#include <math.h>
 *	umfpack_hypot = hypot ;
 *
 * prior to calling any UMFPACK routine.
 *
 * s = hypot (x,y) computes s = sqrt (x*x + y*y) but does so more accurately.
 *
 * The NaN case for the double relops x >= y and x+y == x is safely ignored.
 */

double umf_hypot (double x, double y)
{
    double s, r ;
    x = SCALAR_ABS (x) ;
    y = SCALAR_ABS (y) ;
    if (x >= y)
    {
	if (x + y == x)
	{
	    s = x ;
	}
	else
	{
	    r = y / x ;
	    s = x * sqrt (1.0 + r*r) ;
	}
    }
    else
    {
	if (y + x == y)
	{
	    s = y ;
	}
	else
	{
	    r = x / y ;
	    s = y * sqrt (1.0 + r*r) ;
	}
    } 
    return (s) ;
}


/* ========================================================================== */
/* === umf_divcomplex ======================================================= */
/* ========================================================================== */

/* c = a/b where c, a, and b are complex.  The real and imaginary parts are
 * passed as separate arguments to this routine.  The NaN case is ignored
 * for the double relop br >= bi.  Returns TRUE (1) if the denominator is
 * zero, FALSE (0) otherwise.
 *
 * This uses ACM Algo 116, by R. L. Smith, 1962, which tries to avoid
 * underflow and overflow.
 *
 * c can be the same variable as a or b.
 */

int umf_divcomplex
(
    double ar, double ai,	/* real and imaginary parts of a */
    double br, double bi,	/* real and imaginary parts of b */
    double *cr, double *ci	/* real and imaginary parts of c */
)
{
    double tr, ti, r, den ;
    if (SCALAR_ABS (br) >= SCALAR_ABS (bi))
    {
	r = bi / br ;
	den = br + r * bi ;
	tr = (ar + ai * r) / den ;
	ti = (ai - ar * r) / den ;
    }
    else
    {
	r = br / bi ;
	den = r * br + bi ;
	tr = (ar * r + ai) / den ;
	ti = (ai * r - ar) / den ;
    }
    *cr = tr ;
    *ci = ti ;
    return (SCALAR_IS_ZERO (den)) ;
}
