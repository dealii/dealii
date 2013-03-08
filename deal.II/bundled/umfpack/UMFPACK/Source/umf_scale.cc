/* ========================================================================== */
/* === UMF_scale ============================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* Divide a vector of stride 1 by the pivot value. */

#include "umf_internal.h"

GLOBAL void UMF_scale
(
    Int n,
    Entry pivot,
    Entry X [ ]
)
{
    Entry x ;
    double s ;
    Int i ;

    /* ---------------------------------------------------------------------- */
    /* compute the approximate absolute value of the pivot, and select method */
    /* ---------------------------------------------------------------------- */

    APPROX_ABS (s, pivot) ;

    if (s < RECIPROCAL_TOLERANCE || IS_NAN (pivot))
    {
	/* ------------------------------------------------------------------ */
	/* tiny, or zero, pivot case */
	/* ------------------------------------------------------------------ */

	/* The pivot is tiny, or NaN.  Do not divide zero by the pivot value,
	 * and do not multiply by 1/pivot, either. */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] /= pivot ; */
	    x = X [i] ;

#ifndef NO_DIVIDE_BY_ZERO
	    if (IS_NONZERO (x))
	    {
		DIV (X [i], x, pivot) ;
	    }
#else
	    /* Do not divide by zero */
	    if (IS_NONZERO (x) && IS_NONZERO (pivot))
	    {
		DIV (X [i], x, pivot) ;
	    }
#endif

	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* normal case */
	/* ------------------------------------------------------------------ */

	/* The pivot is not tiny, and is not NaN.   Don't bother to check for
	 * zeros in the pivot column, X.  This is slightly more accurate than
	 * multiplying by 1/pivot (but slightly slower), particularly if the
	 * pivot column consists of only IEEE subnormals. */

	for (i = 0 ; i < n ; i++)
	{
	    /* X [i] /= pivot ; */
	    x = X [i] ;
	    DIV (X [i], x, pivot) ;
	}
    }
}
