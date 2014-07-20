/* ========================================================================== */
/* === UMFPACK_get_determinant ============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* UMFPACK_get_determinant contributed by David Bateman, Motorola, Paris. */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  From the LU factors, scale factor, and permutation vectors
    held in the Numeric object, calculates the determinant of the matrix A.
    See umfpack_get_determinant.h for a more detailed description.

    Dynamic memory usage:  calls UMF_malloc once, for a total space of
    n integers, and then frees all of it via UMF_free when done.

    Contributed by David Bateman, Motorola, Nov. 2004.
    Modified for V4.4, Jan. 2005.
*/

#include "umf_internal.h"
#include "umf_valid_numeric.h"
#include "umf_malloc.h"
#include "umf_free.h"

/* ========================================================================== */
/* === rescale_determinant ================================================== */
/* ========================================================================== */

/* If the mantissa is too big or too small, rescale it and change exponent */

PRIVATE Int rescale_determinant
(
    Entry *d_mantissa,
    double *d_exponent
)
{
    double d_abs ;

    ABS (d_abs, *d_mantissa) ;

    if (SCALAR_IS_ZERO (d_abs))
    {
	/* the determinant is zero */
	*d_exponent = 0 ;
	return (FALSE) ;
    }

    if (SCALAR_IS_NAN (d_abs))
    {
	/* the determinant is NaN */
	return (FALSE) ;
    }

    while (d_abs < 1.)
    {
	SCALE (*d_mantissa, 10.0) ;
	*d_exponent = *d_exponent - 1.0 ;
	ABS (d_abs, *d_mantissa) ;
    }

    while (d_abs >= 10.)
    {
	SCALE (*d_mantissa, 0.1) ;
	*d_exponent = *d_exponent + 1.0 ;
	ABS (d_abs, *d_mantissa) ;
    }

    return (TRUE) ;
}

/* ========================================================================== */
/* === UMFPACK_get_determinant ============================================== */
/* ========================================================================== */

GLOBAL Int UMFPACK_get_determinant
(
    double *Mx,
#ifdef COMPLEX
    double *Mz,
#endif
    double *Ex,
    void *NumericHandle,
    double User_Info [UMFPACK_INFO]
)
{
    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    Entry d_mantissa, d_tmp ;
    double d_exponent, Info2 [UMFPACK_INFO], one [2] = {1.0, 0.0}, d_sign ;
    Entry *D ;
    double *Info, *Rs ;
    NumericType *Numeric ;
    Int i, n, itmp, npiv, *Wi, *Rperm, *Cperm, do_scale ;

#ifndef NRECIPROCAL
    Int do_recip ;
#endif

    /* ---------------------------------------------------------------------- */
    /* check input parameters */
    /* ---------------------------------------------------------------------- */

    if (User_Info != (double *) NULL)
    {
	/* return Info in user's array */
	Info = User_Info ;
    }
    else
    {
	/* no Info array passed - use local one instead */
	Info = Info2 ;
	for (i = 0 ; i < UMFPACK_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
    }

    Info [UMFPACK_STATUS] = UMFPACK_OK ;

    Numeric = (NumericType *) NumericHandle ;
    if (!UMF_valid_numeric (Numeric))
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_Numeric_object ;
	return (UMFPACK_ERROR_invalid_Numeric_object) ;
    }

    if (Numeric->n_row != Numeric->n_col)
    {
	/* only square systems can be handled */
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_invalid_system ;
	return (UMFPACK_ERROR_invalid_system) ;
    }

    if (Mx == (double *) NULL)
    {
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_argument_missing ;
	return (UMFPACK_ERROR_argument_missing) ;
    }

    n = Numeric->n_row ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    Wi = (Int *) UMF_malloc (n, sizeof (Int)) ;

    if (!Wi)
    {
	DEBUGm4 (("out of memory: get determinant\n")) ;
	Info [UMFPACK_STATUS] = UMFPACK_ERROR_out_of_memory ;
	return (UMFPACK_ERROR_out_of_memory) ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute the determinant */
    /* ---------------------------------------------------------------------- */

    Rs = Numeric->Rs ;		/* row scale factors */
    do_scale = (Rs != (double *) NULL) ;

#ifndef NRECIPROCAL
    do_recip = Numeric->do_recip ;
#endif

    d_mantissa = ((Entry *) one) [0] ;
    d_exponent = 0.0 ;
    D = Numeric->D ;

    /* compute product of diagonal entries of U */
    for (i = 0 ; i < n ; i++)
    {
	MULT (d_tmp, d_mantissa, D [i]) ;
	d_mantissa = d_tmp ;

	if (!rescale_determinant (&d_mantissa, &d_exponent))
	{
	    /* the determinant is zero or NaN */
	    Info [UMFPACK_STATUS] = UMFPACK_WARNING_singular_matrix ;
	    /* no need to compute the determinant of R */
	    do_scale = FALSE ;
	    break ;
	}
    }

    /* compute product of diagonal entries of R (or its inverse) */
    if (do_scale)
    {
	for (i = 0 ; i < n ; i++)
	{
#ifndef NRECIPROCAL
	    if (do_recip)
	    {
		/* compute determinant of R inverse */
		SCALE_DIV (d_mantissa, Rs [i]) ;
	    }
	    else
#endif
	    {
		/* compute determinant of R */
		SCALE (d_mantissa, Rs [i]) ;
	    }
	    if (!rescale_determinant (&d_mantissa, &d_exponent))
	    {
		/* the determinant is zero or NaN.  This is very unlikey to
		 * occur here, since the scale factors for a tiny or zero row
		 * are set to 1. */
		Info [UMFPACK_STATUS] = UMFPACK_WARNING_singular_matrix ;
		break ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* determine if P and Q are odd or even permutations */
    /* ---------------------------------------------------------------------- */

    npiv = 0 ;
    Rperm = Numeric->Rperm ;

    for (i = 0 ; i < n ; i++)
    {
	Wi [i] = Rperm [i] ;
    }

    for (i = 0 ; i < n ; i++)
    {
	while (Wi [i] != i)
	{
	    itmp = Wi [Wi [i]] ;
	    Wi [Wi [i]] = Wi [i] ;
	    Wi [i] = itmp ;
	    npiv++ ;
	}
    }

    Cperm = Numeric->Cperm ;

    for (i = 0 ; i < n ; i++)
    {
	Wi [i] = Cperm [i] ;
    }

    for (i = 0 ; i < n ; i++)
    {
	while (Wi [i] != i)
	{
	    itmp = Wi [Wi [i]] ;
	    Wi [Wi [i]] = Wi [i] ;
	    Wi [i] = itmp ;
	    npiv++ ;
	}
    }

    /* if npiv is odd, the sign is -1.  if it is even, the sign is +1 */
    d_sign = (npiv % 2) ? -1. : 1. ;

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    (void) UMF_free ((void *) Wi) ;

    /* ---------------------------------------------------------------------- */
    /* compute the magnitude and exponent of the determinant */
    /* ---------------------------------------------------------------------- */

    if (Ex == (double *) NULL)
    {
	/* Ex is not provided, so return the entire determinant in d_mantissa */
	SCALE (d_mantissa, pow (10.0, d_exponent)) ;
    }
    else
    {
	Ex [0] = d_exponent ;
    }

    Mx [0] = d_sign * REAL_COMPONENT (d_mantissa) ;

#ifdef COMPLEX
    if (SPLIT (Mz))
    {
	Mz [0] = d_sign * IMAG_COMPONENT (d_mantissa) ;
    }
    else
    {
	Mx [1] = d_sign * IMAG_COMPONENT (d_mantissa) ;
    }
#endif

    /* determine if the determinant has (or will) overflow or underflow */
    if (d_exponent + 1.0 > log10 (DBL_MAX))
    {
	Info [UMFPACK_STATUS] = UMFPACK_WARNING_determinant_overflow ;
    }
    else if (d_exponent - 1.0 < log10 (DBL_MIN))
    {
	Info [UMFPACK_STATUS] = UMFPACK_WARNING_determinant_underflow ;
    }

    return (UMFPACK_OK) ;
}
