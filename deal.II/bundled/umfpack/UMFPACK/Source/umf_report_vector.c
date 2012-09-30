/* ========================================================================== */
/* === UMF_report_vector ==================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#include "umf_internal.h"

/* ========================================================================== */
/* === print_value ========================================================== */
/* ========================================================================== */

PRIVATE void print_value
(
    Int i,
    const double Xx [ ],
    const double Xz [ ],    /* used for complex case only */
    Int scalar		    /* if true, then print real part only */
)
{
    Entry xi ;
    /* if Xz is null, then X is in "merged" format (compatible with Entry, */
    /* and ANSI C99 double _Complex type). */
    PRINTF (("    "ID" :", INDEX (i))) ;
    if (scalar)
    {
	PRINT_SCALAR (Xx [i]) ;
    }
    else
    {
	ASSIGN (xi, Xx, Xz, i, SPLIT (Xz)) ;
	PRINT_ENTRY (xi) ;
    }
    PRINTF (("\n")) ;
}

/* ========================================================================== */
/* === UMF_report_vector ==================================================== */
/* ========================================================================== */

GLOBAL Int UMF_report_vector
(
    Int n,
    const double Xx [ ],
    const double Xz [ ],
    Int prl,
    Int user,
    Int scalar
)
{
    Int n2, i ;

    if (user || prl >= 4)
    {
	PRINTF (("dense vector, n = "ID". ", n)) ;
    }

    if (user)
    {
	if (!Xx)
	{
	    PRINTF (("ERROR: vector not present\n\n")) ;
	    return (UMFPACK_ERROR_argument_missing) ;
	}
	if (n < 0)
	{
	    PRINTF (("ERROR: length of vector is < 0\n\n")) ;
	    return (UMFPACK_ERROR_n_nonpositive) ;
	}
    }

    if (user || prl >= 4)
    {
	PRINTF4 (("\n")) ;
    }

    if (prl == 4)
    {
	/* print level of 4 */
	n2 = MIN (10, n) ;
	for (i = 0 ; i < n2 ; i++)
	{
	    print_value (i, Xx, Xz, scalar) ;
	}
	if (n2 < n)
	{
	    PRINTF (("    ...\n")) ;
	    print_value (n-1, Xx, Xz, scalar) ;
	}
    }
    else if (prl > 4)
    {
	/* print level 4 or more */
	for (i = 0 ; i < n ; i++)
	{
	    print_value  (i, Xx, Xz, scalar) ;
	}
    }

    PRINTF4 (("    dense vector ")) ;
    if (user || prl >= 4)
    {
	PRINTF (("OK\n\n")) ;
    }
    return (UMFPACK_OK) ;
}
