/* ========================================================================== */
/* === umfpack_scale ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

int umfpack_di_scale
(
    double X [ ],
    const double B [ ],
    void *Numeric
) ;

UF_long umfpack_dl_scale
(
    double X [ ],
    const double B [ ],
    void *Numeric
) ;

int umfpack_zi_scale
(
    double Xx [ ],	 double Xz [ ],
    const double Bx [ ], const double Bz [ ],
    void *Numeric
) ;

UF_long umfpack_zl_scale
(
    double Xx [ ],	 double Xz [ ],
    const double Bx [ ], const double Bz [ ],
    void *Numeric
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    double *B, *X ;
    status = umfpack_di_scale (X, B, Numeric) ;

double UF_long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    double *B, *X ;
    status = umfpack_dl_scale (X, B, Numeric) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    double *Bx, *Bz, *Xx, *Xz ;
    status = umfpack_zi_scale (Xx, Xz, Bx, Bz, Numeric) ;

complex UF_long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    double *Bx, *Bz, *Xx, *Xz ;
    status = umfpack_zl_scale (Xx, Xz, Bx, Bz, Numeric) ;

packed complex Syntax:

    Same as above, except both Xz and Bz are NULL.

Purpose:

    Given LU factors computed by umfpack_*_numeric (PAQ=LU, PRAQ=LU, or
    P(R\A)Q=LU), and a vector B, this routine computes X = B, X = R*B, or
    X = R\B, as appropriate.  X and B must be vectors equal in length to the
    number of rows of A.

Returns:

    The status code is returned.  UMFPACK_OK is returned if successful.
    UMFPACK_ERROR_invalid_Numeric_object is returned in the Numeric
    object is invalid.  UMFPACK_ERROR_argument_missing is returned if
    any of the input vectors are missing (X and B for the real version,
    and Xx and Bx for the complex version).

Arguments:

    double X [n_row] ;	Output argument.
    or:
    double Xx [n_row] ;	Output argument, real part.
			Size 2*n_row for packed complex case.
    double Xz [n_row] ;	Output argument, imaginary part.

	The output vector X.  If either Xz or Bz are NULL, the vector
	X is in packed complex form, with the kth entry in Xx [2*k] and
	Xx [2*k+1], and likewise for B.

    double B [n_row] ;	Input argument, not modified.
    or:
    double Bx [n_row] ;	Input argument, not modified, real part.
			Size 2*n_row for packed complex case.
    double Bz [n_row] ;	Input argument, not modified, imaginary part.

	The input vector B.  See above if either Xz or Bz are NULL.

    void *Numeric ;		Input argument, not modified.

	Numeric must point to a valid Numeric object, computed by
	umfpack_*_numeric.

*/
