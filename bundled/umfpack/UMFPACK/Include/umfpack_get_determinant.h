/* ========================================================================== */
/* === UMFPACK_get_determinant ============================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* UMFPACK_get_determinant contributed by David Bateman, Motorola, Paris. */
/* -------------------------------------------------------------------------- */

int umfpack_di_get_determinant
(
    double *Mx,
    double *Ex,
    void *NumericHandle,
    double User_Info [UMFPACK_INFO]
) ;

UF_long umfpack_dl_get_determinant
(
    double *Mx,
    double *Ex,
    void *NumericHandle,
    double User_Info [UMFPACK_INFO]
) ;

int umfpack_zi_get_determinant
(
    double *Mx,
    double *Mz,
    double *Ex,
    void *NumericHandle,
    double User_Info [UMFPACK_INFO]
) ;

UF_long umfpack_zl_get_determinant
(
    double *Mx,
    double *Mz,
    double *Ex,
    void *NumericHandle,
    double User_Info [UMFPACK_INFO]
) ;

/*
double int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status ;
    double Mx, Ex, Info [UMFPACK_INFO] ;
    status = umfpack_di_get_determinant (&Mx, &Ex, Numeric, Info) ;

double UF_long Syntax:

    #include "umfpack.h"
    void *Numeric ;
    UF_long status ;
    double Mx, Ex, Info [UMFPACK_INFO] ;
    status = umfpack_dl_get_determinant (&Mx, &Ex, Numeric, Info) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    int status ;
    double Mx, Mz, Ex, Info [UMFPACK_INFO] ;
    status = umfpack_zi_get_determinant (&Mx, &Mz, &Ex, Numeric, Info) ;

complex int Syntax:

    #include "umfpack.h"
    void *Numeric ;
    UF_long status ;
    double *Mx, *Mz, *Ex, Info [UMFPACK_INFO] ;
    status = umfpack_zl_get_determinant (&Mx, &Mz, &Ex, Numeric, Info) ;

packed complex int Syntax:

    Same as above, except Mz is NULL.

Author: Contributed by David Bateman, Motorola, Paris

Purpose:

    Using the LU factors and the permutation vectors contained in the Numeric
    object, calculate the determinant of the matrix A.

    The value of the determinant can be returned in two forms, depending on
    whether Ex is NULL or not.  If Ex is NULL then the value of the determinant
    is returned on Mx and Mz for the real and imaginary parts.  However, to
    avoid over- or underflows, the determinant can be split into a mantissa
    and exponent, and the parts returned separately, in which case Ex is not
    NULL.  The actual determinant is then given by

      double det ;
      det = Mx * pow (10.0, Ex) ;

    for the double case, or

      double det [2] ;
      det [0] = Mx * pow (10.0, Ex) ;	    // real part
      det [1] = Mz * pow (10.0, Ex) ;	    // imaginary part

    for the complex case.  Information on if the determinant will or has
    over or under-flowed is given by Info [UMFPACK_STATUS].

    In the "packed complex" syntax, Mx [0] holds the real part and Mx [1]
    holds the imaginary part.  Mz is not used (it is NULL).

Returns:

    Returns UMFPACK_OK if sucessful.  Returns UMFPACK_ERROR_out_of_memory if
    insufficient memory is available for the n_row integer workspace that
    umfpack_*_get_determinant allocates to construct pivots from the
    permutation vectors.  Returns UMFPACK_ERROR_invalid_Numeric_object if the
    Numeric object provided as input is invalid.  Returns
    UMFPACK_WARNING_singular_matrix if the determinant is zero.  Returns
    UMFPACK_WARNING_determinant_underflow or
    UMFPACK_WARNING_determinant_overflow if the determinant has underflowed
    overflowed (for the case when Ex is NULL), or will overflow if Ex is not
    NULL and det is computed (see above) in the user program.

Arguments:

    double *Mx ;   Output argument (array of size 1, or size 2 if Mz is NULL)
    double *Mz ;   Output argument (optional)
    double *Ex ;   Output argument (optional)

        The determinant returned in mantissa/exponent form, as discussed above.
	If Mz is NULL, then both the original and imaginary parts will be
	returned in Mx. If Ex is NULL then the determinant is returned directly
	in Mx and Mz (or Mx [0] and Mx [1] if Mz is NULL), rather than in
	mantissa/exponent form.

    void *Numeric ;	Input argument, not modified.

	Numeric must point to a valid Numeric object, computed by
	umfpack_*_numeric.

    double Info [UMFPACK_INFO] ;	Output argument.

	Contains information about the calculation of the determinant. If a
	(double *) NULL pointer is passed, then no statistics are returned in
	Info (this is not an error condition).  The following statistics are
	computed in umfpack_*_determinant:

	Info [UMFPACK_STATUS]: status code.  This is also the return value,
	    whether or not Info is present.

	    UMFPACK_OK

	        The determinant was successfully found.

	    UMFPACK_ERROR_out_of_memory

		Insufficient memory to solve the linear system.

	    UMFPACK_ERROR_argument_missing

		Mx is missing (NULL).

	    UMFPACK_ERROR_invalid_Numeric_object

		The Numeric object is not valid.

	    UMFPACK_ERROR_invalid_system

		The matrix is rectangular.  Only square systems can be
		handled.

	    UMFPACK_WARNING_singluar_matrix

		The determinant is zero or NaN.  The matrix is singular.

	    UMFPACK_WARNING_determinant_underflow

	        When passing from mantissa/exponent form to the determinant
		an underflow has or will occur.  If the mantissa/exponent from
		of obtaining the determinant is used, the underflow will occur
		in the user program.  If the single argument method of
		obtaining the determinant is used, the underflow has already
		occurred.

	    UMFPACK_WARNING_determinant_overflow

	        When passing from mantissa/exponent form to the determinant
		an overflow has or will occur.  If the mantissa/exponent from
		of obtaining the determinant is used, the overflow will occur
		in the user program.  If the single argument method of
		obtaining the determinant is used, the overflow has already
		occurred.


*/
