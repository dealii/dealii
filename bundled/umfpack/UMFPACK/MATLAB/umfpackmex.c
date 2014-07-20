/* ========================================================================== */
/* === UMFPACK mexFunction ================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    MATLAB interface for umfpack.

    Factor or solve a sparse linear system, returning either the solution
    x to Ax=b or A'x'=b', or the factorization LU=P(R\A)Q or LU=PAQ.  A must be
    sparse, with nonzero dimensions, but it may be complex, singular, and/or
    rectangular.  b must be a dense n-by-1 vector (real or complex).
    L is unit lower triangular, U is upper triangular, and R is diagonal.
    P and Q are permutation matrices (permutations of an identity matrix).

    The matrix A is scaled, by default.  Each row i is divided by r (i), where
    r (i) is the sum of the absolute values of the entries in that row.  The
    scaled matrix has an infinity norm of 1.  The scale factors r (i) are
    returned in a diagonal sparse matrix.  If the factorization is:

	[L, U, P, Q, R] = umfpack (A) ;

    then the factorization is

	L*U = P * (R \ A) * Q

    This is safer than returning a matrix R such that L*U = P*R*A*Q, because
    it avoids the division by small entries.  If r(i) is subnormal, multiplying
    by 1/r(i) would result in an IEEE Infinity, but dividing by r(i) is safe.

    The factorization

	[L, U, P, Q] = umfpack (A) ;

    returns LU factors such that L*U = P*A*Q, with no scaling.

    See umfpack.m, umfpack_details.m, and umfpack.h for details.

    Note that this mexFunction accesses only the user-callable UMFPACK routines.
    Thus, is also provides another example of how user C code can access
    UMFPACK.

    If NO_TRANSPOSE_FORWARD_SLASH is not defined at compile time, then the
    forward slash (/) operator acts almost like x = b/A in MATLAB 6.1.  It is
    solved by factorizing the array transpose, and then x = (A.'\b.').' is
    solved.  This is the default behavior (for historical reasons), since
    factorizing A can behave perform much differently than factorizing its
    transpose.

    If NO_TRANSPOSE_FORWARD_SLASH is defined at compile time, then the forward
    slash operator does not act like x=b/A in MATLAB 6.1.  It is solved by
    factorizing A, and then solving via the transposed L and U matrices.
    The solution is still x = (A.'\b.').', except that A is factorized instead
    of A.'.

    Modified for v4.3.1, Jan 10, 2005: default has been changed to
    NO_TRANSPOSE_FORWARD_SLASH, to test iterative refinement for b/A.
    v4.4: added method for computing the determinant.

    v5.1: port to 64-bit MATLAB
*/

#define NO_TRANSPOSE_FORWARD_SLASH  /* default has changed for v4.3.1 */

#include "UFconfig.h"
#include "umfpack.h"
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <math.h>
#include <float.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define STRING_MATCH(s1,s2) (strcmp ((s1), (s2)) == 0)
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif

/* ========================================================================== */
/* === error ================================================================ */
/* ========================================================================== */

/* Return an error message */

static void error
(
    char *s,
    UF_long A_is_complex,
    int nargout,
    mxArray *pargout [ ],
    double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO],
    UF_long status,
    UF_long do_info
)
{
    UF_long i ;
    double *Out_Info ;
    if (A_is_complex)
    {
	umfpack_zl_report_status (Control, status) ;
	umfpack_zl_report_info (Control, Info) ;
    }
    else
    {
	umfpack_dl_report_status (Control, status) ;
	umfpack_dl_report_info (Control, Info) ;
    }
    if (do_info > 0)
    {
	/* return Info */
	pargout [do_info] = mxCreateDoubleMatrix (1, UMFPACK_INFO, mxREAL) ;
	Out_Info = mxGetPr (pargout [do_info]) ;
	for (i = 0 ; i < UMFPACK_INFO ; i++)
	{
	    Out_Info [i] = Info [i] ;
	}
    }
    mexErrMsgTxt (s) ;
}


/* ========================================================================== */
/* === UMFPACK ============================================================== */
/* ========================================================================== */

void mexFunction
(
    int nargout,		/* number of outputs */
    mxArray *pargout [ ],	/* output arguments */
    int nargin,			/* number of inputs */
    const mxArray *pargin [ ]	/* input arguments */
)
{

    /* ---------------------------------------------------------------------- */
    /* local variables */
    /* ---------------------------------------------------------------------- */

    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL], dx, dz, dexp ;
    double *Lx, *Lz, *Ux, *Uz, *Ax, *Az, *Bx, *Bz, *Xx, *Xz, *User_Control,
	*p, *q, *Out_Info, *p1, *p2, *p3, *p4, *Ltx, *Ltz, *Rs, *Px, *Qx ;
    void *Symbolic, *Numeric ;
    UF_long *Lp, *Li, *Up, *Ui, *Ap, *Ai, *P, *Q, do_solve, lnz, unz, nn, i,
	transpose, size, do_info, do_numeric, *Front_npivcol, op, k, *Rp, *Ri,
	*Front_parent, *Chain_start, *Chain_maxrows, *Chain_maxcols, nz, status,
	nfronts, nchains, *Ltp, *Ltj, *Qinit, print_level, status2, no_scale,
	*Front_1strow, *Front_leftmostdesc, n_row, n_col, n_inner, sys,
	ignore1, ignore2, ignore3, A_is_complex, B_is_complex, X_is_complex,
	*Pp, *Pi, *Qp, *Qi, do_recip, do_det ;
    mxArray *Amatrix, *Bmatrix, *User_Control_matrix, *User_Qinit ;
    char *operator, *operation ;
    mxComplexity Atype, Xtype ;
    char warning [200] ;

#ifndef NO_TRANSPOSE_FORWARD_SLASH
    UF_long *Cp, *Ci ;
    double *Cx, *Cz ;
#endif

    /* ---------------------------------------------------------------------- */
    /* define the memory manager and printf functions for UMFPACK and AMD */ 
    /* ---------------------------------------------------------------------- */

    /* with these settings, the UMFPACK mexFunction can use ../Lib/libumfpack.a
     * and ../Lib/libamd.a, instead compiling UMFPACK and AMD specifically for
     * the MATLAB mexFunction. */
    amd_malloc = mxMalloc ;
    amd_free = mxFree ;
    amd_calloc = mxCalloc ;
    amd_realloc = mxRealloc ;
    amd_printf = mexPrintf ;

    /* The default values for these function pointers are fine.
    umfpack_hypot = umf_hypot ;
    umfpack_divcomplex = umf_divcomplex ;
    */

    /* ---------------------------------------------------------------------- */
    /* get inputs A, b, and the operation to perform */
    /* ---------------------------------------------------------------------- */

    User_Control_matrix = (mxArray *) NULL ;
    User_Qinit = (mxArray *) NULL ;

    do_info = 0 ;
    do_solve = FALSE ;
    do_numeric = TRUE ;
    transpose = FALSE ;
    no_scale = FALSE ;
    do_det = FALSE ;

    /* find the operator */
    op = 0 ;
    for (i = 0 ; i < nargin ; i++)
    {
	if (mxIsChar (pargin [i]))
	{
	    op = i ;
	    break ;
	}
    }

    if (op > 0)
    {
	operator = mxArrayToString (pargin [op]) ;

	if (STRING_MATCH (operator, "\\"))
	{

	    /* -------------------------------------------------------------- */
	    /* matrix left divide, x = A\b */
	    /* -------------------------------------------------------------- */

	    /*
		[x, Info] = umfpack (A, '\', b) ;
		[x, Info] = umfpack (A, '\', b, Control) ;
		[x, Info] = umfpack (A, Qinit, '\', b, Control) ;
		[x, Info] = umfpack (A, Qinit, '\', b) ;
	    */

	    operation = "x = A\\b" ;
	    do_solve = TRUE ;
	    Amatrix = (mxArray *) pargin [0] ;
	    Bmatrix = (mxArray *) pargin [op+1] ;

	    if (nargout == 2)
	    {
		do_info = 1 ;
	    }
	    if (op == 2)
	    {
		User_Qinit = (mxArray *) pargin [1] ;
	    }
	    if ((op == 1 && nargin == 4) || (op == 2 && nargin == 5))
	    {
		User_Control_matrix = (mxArray *) pargin [nargin-1] ;
	    }
	    if (nargin < 3 || nargin > 5 || nargout > 2)
	    {
		mexErrMsgTxt ("wrong number of arguments") ;
	    }

	}
	else if (STRING_MATCH (operator, "/"))
	{

	    /* -------------------------------------------------------------- */
	    /* matrix right divide, x = b/A */
	    /* -------------------------------------------------------------- */

	    /*
		[x, Info] = umfpack (b, '/', A) ;
		[x, Info] = umfpack (b, '/', A, Control) ;
		[x, Info] = umfpack (b, '/', A, Qinit) ;
		[x, Info] = umfpack (b, '/', A, Qinit, Control) ;
	    */

	    operation = "x = b/A" ;
	    do_solve = TRUE ;
	    transpose = TRUE ;
	    Amatrix = (mxArray *) pargin [2] ;
	    Bmatrix = (mxArray *) pargin [0] ;

	    if (nargout == 2)
	    {
		do_info = 1 ;
	    }
	    if (nargin == 5)
	    {
		User_Qinit = (mxArray *) pargin [3] ;
		User_Control_matrix = (mxArray *) pargin [4] ;
	    }
	    else if (nargin == 4)
	    {
		/* Control is k-by-1 where k > 1, Qinit is 1-by-n */
		if (mxGetM (pargin [3]) == 1)
		{
		    User_Qinit = (mxArray *) pargin [3] ;
		}
		else
		{
		    User_Control_matrix = (mxArray *) pargin [3] ;
		}
	    }
	    else if (nargin < 3 || nargin > 5 || nargout > 2)
	    {
		mexErrMsgTxt ("wrong number of arguments") ;
	    }

	}
	else if (STRING_MATCH (operator, "symbolic"))
	{

	    /* -------------------------------------------------------------- */
	    /* symbolic factorization only */
	    /* -------------------------------------------------------------- */

	    /*
	    [P Q Fr Ch Info] = umfpack (A, 'symbolic') ;
	    [P Q Fr Ch Info] = umfpack (A, 'symbolic', Control) ;
	    [P Q Fr Ch Info] = umfpack (A, Qinit, 'symbolic') ;
	    [P Q Fr Ch Info] = umfpack (A, Qinit, 'symbolic', Control) ;
	    */

	    operation = "symbolic factorization" ;
	    do_numeric = FALSE ;
	    Amatrix = (mxArray *) pargin [0] ;

	    if (nargout == 5)
	    {
		do_info = 4 ;
	    }
	    if (op == 2)
	    {
		User_Qinit = (mxArray *) pargin [1] ;
	    }
	    if ((op == 1 && nargin == 3) || (op == 2 && nargin == 4))
	    {
		User_Control_matrix = (mxArray *) pargin [nargin-1] ;
	    }
	    if (nargin < 2 || nargin > 4 || nargout > 5 || nargout < 4)
	    {
		mexErrMsgTxt ("wrong number of arguments") ;
	    }

	}
	else if (STRING_MATCH (operator, "det"))
	{

	    /* -------------------------------------------------------------- */
	    /* compute the determinant */
	    /* -------------------------------------------------------------- */

	    /*
	     * [det] = umfpack (A, 'det') ;
	     * [dmantissa dexp] = umfpack (A, 'det') ;
	     */

	    operation = "determinant" ;
	    do_det = TRUE ;
	    Amatrix = (mxArray *) pargin [0] ;
	    if (nargin > 2 || nargout > 2)
	    {
		mexErrMsgTxt ("wrong number of arguments") ;
	    }

	}
	else
	{
	    mexErrMsgTxt ("operator must be '/', '\\', or 'symbolic'") ;
	}
	mxFree (operator) ;

    }
    else if (nargin > 0)
    {

	/* ------------------------------------------------------------------ */
	/* LU factorization */
	/* ------------------------------------------------------------------ */

	/*
	    with scaling:
	    [L, U, P, Q, R, Info] = umfpack (A) ;
	    [L, U, P, Q, R, Info] = umfpack (A, Qinit) ;

	    scaling determined by Control settings:
	    [L, U, P, Q, R, Info] = umfpack (A, Control) ;
	    [L, U, P, Q, R, Info] = umfpack (A, Qinit, Control) ;

	    with no scaling:
	    [L, U, P, Q] = umfpack (A) ;
	    [L, U, P, Q] = umfpack (A, Control) ;
	    [L, U, P, Q] = umfpack (A, Qinit) ;
	    [L, U, P, Q] = umfpack (A, Qinit, Control) ;
	*/

	operation = "numeric factorization" ;
	Amatrix = (mxArray *) pargin [0] ;

	no_scale = nargout <= 4 ;

	if (nargout == 6)
	{
	    do_info = 5 ;
	}
	if (nargin == 3)
	{
	    User_Qinit = (mxArray *) pargin [1] ;
	    User_Control_matrix = (mxArray *) pargin [2] ;
	}
	else if (nargin == 2)
	{
	    /* Control is k-by-1 where k > 1, Qinit is 1-by-n */
	    if (mxGetM (pargin [1]) == 1)
	    {
		User_Qinit = (mxArray *) pargin [1] ;
	    }
	    else
	    {
		User_Control_matrix = (mxArray *) pargin [1] ;
	    }
	}
	else if (nargin > 3 || nargout > 6 || nargout < 4)
	{
	    mexErrMsgTxt ("wrong number of arguments") ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* return default control settings */
	/* ------------------------------------------------------------------ */

	/*
	    Control = umfpack ;
	    umfpack ;
	*/

	if (nargout > 1)
	{
	    mexErrMsgTxt ("wrong number of arguments") ;
	}

	pargout [0] = mxCreateDoubleMatrix (UMFPACK_CONTROL, 1, mxREAL) ;
	User_Control = mxGetPr (pargout [0]) ;
	umfpack_dl_defaults (User_Control) ;

	return ;
    }

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (mxGetNumberOfDimensions (Amatrix) != 2)
    {
	mexErrMsgTxt ("input matrix A must be 2-dimensional") ;
    }
    n_row = mxGetM (Amatrix) ;
    n_col = mxGetN (Amatrix) ;
    nn = MAX (n_row, n_col) ;
    n_inner = MIN (n_row, n_col) ;
    if (do_solve && n_row != n_col)
    {
	mexErrMsgTxt ("input matrix A must square for '\\' or '/'") ;
    }
    if (!mxIsSparse (Amatrix))
    {
	mexErrMsgTxt ("input matrix A must be sparse") ;
    }
    if (n_row == 0 || n_col == 0)
    {
	mexErrMsgTxt ("input matrix A cannot have zero rows or zero columns") ;
    }

    /* The real/complex status of A determines which version to use, */
    /* (umfpack_dl_* or umfpack_zl_*). */
    A_is_complex = mxIsComplex (Amatrix) ;
    Atype = A_is_complex ? mxCOMPLEX : mxREAL ;
    Ap = (UF_long *) mxGetJc (Amatrix) ;
    Ai = (UF_long *) mxGetIr (Amatrix) ;
    Ax = mxGetPr (Amatrix) ;
    Az = mxGetPi (Amatrix) ;

    if (do_solve)
    {

	if (n_row != n_col)
	{
	    mexErrMsgTxt ("A must be square for \\ or /") ;
	}
	if (transpose)
	{
	    if (mxGetM (Bmatrix) != 1 || mxGetN (Bmatrix) != nn)
	    {
		mexErrMsgTxt ("b has the wrong dimensions") ;
	    }
	}
	else
	{
	    if (mxGetM (Bmatrix) != nn || mxGetN (Bmatrix) != 1)
	    {
		mexErrMsgTxt ("b has the wrong dimensions") ;
	    }
	}
	if (mxGetNumberOfDimensions (Bmatrix) != 2)
	{
	    mexErrMsgTxt ("input matrix b must be 2-dimensional") ;
	}
	if (mxIsSparse (Bmatrix))
	{
	    mexErrMsgTxt ("input matrix b cannot be sparse") ;
	}
	if (mxGetClassID (Bmatrix) != mxDOUBLE_CLASS)
	{
	    mexErrMsgTxt ("input matrix b must double precision matrix") ;
	}

	B_is_complex = mxIsComplex (Bmatrix) ;
	Bx = mxGetPr (Bmatrix) ;
	Bz = mxGetPi (Bmatrix) ;

	X_is_complex = A_is_complex || B_is_complex ;
	Xtype = X_is_complex ? mxCOMPLEX : mxREAL ;
    }

    /* ---------------------------------------------------------------------- */
    /* set the Control parameters */
    /* ---------------------------------------------------------------------- */

    if (A_is_complex)
    {
	umfpack_zl_defaults (Control) ;
    }
    else
    {
	umfpack_dl_defaults (Control) ;
    }
    if (User_Control_matrix)
    {
	if (mxGetClassID (User_Control_matrix) != mxDOUBLE_CLASS ||
	    mxIsSparse (User_Control_matrix))
	{
	    mexErrMsgTxt ("Control must be a dense real matrix") ;
	}
	size = UMFPACK_CONTROL ;
	size = MIN (size, mxGetNumberOfElements (User_Control_matrix)) ;
	User_Control = mxGetPr (User_Control_matrix) ;
	for (i = 0 ; i < size ; i++)
	{
	    Control [i] = User_Control [i] ;
	}
    }

    if (no_scale)
    {
	/* turn off scaling for [L, U, P, Q] = umfpack (A) ;
	 * ignoring the input value of Control (24) for the usage
	 * [L, U, P, Q] = umfpack (A, Control) ; */
	Control [UMFPACK_SCALE] = UMFPACK_SCALE_NONE ;
    }

    if (mxIsNaN (Control [UMFPACK_PRL]))
    {
	print_level = UMFPACK_DEFAULT_PRL ;
    }
    else
    {
	print_level = (UF_long) Control [UMFPACK_PRL] ;
    }

    Control [UMFPACK_PRL] = print_level ;

    /* ---------------------------------------------------------------------- */
    /* get Qinit, if present */
    /* ---------------------------------------------------------------------- */

    if (User_Qinit)
    {
	if (mxGetM (User_Qinit) != 1 || mxGetN (User_Qinit) != n_col)
	{
	    mexErrMsgTxt ("Qinit must be 1-by-n_col") ;
	}
	if (mxGetNumberOfDimensions (User_Qinit) != 2)
	{
	    mexErrMsgTxt ("input Qinit must be 2-dimensional") ;
	}
	if (mxIsComplex (User_Qinit))
	{
	    mexErrMsgTxt ("input Qinit must not be complex") ;
	}
	if (mxGetClassID (User_Qinit) != mxDOUBLE_CLASS)
	{
	    mexErrMsgTxt ("input Qinit must be a double matrix") ;
	}
	if (mxIsSparse (User_Qinit))
	{
	    mexErrMsgTxt ("input Qinit must be dense") ;
	}
	Qinit = (UF_long *) mxMalloc (n_col * sizeof (UF_long)) ;
	p = mxGetPr (User_Qinit) ;
	for (k = 0 ; k < n_col ; k++)
	{
	    /* convert from 1-based to 0-based indexing */
	    Qinit [k] = ((UF_long) (p [k])) - 1 ;
	}

    }
    else
    {
	/* umfpack_*_qsymbolic will call colamd to get Qinit. This is the */
	/* same as calling umfpack_*_symbolic with Qinit set to NULL*/
	Qinit = (UF_long *) NULL ;
    }

    /* ---------------------------------------------------------------------- */
    /* report the inputs A and Qinit */
    /* ---------------------------------------------------------------------- */

    if (print_level >= 2)
    {
	/* print the operation */
	mexPrintf ("\numfpack: %s\n", operation) ;
    }

    if (A_is_complex)
    {
	umfpack_zl_report_control (Control) ;
	if (print_level >= 3) mexPrintf ("\nA: ") ;
	(void) umfpack_zl_report_matrix (n_row, n_col, Ap, Ai, Ax, Az,
	    1, Control) ;
	if (Qinit)
	{
	    if (print_level >= 3) mexPrintf ("\nQinit: ") ;
	    (void) umfpack_zl_report_perm (n_col, Qinit, Control) ;
	}
    }
    else
    {
	umfpack_dl_report_control (Control) ;
	if (print_level >= 3) mexPrintf ("\nA: ") ;
	(void) umfpack_dl_report_matrix (n_row, n_col, Ap, Ai, Ax,
	    1, Control) ;
	if (Qinit)
	{
	    if (print_level >= 3) mexPrintf ("\nQinit: ") ;
	    (void) umfpack_dl_report_perm (n_col, Qinit, Control) ;
	}
    }

#ifndef NO_TRANSPOSE_FORWARD_SLASH
    /* ---------------------------------------------------------------------- */
    /* create the array transpose for x = b/A */
    /* ---------------------------------------------------------------------- */

    if (transpose)
    {
	/* note that in this case A will be square (nn = n_row = n_col) */
	/* x = (A.'\b.').' will be computed */

	/* make sure Ci and Cx exist, avoid malloc of zero-sized arrays. */
	nz = MAX (Ap [nn], 1) ;

	Cp = (UF_long *) mxMalloc ((nn+1) * sizeof (UF_long)) ;
	Ci = (UF_long *) mxMalloc (nz * sizeof (UF_long)) ;
	Cx = (double *) mxMalloc (nz * sizeof (double)) ;
	if (A_is_complex)
	{
	    Cz = (double *) mxMalloc (nz * sizeof (double)) ;
	    status = umfpack_zl_transpose (nn, nn, Ap, Ai, Ax, Az,
	        (UF_long *) NULL, (UF_long *) NULL, Cp, Ci, Cx, Cz, FALSE) ;
	}
	else
	{
	    status = umfpack_dl_transpose (nn, nn, Ap, Ai, Ax,
	        (UF_long *) NULL, (UF_long *) NULL, Cp, Ci, Cx) ;
	}

	if (status != UMFPACK_OK)
	{
	    error ("transpose of A failed", A_is_complex, nargout, pargout,
		Control, Info, status, do_info);
	    return ;
	}

	/* modify pointers so that C will be factorized and solved, not A */
	Ap = Cp ;
	Ai = Ci ;
	Ax = Cx ;
	if (A_is_complex)
	{
	    Az = Cz ;
	}
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* perform the symbolic factorization */
    /* ---------------------------------------------------------------------- */

    if (A_is_complex)
    {
	status = umfpack_zl_qsymbolic (n_row, n_col, Ap, Ai, Ax, Az,
	    Qinit, &Symbolic, Control, Info) ;
    }
    else
    {
	status = umfpack_dl_qsymbolic (n_row, n_col, Ap, Ai, Ax,
	    Qinit, &Symbolic, Control, Info) ;
    }

    if (Qinit)
    {
	mxFree (Qinit) ;
    }

    if (status < 0)
    {
	error ("symbolic factorization failed", A_is_complex, nargout, pargout,
	    Control, Info, status, do_info) ;
	return ;
    }

    /* ---------------------------------------------------------------------- */
    /* report the Symbolic object */
    /* ---------------------------------------------------------------------- */

    if (A_is_complex)
    {
	(void) umfpack_zl_report_symbolic (Symbolic, Control) ;
    }
    else
    {
	(void) umfpack_dl_report_symbolic (Symbolic, Control) ;
    }

    /* ---------------------------------------------------------------------- */
    /* perform numeric factorization, or just return symbolic factorization */
    /* ---------------------------------------------------------------------- */

    if (do_numeric)
    {

	/* ------------------------------------------------------------------ */
	/* perform the numeric factorization */
	/* ------------------------------------------------------------------ */

	if (A_is_complex)
	{
	    status = umfpack_zl_numeric (Ap, Ai, Ax, Az, Symbolic, &Numeric,
		Control, Info) ;
	}
	else
	{
	    status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric,
		Control, Info) ;
	}

	/* ------------------------------------------------------------------ */
	/* free the symbolic factorization */
	/* ------------------------------------------------------------------ */

	if (A_is_complex)
	{
	    umfpack_zl_free_symbolic (&Symbolic) ;
	}
	else
	{
	    umfpack_dl_free_symbolic (&Symbolic) ;
	}

	/* ------------------------------------------------------------------ */
	/* report the Numeric object */
	/* ------------------------------------------------------------------ */

	if (status < 0)
	{
	    error ("numeric factorization failed", A_is_complex, nargout,
		pargout, Control, Info, status, do_info);
	    return ;
	}

	if (A_is_complex)
	{
	    (void) umfpack_zl_report_numeric (Numeric, Control) ;
	}
	else
	{
	    (void) umfpack_dl_report_numeric (Numeric, Control) ;
	}

	/* ------------------------------------------------------------------ */
	/* return the solution, determinant, or the factorization */
	/* ------------------------------------------------------------------ */

	if (do_solve)
	{
	    /* -------------------------------------------------------------- */
	    /* solve Ax=b or A'x'=b', and return just the solution x */
	    /* -------------------------------------------------------------- */

#ifndef NO_TRANSPOSE_FORWARD_SLASH
	    if (transpose)
	    {
		/* A.'x.'=b.' gives the same x=b/A as solving A'x'=b' */
		/* since C=A.' was factorized, solve with sys = UMFPACK_A */
		/* since x and b are vectors, x.' and b.' are implicit */
		pargout [0] = mxCreateDoubleMatrix (1, nn, Xtype) ;
	    }
	    else
	    {
		pargout [0] = mxCreateDoubleMatrix (nn, 1, Xtype) ;
	    }
	    sys = UMFPACK_A ;
#else
	    if (transpose)
	    {
		/* If A is real, A'x=b is the same as A.'x=b. */
		/* x and b are vectors, so x and b are the same as x' and b'. */
		/* If A is complex, then A.'x.'=b.' gives the same solution x */
		/* as the complex conjugate transpose.  If we used the A'x=b */
		/* option in umfpack_*_solve, we would have to form b' on */
		/* input and x' on output (negating the imaginary part). */
		/* We can save this work by just using the A.'x=b option in */
		/* umfpack_*_solve.  Then, forming x.' and b.' is implicit, */
		/* since x and b are just vectors anyway. */
		/* In both cases, the system to solve is A.'x=b */
		pargout [0] = mxCreateDoubleMatrix (1, nn, Xtype) ;
		sys = UMFPACK_Aat ;
	    }
	    else
	    {
		pargout [0] = mxCreateDoubleMatrix (nn, 1, Xtype) ;
		sys = UMFPACK_A ;
	    }
#endif

	    /* -------------------------------------------------------------- */
	    /* print the right-hand-side, B */
	    /* -------------------------------------------------------------- */

	    if (print_level >= 3) mexPrintf ("\nright-hand side, b: ") ;
	    if (B_is_complex)
	    {
		(void) umfpack_zl_report_vector (nn, Bx, Bz, Control) ;
	    }
	    else
	    {
		(void) umfpack_dl_report_vector (nn, Bx, Control) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* solve the system */
	    /* -------------------------------------------------------------- */

	    Xx = mxGetPr (pargout [0]) ;
	    Xz = mxGetPi (pargout [0]) ;
	    status2 = UMFPACK_OK ;

	    if (A_is_complex)
	    {
		if (!B_is_complex)
		{
		    /* umfpack_zl_solve expects a complex B */
		    Bz = (double *) mxCalloc (nn, sizeof (double)) ;
		}
		status = umfpack_zl_solve (sys, Ap, Ai, Ax, Az, Xx, Xz, Bx, Bz,
		    Numeric, Control, Info) ;
		if (!B_is_complex)
		{
		    mxFree (Bz) ;
		}
	    }
	    else
	    {
		if (B_is_complex)
		{
		    /* Ax=b when b is complex and A is sparse can be split */
		    /* into two systems, A*xr=br and A*xi=bi, where r denotes */
		    /* the real part and i the imaginary part of x and b. */
		    status2 = umfpack_dl_solve (sys, Ap, Ai, Ax, Xz, Bz,
		    Numeric, Control, Info) ;
		}
		status = umfpack_dl_solve (sys, Ap, Ai, Ax, Xx, Bx,
		    Numeric, Control, Info) ;
	    }

#ifndef NO_TRANSPOSE_FORWARD_SLASH
	    /* -------------------------------------------------------------- */
	    /* free the transposed matrix C */
	    /* -------------------------------------------------------------- */

	    if (transpose)
	    {
	        mxFree (Cp) ;
	        mxFree (Ci) ;
	        mxFree (Cx) ;
	        if (A_is_complex)
	        {
	            mxFree (Cz) ;
	        }
	    }
#endif

	    /* -------------------------------------------------------------- */
	    /* free the Numeric object */
	    /* -------------------------------------------------------------- */

	    if (A_is_complex)
	    {
		umfpack_zl_free_numeric (&Numeric) ;
	    }
	    else
	    {
		umfpack_dl_free_numeric (&Numeric) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* check error status */
	    /* -------------------------------------------------------------- */

	    if (status < 0 || status2 < 0)
	    {
		mxDestroyArray (pargout [0]) ;
		error ("solve failed", A_is_complex, nargout, pargout, Control,
			Info, status, do_info) ;
		return ;
	    }

	    /* -------------------------------------------------------------- */
	    /* print the solution, X */
	    /* -------------------------------------------------------------- */

	    if (print_level >= 3) mexPrintf ("\nsolution, x: ") ;
	    if (X_is_complex)
	    {
		(void) umfpack_zl_report_vector (nn, Xx, Xz, Control) ;
	    }
	    else
	    {
		(void) umfpack_dl_report_vector (nn, Xx, Control) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* warn about singular or near-singular matrices */
	    /* -------------------------------------------------------------- */

	    /* no warning is given if Control (1) is zero */

	    if (Control [UMFPACK_PRL] >= 1)
	    {
		if (status == UMFPACK_WARNING_singular_matrix)
		{
		    sprintf (warning, "matrix is singular\n"
			"Try increasing Control (%d) and Control (%d).\n"
			"(Suppress this warning with Control (%d) = 0.)\n",
			1+UMFPACK_PIVOT_TOLERANCE,
			1+UMFPACK_SYM_PIVOT_TOLERANCE,
			1+UMFPACK_PRL) ;
		    mexWarnMsgTxt (warning) ;
		}
		else if (Info [UMFPACK_RCOND] < DBL_EPSILON)
		{
		    sprintf (warning, "matrix is nearly singular, rcond = %g\n"
			"Try increasing Control (%d) and Control (%d).\n"
			"(Suppress this warning with Control (%d) = 0.)\n",
			Info [UMFPACK_RCOND],
			1+UMFPACK_PIVOT_TOLERANCE,
			1+UMFPACK_SYM_PIVOT_TOLERANCE,
			1+UMFPACK_PRL) ;
		    mexWarnMsgTxt (warning) ;
		}
	    }

	}
	else if (do_det)
	{

	    /* -------------------------------------------------------------- */
	    /* get the determinant */
	    /* -------------------------------------------------------------- */

	    if (nargout == 2)
	    {
		/* [det dexp] = umfpack (A, 'det') ;
		 * return determinant in the form det * 10^dexp */
		p = &dexp ;
	    }
	    else
	    {
		/* [det] = umfpack (A, 'det') ;
		 * return determinant as a single scalar (overflow or
		 * underflow is much more likely) */
		p = (double *) NULL ;
	    }
	    if (A_is_complex)
	    {
		status = umfpack_zl_get_determinant (&dx, &dz, p,
			Numeric, Info) ;
		umfpack_zl_free_numeric (&Numeric) ;
	    }
	    else
	    {
		status = umfpack_dl_get_determinant (&dx, p,
			Numeric, Info) ;
		umfpack_dl_free_numeric (&Numeric) ;
		dz = 0 ;
	    }
	    if (status < 0)
	    {
		error ("extracting LU factors failed", A_is_complex, nargout,
		    pargout, Control, Info, status, do_info) ;
	    }
	    if (A_is_complex)
	    {
		pargout [0] = mxCreateDoubleMatrix (1, 1, mxCOMPLEX) ;
		p = mxGetPr (pargout [0]) ;
		*p = dx ;
		p = mxGetPi (pargout [0]) ;
		*p = dz ;
	    }
	    else
	    {
		pargout [0] = mxCreateDoubleMatrix (1, 1, mxREAL) ;
		p = mxGetPr (pargout [0]) ;
		*p = dx ;
	    }
	    if (nargout == 2)
	    {
		pargout [1] = mxCreateDoubleMatrix (1, 1, mxREAL) ;
		p = mxGetPr (pargout [1]) ;
		*p = dexp ;
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* get L, U, P, Q, and r */
	    /* -------------------------------------------------------------- */

	    if (A_is_complex)
	    {
	        status = umfpack_zl_get_lunz (&lnz, &unz, &ignore1, &ignore2,
		    &ignore3, Numeric) ;
	    }
	    else
	    {
	        status = umfpack_dl_get_lunz (&lnz, &unz, &ignore1, &ignore2,
		    &ignore3, Numeric) ;
	    }

	    if (status < 0)
	    {
		if (A_is_complex)
		{
		    umfpack_zl_free_numeric (&Numeric) ;
		}
		else
		{
		    umfpack_dl_free_numeric (&Numeric) ;
		}
		error ("extracting LU factors failed", A_is_complex, nargout,
		    pargout, Control, Info, status, do_info) ;
		return ;
	    }

	    /* avoid malloc of zero-sized arrays */
	    lnz = MAX (lnz, 1) ;
	    unz = MAX (unz, 1) ;

	    /* get temporary space, for the *** ROW *** form of L */
	    Ltp = (UF_long *) mxMalloc ((n_row+1) * sizeof (UF_long)) ;
	    Ltj = (UF_long *) mxMalloc (lnz * sizeof (UF_long)) ;
	    Ltx = (double *) mxMalloc (lnz * sizeof (double)) ;
	    if (A_is_complex)
	    {
	        Ltz = (double *) mxMalloc (lnz * sizeof (double)) ;
	    }
	    else
	    {
	        Ltz = (double *) NULL ;
	    }

	    /* create permanent copy of the output matrix U */
	    pargout [1] = mxCreateSparse (n_inner, n_col, unz, Atype) ;
	    Up = (UF_long *) mxGetJc (pargout [1]) ;
	    Ui = (UF_long *) mxGetIr (pargout [1]) ;
	    Ux = mxGetPr (pargout [1]) ;
	    Uz = mxGetPi (pargout [1]) ;

	    /* temporary space for the integer permutation vectors */
	    P = (UF_long *) mxMalloc (n_row * sizeof (UF_long)) ;
	    Q = (UF_long *) mxMalloc (n_col * sizeof (UF_long)) ;

	    /* get scale factors, if requested */
	    status2 = UMFPACK_OK ;
	    if (!no_scale)
	    {
		/* create a diagonal sparse matrix for the scale factors */
		pargout [4] = mxCreateSparse (n_row, n_row, n_row, mxREAL) ;
		Rp = (UF_long *) mxGetJc (pargout [4]) ;
		Ri = (UF_long *) mxGetIr (pargout [4]) ;
		for (i = 0 ; i < n_row ; i++)
		{
		    Rp [i] = i ;
		    Ri [i] = i ;
		}
		Rp [n_row] = n_row ;
		Rs = mxGetPr (pargout [4]) ;
	    }
	    else
	    {
		Rs = (double *) NULL ;
	    }

	    /* get Lt, U, P, Q, and Rs from the numeric object */
	    if (A_is_complex)
	    {
		status = umfpack_zl_get_numeric (Ltp, Ltj, Ltx, Ltz, Up, Ui, Ux,
		    Uz, P, Q, (double *) NULL, (double *) NULL,
		    &do_recip, Rs, Numeric) ;
		umfpack_zl_free_numeric (&Numeric) ;
	    }
	    else
	    {
		status = umfpack_dl_get_numeric (Ltp, Ltj, Ltx, Up, Ui,
		    Ux, P, Q, (double *) NULL,
		    &do_recip, Rs, Numeric) ;
		umfpack_dl_free_numeric (&Numeric) ;
	    }

	    /* for the mexFunction, -DNRECIPROCAL must be set,
	     * so do_recip must be FALSE */

	    if (status < 0 || status2 < 0 || do_recip)
	    {
		mxFree (Ltp) ;
		mxFree (Ltj) ;
		mxFree (Ltx) ;
		if (Ltz) mxFree (Ltz) ;
		mxFree (P) ;
		mxFree (Q) ;
		mxDestroyArray (pargout [1]) ;
		error ("extracting LU factors failed", A_is_complex, nargout,
		    pargout, Control, Info, status, do_info) ;
		return ;
	    }

	    /* create sparse permutation matrix for P */
	    pargout [2] = mxCreateSparse (n_row, n_row, n_row, mxREAL) ;
	    Pp = (UF_long *) mxGetJc (pargout [2]) ;
	    Pi = (UF_long *) mxGetIr (pargout [2]) ;
	    Px = mxGetPr (pargout [2]) ;
	    for (k = 0 ; k < n_row ; k++)
	    {
		Pp [k] = k ;
		Px [k] = 1 ;
		Pi [P [k]] = k ;
	    }
	    Pp [n_row] = n_row ;

	    /* create sparse permutation matrix for Q */
	    pargout [3] = mxCreateSparse (n_col, n_col, n_col, mxREAL) ;
	    Qp = (UF_long *) mxGetJc (pargout [3]) ;
	    Qi = (UF_long *) mxGetIr (pargout [3]) ;
	    Qx = mxGetPr (pargout [3]) ;
	    for (k = 0 ; k < n_col ; k++)
	    {
		Qp [k] = k ;
		Qx [k] = 1 ;
		Qi [k] = Q [k] ;
	    }
	    Qp [n_col] = n_col ;

	    /* permanent copy of L */
	    pargout [0] = mxCreateSparse (n_row, n_inner, lnz, Atype) ;
	    Lp = (UF_long *) mxGetJc (pargout [0]) ;
	    Li = (UF_long *) mxGetIr (pargout [0]) ;
	    Lx = mxGetPr (pargout [0]) ;
	    Lz = mxGetPi (pargout [0]) ;

	    /* convert L from row form to column form */
	    if (A_is_complex)
	    {
		/* non-conjugate array transpose */
	        status = umfpack_zl_transpose (n_inner, n_row, Ltp, Ltj, Ltx,
		    Ltz, (UF_long *) NULL, (UF_long *) NULL, Lp, Li, Lx, Lz,
		    FALSE) ;
	    }
	    else
	    {
	        status = umfpack_dl_transpose (n_inner, n_row, Ltp, Ltj, Ltx,
		    (UF_long *) NULL, (UF_long *) NULL, Lp, Li, Lx) ;
	    }

	    mxFree (Ltp) ;
	    mxFree (Ltj) ;
	    mxFree (Ltx) ;
	    if (Ltz) mxFree (Ltz) ;

	    if (status < 0)
	    {
		mxFree (P) ;
		mxFree (Q) ;
		mxDestroyArray (pargout [0]) ;
		mxDestroyArray (pargout [1]) ;
		mxDestroyArray (pargout [2]) ;
		mxDestroyArray (pargout [3]) ;
		error ("constructing L failed", A_is_complex, nargout, pargout,
		    Control, Info, status, do_info) ;
		return ;
	    }

	    /* -------------------------------------------------------------- */
	    /* print L, U, P, and Q */
	    /* -------------------------------------------------------------- */

	    if (A_is_complex)
	    {
		if (print_level >= 3) mexPrintf ("\nL: ") ;
	        (void) umfpack_zl_report_matrix (n_row, n_inner, Lp, Li,
		    Lx, Lz, 1, Control) ;
		if (print_level >= 3) mexPrintf ("\nU: ") ;
	        (void) umfpack_zl_report_matrix (n_inner, n_col,  Up, Ui,
		    Ux, Uz, 1, Control) ;
		if (print_level >= 3) mexPrintf ("\nP: ") ;
	        (void) umfpack_zl_report_perm (n_row, P, Control) ;
		if (print_level >= 3) mexPrintf ("\nQ: ") ;
	        (void) umfpack_zl_report_perm (n_col, Q, Control) ;
	    }
	    else
	    {
		if (print_level >= 3) mexPrintf ("\nL: ") ;
	        (void) umfpack_dl_report_matrix (n_row, n_inner, Lp, Li,
		    Lx, 1, Control) ;
		if (print_level >= 3) mexPrintf ("\nU: ") ;
	        (void) umfpack_dl_report_matrix (n_inner, n_col,  Up, Ui,
		    Ux, 1, Control) ;
		if (print_level >= 3) mexPrintf ("\nP: ") ;
	        (void) umfpack_dl_report_perm (n_row, P, Control) ;
		if (print_level >= 3) mexPrintf ("\nQ: ") ;
	        (void) umfpack_dl_report_perm (n_col, Q, Control) ;
	    }

	    mxFree (P) ;
	    mxFree (Q) ;

	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* return the symbolic factorization */
	/* ------------------------------------------------------------------ */

	Q = (UF_long *) mxMalloc (n_col * sizeof (UF_long)) ;
	P = (UF_long *) mxMalloc (n_row * sizeof (UF_long)) ;
	Front_npivcol = (UF_long *) mxMalloc ((nn+1) * sizeof (UF_long)) ;
	Front_parent = (UF_long *) mxMalloc ((nn+1) * sizeof (UF_long)) ;
	Front_1strow = (UF_long *) mxMalloc ((nn+1) * sizeof (UF_long)) ;
	Front_leftmostdesc = (UF_long *) mxMalloc ((nn+1) * sizeof (UF_long)) ;
	Chain_start = (UF_long *) mxMalloc ((nn+1) * sizeof (UF_long)) ;
	Chain_maxrows = (UF_long *) mxMalloc ((nn+1) * sizeof (UF_long)) ;
	Chain_maxcols = (UF_long *) mxMalloc ((nn+1) * sizeof (UF_long)) ;

	if (A_is_complex)
	{
	    status = umfpack_zl_get_symbolic (&ignore1, &ignore2, &ignore3,
	        &nz, &nfronts, &nchains, P, Q, Front_npivcol,
	        Front_parent, Front_1strow, Front_leftmostdesc,
	        Chain_start, Chain_maxrows, Chain_maxcols, Symbolic) ;
	    umfpack_zl_free_symbolic (&Symbolic) ;
	}
	else
	{
	    status = umfpack_dl_get_symbolic (&ignore1, &ignore2, &ignore3,
	        &nz, &nfronts, &nchains, P, Q, Front_npivcol,
	        Front_parent, Front_1strow, Front_leftmostdesc,
	        Chain_start, Chain_maxrows, Chain_maxcols, Symbolic) ;
	    umfpack_dl_free_symbolic (&Symbolic) ;
	}

	if (status < 0)
	{
	    mxFree (P) ;
	    mxFree (Q) ;
	    mxFree (Front_npivcol) ;
	    mxFree (Front_parent) ;
	    mxFree (Front_1strow) ;
	    mxFree (Front_leftmostdesc) ;
	    mxFree (Chain_start) ;
	    mxFree (Chain_maxrows) ;
	    mxFree (Chain_maxcols) ;
	    error ("extracting symbolic factors failed", A_is_complex, nargout,
		pargout, Control, Info, status, do_info) ;
	    return ;
	}

	/* create sparse permutation matrix for P */
	pargout [0] = mxCreateSparse (n_row, n_row, n_row, mxREAL) ;
	Pp = (UF_long *) mxGetJc (pargout [0]) ;
	Pi = (UF_long *) mxGetIr (pargout [0]) ;
	Px = mxGetPr (pargout [0]) ;
	for (k = 0 ; k < n_row ; k++)
	{
	    Pp [k] = k ;
	    Px [k] = 1 ;
	    Pi [P [k]] = k ;
	}
	Pp [n_row] = n_row ;

	/* create sparse permutation matrix for Q */
	pargout [1] = mxCreateSparse (n_col, n_col, n_col, mxREAL) ;
	Qp = (UF_long *) mxGetJc (pargout [1]) ;
	Qi = (UF_long *) mxGetIr (pargout [1]) ;
	Qx = mxGetPr (pargout [1]) ;
	for (k = 0 ; k < n_col ; k++)
	{
	    Qp [k] = k ;
	    Qx [k] = 1 ;
	    Qi [k] = Q [k] ;
	}
	Qp [n_col] = n_col ;

	/* create Fr */
	pargout [2] = mxCreateDoubleMatrix (nfronts+1, 4, mxREAL) ;

	p1 = mxGetPr (pargout [2]) ;
	p2 = p1 + nfronts + 1 ;
	p3 = p2 + nfronts + 1 ;
	p4 = p3 + nfronts + 1 ;
	for (i = 0 ; i <= nfronts ; i++)
	{
	    /* convert parent, 1strow, and leftmostdesc to 1-based */
	    p1 [i] = (double) (Front_npivcol [i]) ;
	    p2 [i] = (double) (Front_parent [i] + 1) ;
	    p3 [i] = (double) (Front_1strow [i] + 1) ;
	    p4 [i] = (double) (Front_leftmostdesc [i] + 1) ;
	}

	/* create Ch */
	pargout [3] = mxCreateDoubleMatrix (nchains+1, 3, mxREAL) ;
	p1 = mxGetPr (pargout [3]) ;
	p2 = p1 + nchains + 1 ;
	p3 = p2 + nchains + 1 ;
	for (i = 0 ; i < nchains ; i++)
	{
	    p1 [i] = (double) (Chain_start [i] + 1) ;	/* convert to 1-based */
	    p2 [i] = (double) (Chain_maxrows [i]) ;
	    p3 [i] = (double) (Chain_maxcols [i]) ;
	}
	p1 [nchains] = Chain_start [nchains] + 1 ;
	p2 [nchains] = 0 ;
	p3 [nchains] = 0 ;

	mxFree (P) ;
	mxFree (Q) ;
	mxFree (Front_npivcol) ;
	mxFree (Front_parent) ;
	mxFree (Front_1strow) ;
	mxFree (Front_leftmostdesc) ;
	mxFree (Chain_start) ;
	mxFree (Chain_maxrows) ;
	mxFree (Chain_maxcols) ;
    }

    /* ---------------------------------------------------------------------- */
    /* report Info */
    /* ---------------------------------------------------------------------- */

    if (A_is_complex)
    {
	umfpack_zl_report_info (Control, Info) ;
    }
    else
    {
	umfpack_dl_report_info (Control, Info) ;
    }

    if (do_info > 0)
    {
	/* return Info */
	pargout [do_info] = mxCreateDoubleMatrix (1, UMFPACK_INFO, mxREAL) ;
	Out_Info = mxGetPr (pargout [do_info]) ;
	for (i = 0 ; i < UMFPACK_INFO ; i++)
	{
	    Out_Info [i] = Info [i] ;
	}
    }
}
