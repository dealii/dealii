/* ========================================================================== */
/* === umfpack.h ============================================================ */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    This is the umfpack.h include file, and should be included in all user code
    that uses UMFPACK.  Do not include any of the umf_* header files in user
    code.  All routines in UMFPACK starting with "umfpack_" are user-callable.
    All other routines are prefixed "umf_XY_", (where X is d or z, and Y is
    i or l) and are not user-callable.
*/

#ifndef UMFPACK_H
#define UMFPACK_H

/* -------------------------------------------------------------------------- */
/* Make it easy for C++ programs to include UMFPACK */
/* -------------------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

/* define UF_long */
#include "UFconfig.h"

/* -------------------------------------------------------------------------- */
/* size of Info and Control arrays */
/* -------------------------------------------------------------------------- */

/* These might be larger in future versions, since there are only 3 unused
 * entries in Info, and no unused entries in Control. */

#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20

/* -------------------------------------------------------------------------- */
/* User-callable routines */
/* -------------------------------------------------------------------------- */

/* Primary routines: */
#include "umfpack_symbolic.h"
#include "umfpack_numeric.h"
#include "umfpack_solve.h"
#include "umfpack_free_symbolic.h"
#include "umfpack_free_numeric.h"

/* Alternative routines: */
#include "umfpack_defaults.h"
#include "umfpack_qsymbolic.h"
#include "umfpack_wsolve.h"

/* Matrix manipulation routines: */
#include "umfpack_triplet_to_col.h"
#include "umfpack_col_to_triplet.h"
#include "umfpack_transpose.h"
#include "umfpack_scale.h"

/* Getting the contents of the Symbolic and Numeric opaque objects: */
#include "umfpack_get_lunz.h"
#include "umfpack_get_numeric.h"
#include "umfpack_get_symbolic.h"
#include "umfpack_save_numeric.h"
#include "umfpack_load_numeric.h"
#include "umfpack_save_symbolic.h"
#include "umfpack_load_symbolic.h"
#include "umfpack_get_determinant.h"

/* Reporting routines (the above 14 routines print nothing): */
#include "umfpack_report_status.h"
#include "umfpack_report_info.h"
#include "umfpack_report_control.h"
#include "umfpack_report_matrix.h"
#include "umfpack_report_triplet.h"
#include "umfpack_report_vector.h"
#include "umfpack_report_symbolic.h"
#include "umfpack_report_numeric.h"
#include "umfpack_report_perm.h"

/* Utility routines: */
#include "umfpack_timer.h"
#include "umfpack_tictoc.h"

/* AMD */
#include <amd.h>

/* global function pointers */
#include "umfpack_global.h"

/* -------------------------------------------------------------------------- */
/* Version, copyright, and license */
/* -------------------------------------------------------------------------- */

#define UMFPACK_VERSION "UMFPACK V5.1.0 (May 31, 2007)"

#define UMFPACK_COPYRIGHT \
"UMFPACK:  Copyright (c) 2005-2006 by Timothy A. Davis.  All Rights Reserved.\n"

#define UMFPACK_LICENSE_PART1 \
"\nUMFPACK License:\n" \
"\n" \
"   UMFPACK is available under alternate licenses,\n" \
"   contact T. Davis for details.\n" \
"\n" \
"   Your use or distribution of UMFPACK or any modified version of\n" \
"   UMFPACK implies that you agree to this License.\n" \
"\n" \
"   This library is free software; you can redistribute it and/or\n" \
"   modify it under the terms of the GNU Lesser General Public\n" \
"   License as published by the Free Software Foundation; either\n" \
"   version 2.1 of the License, or (at your option) any later version.\n" \
"\n" \
"   This library is distributed in the hope that it will be useful,\n" \
"   but WITHOUT ANY WARRANTY; without even the implied warranty of\n" \
"   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n" \
"   Lesser General Public License for more details.\n" \
"\n" \
"   You should have received a copy of the GNU Lesser General Public\n" \
"   License along with this library; if not, write to the Free Software\n" \
"   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301\n" \
"   USA\n" \

#define UMFPACK_LICENSE_PART2 \
"\n" \
"   Permission is hereby granted to use or copy this program under the\n" \
"   terms of the GNU LGPL, provided that the Copyright, this License,\n" \
"   and the Availability of the original version is retained on all copies.\n" \
"   User documentation of any code that uses this code or any modified\n" \
"   version of this code must cite the Copyright, this License, the\n" \
"   Availability note, and \"Used by permission.\" Permission to modify\n" \
"   the code and to distribute modified code is granted, provided the\n" \
"   Copyright, this License, and the Availability note are retained,\n" \
"   and a notice that the code was modified is included.\n"

#define UMFPACK_LICENSE_PART3 \
"\n" \
"Availability: http://www.cise.ufl.edu/research/sparse/umfpack\n" \
"\n"

/* UMFPACK Version 4.5 and later will include the following definitions.
 * As an example, to test if the version you are using is 4.5 or later:
 *
 * #ifdef UMFPACK_VER
 *	if (UMFPACK_VER >= UMFPACK_VER_CODE (4,5)) ...
 * #endif
 *
 * This also works during compile-time:
 *
 *	#if defined(UMFPACK_VER) && (UMFPACK >= UMFPACK_VER_CODE (4,5))
 *	    printf ("This is version 4.5 or later\n") ;
 *	#else
 *	    printf ("This is an early version\n") ;
 *	#endif
 *
 * Versions 4.4 and earlier of UMFPACK do not include a #define'd version
 * number, although they do include the UMFPACK_VERSION string, defined
 * above.
 */

#define UMFPACK_DATE "May 31, 2007"
#define UMFPACK_VER_CODE(main,sub) ((main) * 1000 + (sub))
#define UMFPACK_MAIN_VERSION 5
#define UMFPACK_SUB_VERSION 1
#define UMFPACK_SUBSUB_VERSION 0
#define UMFPACK_VER UMFPACK_VER_CODE(UMFPACK_MAIN_VERSION,UMFPACK_SUB_VERSION)

/* -------------------------------------------------------------------------- */
/* contents of Info */
/* -------------------------------------------------------------------------- */

/* Note that umfpack_report.m must coincide with these definitions.  S is
 * the submatrix of A after removing row/col singletons and empty rows/cols. */

/* returned by all routines that use Info: */
#define UMFPACK_STATUS 0	/* UMFPACK_OK, or other result */
#define UMFPACK_NROW 1		/* n_row input value */
#define UMFPACK_NCOL 16		/* n_col input value */
#define UMFPACK_NZ 2		/* # of entries in A */

/* computed in UMFPACK_*symbolic and UMFPACK_numeric: */
#define UMFPACK_SIZE_OF_UNIT 3		/* sizeof (Unit) */

/* computed in UMFPACK_*symbolic: */
#define UMFPACK_SIZE_OF_INT 4		/* sizeof (int) */
#define UMFPACK_SIZE_OF_LONG 5		/* sizeof (UF_long) */
#define UMFPACK_SIZE_OF_POINTER 6	/* sizeof (void *) */
#define UMFPACK_SIZE_OF_ENTRY 7		/* sizeof (Entry), real or complex */
#define UMFPACK_NDENSE_ROW 8		/* number of dense rows */
#define UMFPACK_NEMPTY_ROW 9		/* number of empty rows */
#define UMFPACK_NDENSE_COL 10		/* number of dense rows */
#define UMFPACK_NEMPTY_COL 11		/* number of empty rows */
#define UMFPACK_SYMBOLIC_DEFRAG 12	/* # of memory compactions */
#define UMFPACK_SYMBOLIC_PEAK_MEMORY 13	/* memory used by symbolic analysis */
#define UMFPACK_SYMBOLIC_SIZE 14	/* size of Symbolic object, in Units */
#define UMFPACK_SYMBOLIC_TIME 15	/* time (sec.) for symbolic analysis */
#define UMFPACK_SYMBOLIC_WALLTIME 17	/* wall clock time for sym. analysis */
#define UMFPACK_STRATEGY_USED 18	/* strategy used: sym, unsym, 2by2 */
#define UMFPACK_ORDERING_USED 19	/* ordering used: colamd, amd, given */
#define UMFPACK_QFIXED 31		/* whether Q is fixed or refined */
#define UMFPACK_DIAG_PREFERRED 32	/* whether diagonal pivoting attempted*/
#define UMFPACK_PATTERN_SYMMETRY 33	/* symmetry of pattern of S */
#define UMFPACK_NZ_A_PLUS_AT 34		/* nnz (S+S'), excl. diagonal */
#define UMFPACK_NZDIAG 35		/* nnz (diag (S)) */

/* AMD statistics, computed in UMFPACK_*symbolic: */
#define UMFPACK_SYMMETRIC_LUNZ 36	/* nz in L+U, if AMD ordering used */
#define UMFPACK_SYMMETRIC_FLOPS 37	/* flops for LU, if AMD ordering used */
#define UMFPACK_SYMMETRIC_NDENSE 38	/* # of "dense" rows/cols in S+S' */
#define UMFPACK_SYMMETRIC_DMAX 39	/* max nz in cols of L, for AMD */

/* statistics for 2-by-2 strategy */
#define UMFPACK_2BY2_NWEAK 51		    /* number of weak diagonal entries*/
#define UMFPACK_2BY2_UNMATCHED 52	    /* # of weak diagonals not matched*/
#define UMFPACK_2BY2_PATTERN_SYMMETRY 53    /* symmetry of pattern of P*S */
#define UMFPACK_2BY2_NZ_PA_PLUS_PAT 54	    /* nz in PS+(PS)' */
#define UMFPACK_2BY2_NZDIAG 55		    /* nz on diagonal of PS+(PS)' */

/* statistcs for singleton pruning */
#define UMFPACK_COL_SINGLETONS 56	/* # of column singletons */
#define UMFPACK_ROW_SINGLETONS 57	/* # of row singletons */
#define UMFPACK_N2 58			/* size of S */
#define UMFPACK_S_SYMMETRIC 59		/* 1 if S square and symmetricly perm.*/

/* estimates computed in UMFPACK_*symbolic: */
#define UMFPACK_NUMERIC_SIZE_ESTIMATE 20    /* final size of Numeric->Memory */
#define UMFPACK_PEAK_MEMORY_ESTIMATE 21	    /* for symbolic & numeric */
#define UMFPACK_FLOPS_ESTIMATE 22	    /* flop count */
#define UMFPACK_LNZ_ESTIMATE 23		    /* nz in L, incl. diagonal */
#define UMFPACK_UNZ_ESTIMATE 24		    /* nz in U, incl. diagonal */
#define UMFPACK_VARIABLE_INIT_ESTIMATE 25   /* initial size of Numeric->Memory*/
#define UMFPACK_VARIABLE_PEAK_ESTIMATE 26   /* peak size of Numeric->Memory */
#define UMFPACK_VARIABLE_FINAL_ESTIMATE 27  /* final size of Numeric->Memory */
#define UMFPACK_MAX_FRONT_SIZE_ESTIMATE 28  /* max frontal matrix size */
#define UMFPACK_MAX_FRONT_NROWS_ESTIMATE 29 /* max # rows in any front */
#define UMFPACK_MAX_FRONT_NCOLS_ESTIMATE 30 /* max # columns in any front */

/* exact values, (estimates shown above) computed in UMFPACK_numeric: */
#define UMFPACK_NUMERIC_SIZE 40		    /* final size of Numeric->Memory */
#define UMFPACK_PEAK_MEMORY 41		    /* for symbolic & numeric */
#define UMFPACK_FLOPS 42		    /* flop count */
#define UMFPACK_LNZ 43			    /* nz in L, incl. diagonal */
#define UMFPACK_UNZ 44			    /* nz in U, incl. diagonal */
#define UMFPACK_VARIABLE_INIT 45	    /* initial size of Numeric->Memory*/
#define UMFPACK_VARIABLE_PEAK 46	    /* peak size of Numeric->Memory */
#define UMFPACK_VARIABLE_FINAL 47	    /* final size of Numeric->Memory */
#define UMFPACK_MAX_FRONT_SIZE 48	    /* max frontal matrix size */
#define UMFPACK_MAX_FRONT_NROWS 49	    /* max # rows in any front */
#define UMFPACK_MAX_FRONT_NCOLS 50	    /* max # columns in any front */

/* computed in UMFPACK_numeric: */
#define UMFPACK_NUMERIC_DEFRAG 60	    /* # of garbage collections */
#define UMFPACK_NUMERIC_REALLOC 61	    /* # of memory reallocations */
#define UMFPACK_NUMERIC_COSTLY_REALLOC 62   /* # of costlly memory realloc's */
#define UMFPACK_COMPRESSED_PATTERN 63	    /* # of integers in LU pattern */
#define UMFPACK_LU_ENTRIES 64		    /* # of reals in LU factors */
#define UMFPACK_NUMERIC_TIME 65		    /* numeric factorization time */
#define UMFPACK_UDIAG_NZ 66		    /* nz on diagonal of U */
#define UMFPACK_RCOND 67		    /* est. reciprocal condition # */
#define UMFPACK_WAS_SCALED 68		    /* none, max row, or sum row */
#define UMFPACK_RSMIN 69		    /* min (max row) or min (sum row) */
#define UMFPACK_RSMAX 70		    /* max (max row) or max (sum row) */
#define UMFPACK_UMIN 71			    /* min abs diagonal entry of U */
#define UMFPACK_UMAX 72			    /* max abs diagonal entry of U */
#define UMFPACK_ALLOC_INIT_USED 73	    /* alloc_init parameter used */
#define UMFPACK_FORCED_UPDATES 74	    /* # of forced updates */
#define UMFPACK_NUMERIC_WALLTIME 75	    /* numeric wall clock time */
#define UMFPACK_NOFF_DIAG 76		    /* number of off-diagonal pivots */

#define UMFPACK_ALL_LNZ 77		    /* nz in L, if no dropped entries */
#define UMFPACK_ALL_UNZ 78		    /* nz in U, if no dropped entries */
#define UMFPACK_NZDROPPED 79		    /* # of dropped small entries */

/* computed in UMFPACK_solve: */
#define UMFPACK_IR_TAKEN 80	    /* # of iterative refinement steps taken */
#define UMFPACK_IR_ATTEMPTED 81	    /* # of iter. refinement steps attempted */
#define UMFPACK_OMEGA1 82	    /* omega1, sparse backward error estimate */
#define UMFPACK_OMEGA2 83	    /* omega2, sparse backward error estimate */
#define UMFPACK_SOLVE_FLOPS 84	    /* flop count for solve */
#define UMFPACK_SOLVE_TIME 85	    /* solve time (seconds) */
#define UMFPACK_SOLVE_WALLTIME 86   /* solve time (wall clock, seconds) */

/* Info [87, 88, 89] unused */

/* Unused parts of Info may be used in future versions of UMFPACK. */

/* -------------------------------------------------------------------------- */

/* Info [UMFPACK_ORDERING_USED] is one of the following: */
#define UMFPACK_ORDERING_COLAMD 0	/* COLAMD(A) */
#define UMFPACK_ORDERING_AMD 1		/* AMD(A+A') */
#define UMFPACK_ORDERING_GIVEN 2	/* Q is provided on input */

/* -------------------------------------------------------------------------- */
/* contents of Control */
/* -------------------------------------------------------------------------- */

/* used in all UMFPACK_report_* routines: */
#define UMFPACK_PRL 0			/* print level */

/* used in UMFPACK_*symbolic only: */
#define UMFPACK_DENSE_ROW 1		/* dense row parameter */
#define UMFPACK_DENSE_COL 2		/* dense col parameter */
#define UMFPACK_BLOCK_SIZE 4		/* BLAS-3 block size */
#define UMFPACK_STRATEGY 5		/* auto, symmetric, unsym., or 2by2 */
#define UMFPACK_2BY2_TOLERANCE 12	/* 2-by-2 pivot tolerance */
#define UMFPACK_FIXQ 13			/* -1: no fixQ, 0: default, 1: fixQ */
#define UMFPACK_AMD_DENSE 14		/* for AMD ordering */
#define UMFPACK_AGGRESSIVE 19		/* whether or not to use aggressive
					 * absorption in AMD and COLAMD */

/* used in UMFPACK_numeric only: */
#define UMFPACK_PIVOT_TOLERANCE 3	/* threshold partial pivoting setting */
#define UMFPACK_ALLOC_INIT 6		/* initial allocation ratio */
#define UMFPACK_SYM_PIVOT_TOLERANCE 15	/* threshold, only for diag. entries */
#define UMFPACK_SCALE 16		/* what row scaling to do */
#define UMFPACK_FRONT_ALLOC_INIT 17	/* frontal matrix allocation ratio */
#define UMFPACK_DROPTOL 18		/* drop tolerance for entries in L,U */

/* used in UMFPACK_*solve only: */
#define UMFPACK_IRSTEP 7		/* max # of iterative refinements */

/* compile-time settings - Control [8..11] cannot be changed at run time: */
#define UMFPACK_COMPILED_WITH_BLAS 8	    /* uses the BLAS */
#define UMFPACK_COMPILED_FOR_MATLAB 9	    /* 1 if MATLAB mexFunction, etc. */
#define UMFPACK_COMPILED_WITH_GETRUSAGE 10  /* uses getrusage timer, or not */
#define UMFPACK_COMPILED_IN_DEBUG_MODE 11   /* debugging enabled (very slow!) */

/* -------------------------------------------------------------------------- */

/* Control [UMFPACK_STRATEGY] is one of the following: */
#define UMFPACK_STRATEGY_AUTO 0		/* use sym. or unsym. strategy */
#define UMFPACK_STRATEGY_UNSYMMETRIC 1	/* COLAMD(A), coletree postorder,
					   not prefer diag*/
#define UMFPACK_STRATEGY_2BY2 2		/* AMD(PA+PA'), no coletree postorder,
					   prefer diag(PA) where P is pseudo
					   max transversal */
#define UMFPACK_STRATEGY_SYMMETRIC 3	/* AMD(A+A'), no coletree postorder,
					   prefer diagonal */

/* Control [UMFPACK_SCALE] is one of the following: */
#define UMFPACK_SCALE_NONE 0	/* no scaling */
#define UMFPACK_SCALE_SUM 1	/* default: divide each row by sum (abs (row))*/
#define UMFPACK_SCALE_MAX 2	/* divide each row by max (abs (row)) */

/* -------------------------------------------------------------------------- */
/* default values of Control: */
/* -------------------------------------------------------------------------- */

#define UMFPACK_DEFAULT_PRL 1
#define UMFPACK_DEFAULT_DENSE_ROW 0.2
#define UMFPACK_DEFAULT_DENSE_COL 0.2
#define UMFPACK_DEFAULT_PIVOT_TOLERANCE 0.1
#define UMFPACK_DEFAULT_2BY2_TOLERANCE 0.01
#define UMFPACK_DEFAULT_SYM_PIVOT_TOLERANCE 0.001
#define UMFPACK_DEFAULT_BLOCK_SIZE 32
#define UMFPACK_DEFAULT_ALLOC_INIT 0.7
#define UMFPACK_DEFAULT_FRONT_ALLOC_INIT 0.5
#define UMFPACK_DEFAULT_IRSTEP 2
#define UMFPACK_DEFAULT_SCALE UMFPACK_SCALE_SUM
#define UMFPACK_DEFAULT_STRATEGY UMFPACK_STRATEGY_AUTO
#define UMFPACK_DEFAULT_AMD_DENSE AMD_DEFAULT_DENSE
#define UMFPACK_DEFAULT_FIXQ 0
#define UMFPACK_DEFAULT_AGGRESSIVE 1
#define UMFPACK_DEFAULT_DROPTOL 0

/* default values of Control may change in future versions of UMFPACK. */

/* -------------------------------------------------------------------------- */
/* status codes */
/* -------------------------------------------------------------------------- */

#define UMFPACK_OK (0)

/* status > 0 means a warning, but the method was successful anyway. */
/* A Symbolic or Numeric object was still created. */
#define UMFPACK_WARNING_singular_matrix (1)

/* The following warnings were added in umfpack_*_get_determinant */
#define UMFPACK_WARNING_determinant_underflow (2)
#define UMFPACK_WARNING_determinant_overflow (3)

/* status < 0 means an error, and the method was not successful. */
/* No Symbolic of Numeric object was created. */
#define UMFPACK_ERROR_out_of_memory (-1)
#define UMFPACK_ERROR_invalid_Numeric_object (-3)
#define UMFPACK_ERROR_invalid_Symbolic_object (-4)
#define UMFPACK_ERROR_argument_missing (-5)
#define UMFPACK_ERROR_n_nonpositive (-6)
#define UMFPACK_ERROR_invalid_matrix (-8)
#define UMFPACK_ERROR_different_pattern (-11)
#define UMFPACK_ERROR_invalid_system (-13)
#define UMFPACK_ERROR_invalid_permutation (-15)
#define UMFPACK_ERROR_internal_error (-911) /* yes, call me if you get this! */
#define UMFPACK_ERROR_file_IO (-17)

/* -------------------------------------------------------------------------- */
/* solve codes */
/* -------------------------------------------------------------------------- */

/* Solve the system ( )x=b, where ( ) is defined below.  "t" refers to the */
/* linear algebraic transpose (complex conjugate if A is complex), or the (') */
/* operator in MATLAB.  "at" refers to the array transpose, or the (.') */
/* operator in MATLAB. */

#define UMFPACK_A	(0)	/* Ax=b    */
#define UMFPACK_At	(1)	/* A'x=b   */
#define UMFPACK_Aat	(2)	/* A.'x=b  */

#define UMFPACK_Pt_L	(3)	/* P'Lx=b  */
#define UMFPACK_L	(4)	/* Lx=b    */
#define UMFPACK_Lt_P	(5)	/* L'Px=b  */
#define UMFPACK_Lat_P	(6)	/* L.'Px=b */
#define UMFPACK_Lt	(7)	/* L'x=b   */
#define UMFPACK_Lat	(8)	/* L.'x=b  */

#define UMFPACK_U_Qt	(9)	/* UQ'x=b  */
#define UMFPACK_U	(10)	/* Ux=b    */
#define UMFPACK_Q_Ut	(11)	/* QU'x=b  */
#define UMFPACK_Q_Uat	(12)	/* QU.'x=b */
#define UMFPACK_Ut	(13)	/* U'x=b   */
#define UMFPACK_Uat	(14)	/* U.'x=b  */

/* -------------------------------------------------------------------------- */

/* Integer constants are used for status and solve codes instead of enum */
/* to make it easier for a Fortran code to call UMFPACK. */

#ifdef __cplusplus
}
#endif

#endif /* UMFPACK_H */
