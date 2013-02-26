/* ========================================================================== */
/* === UMF_multicompile ===================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* This file is not needed if you have the Unix/Linux "make" command for
 * compiling UMFPACK.  Microsoft Visual Studio cannot be configured to compile
 * one file multiple times, with different -D flags.  In this case, you can
 * use this file instead.  To use this file, see the Demo/simple_compile file.
 *
 * This file includes the following source files:
 *
 *	umf_ltsolve.c
 *	umf_utsolve.c
 *	umf_triplet.c
 *	umf_assemble.c
 *	umf_store_lu.c
 *	umfpack_solve.c
 *
 * This file simply compiles the above files with different pre-#define'd flags,
 * by defining the flags and then #include'ing the source files themselves.
 * This is a rather unconventional approach, since by convention #include is
 * supposed to be used with *.h files not *.c.  However, it is one way of
 * working around the limitations of Microsoft Visual Studio.
 *
 * You still need to compile all files separately as well, with none of the
 * pre-#define'd terms listed below.
 */

/* compile the complex conjugate forward/backsolves */
#define CONJUGATE_SOLVE
#include "umf_ltsolve.c"
#include "umf_utsolve.c"

/* compile umf_triplet with DO_MAP, DO_VALUES and DO_MAP, and just DO_VALUES */
#define DO_MAP
#include "umf_triplet.c"
#define DO_VALUES
#include "umf_triplet.c"
#undef DO_MAP
#include "umf_triplet.c"

/* compile the FIXQ version of umf_assemble */
#define FIXQ
#include "umf_assemble.c"

/* compile the DROP version of umf_store_lu */
#define DROP
#include "umf_store_lu.c"

/* compile umfpack_wsolve */
#define WSOLVE
#include "umfpack_solve.c"
