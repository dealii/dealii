/* ========================================================================== */
/* === umfpack_global ======================================================= */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/* prototypes for global variables, and basic operators for complex values  */

#ifndef EXTERN
#define EXTERN extern
#endif

EXTERN double (*umfpack_hypot) (double, double) ;
EXTERN int (*umfpack_divcomplex) (double, double, double, double, double *, double *) ;

double umf_hypot (double x, double y) ;
int umf_divcomplex (double, double, double, double, double *, double *) ;

