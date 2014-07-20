/* ========================================================================== */
/* === umfpack_timer ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

double umfpack_timer ( void ) ;

/*
Syntax (for all versions: di, dl, zi, and zl):

    #include "umfpack.h"
    double t ;
    t = umfpack_timer ( ) ;

Purpose:

    Returns the CPU time used by the process.  Includes both "user" and "system"
    time (the latter is time spent by the system on behalf of the process, and
    is thus charged to the process).  It does not return the wall clock time.
    See umfpack_tic and umfpack_toc (the file umfpack_tictoc.h) for the timer
    used internally by UMFPACK.

    This routine uses the Unix getrusage routine, if available.  It is less
    subject to overflow than the ANSI C clock routine.  If getrusage is not
    available, the portable ANSI C clock routine is used instead.
    Unfortunately, clock ( ) overflows if the CPU time exceeds 2147 seconds
    (about 36 minutes) when sizeof (clock_t) is 4 bytes.  If you have getrusage,
    be sure to compile UMFPACK with the -DGETRUSAGE flag set; see umf_config.h
    and the User Guide for details.  Even the getrusage routine can overlow.

Arguments:

    None.
*/
