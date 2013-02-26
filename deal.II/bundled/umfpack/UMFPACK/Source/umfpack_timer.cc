/* ========================================================================== */
/* === umfpack_timer ======================================================== */
/* ========================================================================== */

/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

/*
    User-callable.  Returns the time in seconds used by the process.  BE
    CAREFUL:  if you compare the run time of UMFPACK with other sparse matrix
    packages, be sure to use the same timer.  See umfpack_timer.h for details.
    See umfpack_tictoc.h, which is the timer used internally by UMFPACK.
*/

#ifdef NO_TIMER

/* -------------------------------------------------------------------------- */
/* no timer used if -DNO_TIMER is defined at compile time */
/* -------------------------------------------------------------------------- */

double umfpack_timer ( void )
{
    return (0) ;
}

#else

#ifdef GETRUSAGE

/* -------------------------------------------------------------------------- */
/* use getrusage for accurate process times (and no overflow) */
/* -------------------------------------------------------------------------- */

/*
    This works under Solaris, SGI Irix, Linux, IBM RS 6000 (AIX), and Compaq
    Alpha.  It might work on other Unix systems, too.  Includes both the "user
    time" and the "system time".  The system time is the time spent by the
    operating system on behalf of the process, and thus should be charged to
    the process.
*/

#include <sys/time.h>
#include <sys/resource.h>

double umfpack_timer ( void )
{
    struct rusage ru ;
    double user_time, sys_time ;

    (void) getrusage (RUSAGE_SELF, &ru) ;

    user_time =
    ru.ru_utime.tv_sec			/* user time (seconds) */
    + 1e-6 * ru.ru_utime.tv_usec ;	/* user time (microseconds) */

    sys_time =
    ru.ru_stime.tv_sec			/* system time (seconds) */
    + 1e-6 * ru.ru_stime.tv_usec ;	/* system time (microseconds) */

    return (user_time + sys_time) ;
}

#else

/* -------------------------------------------------------------------------- */
/* Generic ANSI C: use the ANSI clock function */
/* -------------------------------------------------------------------------- */

/* This is portable, but may overflow.  On Sun Solaris, when compiling in */
/* 32-bit mode, the overflow occurs in only 2147 seconds (about 36 minutes). */

#include <time.h>

double umfpack_timer ( void )
{
    return (((double) (clock ( ))) / ((double) (CLOCKS_PER_SEC))) ;
}

#endif
#endif
