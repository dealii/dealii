/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#ifndef _UMF_MALLOC
#define _UMF_MALLOC

#if defined (UMF_MALLOC_COUNT) || !defined (NDEBUG)

#ifndef EXTERN
#define EXTERN extern
#endif

GLOBAL EXTERN Int UMF_malloc_count ;
#endif

GLOBAL void *UMF_malloc
(
    Int n_objects,
    size_t size_of_object
) ;

#endif
