/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL void UMF_2by2
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    const double Ax [ ],
#ifdef COMPLEX
    const double Az [ ],
#endif
    double tol,
    Int scale,
    Int Cperm1 [ ],
#ifndef NDEBUG
    Int Rperm1 [ ],
#endif
    Int InvRperm [ ],
    Int n1,
    Int nempty,
    Int Degree [ ],
    Int P [ ],
    Int *p_nweak,
    Int *p_nmatched,
    Int Ri [ ],
    Int Rp [ ],
    double Rs [ ],
    Int Head [ ],
    Int Next [ ],
    Int Si [ ],
    Int Sp [ ]
) ;
