/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_singletons
(
    Int n_row,
    Int n_col,
    const Int Ap [ ],
    const Int Ai [ ],
    const Int Quser [ ],
    Int strategy,
    Int Cdeg [ ],
    Int Cperm [ ],
    Int Rdeg [ ],
    Int Rperm [ ],
    Int InvRperm [ ],
    Int *n1,
    Int *n1c,
    Int *n1r,
    Int *nempty_col,
    Int *nempty_row,
    Int *is_sym,
    Int *max_rdeg,
    Int Rp [ ],
    Int Ri [ ],
    Int W [ ],
    Int Next [ ]
) ;
