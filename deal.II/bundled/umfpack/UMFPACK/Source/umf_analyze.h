/* -------------------------------------------------------------------------- */
/* UMFPACK Copyright (c) Timothy A. Davis, CISE,                              */
/* Univ. of Florida.  All Rights Reserved.  See ../Doc/License for License.   */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_analyze
(
    Int n_row,
    Int n_col,
    Int Ai [ ],
    Int Ap [ ],
    Int Up [ ],
    Int fixQ,
    Int Front_ncols [ ],
    Int W [ ],
    Int Link [ ],
    Int Front_nrows [ ],
    Int Front_npivcol [ ],
    Int Front_parent [ ],
    Int *nfr_out,
    Int *p_ncompactions
) ;
