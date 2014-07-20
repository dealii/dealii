function [p, Info] = amd (A, Control)
%AMD Approximate minimum degree permutation.
%    P = AMD (S) returns the approximate minimum degree permutation vector for
%    the sparse matrix C = S+S'.  The Cholesky factorization of C (P,P), or
%    S (P,P), tends to be sparser than that of C or S.  AMD tends to be faster
%    than SYMMMD and SYMAMD, and tends to return better orderings than SYMMMD.
%    S must be square. If S is full, amd (S) is equivalent to amd (sparse (S)).
%
%    Usage:  P = amd (S) ;                   % finds the ordering
%            [P, Info] = amd (S, Control) ;  % optional parameters & statistics
%            Control = amd ;                 % returns default parameters
%            amd ;                           % prints default parameters.
%
%       Control (1); If S is n-by-n, then rows/columns with more than
%           max (16, (Control (1))* sqrt(n)) entries in S+S' are considered
%           "dense", and ignored during ordering.  They are placed last in the
%           output permutation.  The default is 10.0 if Control is not present.
%       Control (2): If nonzero, then aggressive absorption is performed.
%           This is the default if Control is not present.
%       Control (3): If nonzero, print statistics about the ordering.
%
%       Info (1): status (0: ok, -1: out of memory, -2: matrix invalid)
%       Info (2): n = size (A,1)
%       Info (3): nnz (A)
%       Info (4): the symmetry of the matrix S (0.0 means purely unsymmetric,
%           1.0 means purely symmetric).  Computed as:
%           B = tril (S, -1) + triu (S, 1) ; symmetry = nnz (B & B') / nnz (B);
%       Info (5): nnz (diag (S))
%       Info (6): nnz in S+S', excluding the diagonal (= nnz (B+B'))
%       Info (7): number "dense" rows/columns in S+S'
%       Info (8): the amount of memory used by AMD, in bytes
%       Info (9): the number of memory compactions performed by AMD
%
%    The following statistics are slight upper bounds because of the
%    approximate degree in AMD.  The bounds are looser if "dense" rows/columns
%    are ignored during ordering (Info (7) > 0).  The statistics are for a
%    subsequent factorization of the matrix C (P,P).  The LU factorization
%    statistics assume no pivoting.
%
%       Info (10): the number of nonzeros in L, excluding the diagonal
%       Info (11): the number of divide operations for LL', LDL', or LU
%       Info (12): the number of multiply-subtract pairs for LL' or LDL'
%       Info (13): the number of multiply-subtract pairs for LU
%       Info (14): the max # of nonzeros in any column of L (incl. diagonal)
%       Info (15:20): unused, reserved for future use
%
%    An assembly tree post-ordering is performed, which is typically the same
%    as an elimination tree post-ordering.  It is not always identical because
%    of the approximate degree update used, and because "dense" rows/columns
%    do not take part in the post-order.  It well-suited for a subsequent
%    "chol", however.  If you require a precise elimination tree post-ordering,
%    then do:
%
%       P = amd (S) ;
%       C = spones (S) + spones (S') ;  % skip this if S already symmetric
%       [ignore, Q] = sparsfun ('symetree', C (P,P)) ;
%       P = P (Q) ;
%
% --------------------------------------------------------------------------
% AMD Version 1.1 (Jan. 21, 2004), Copyright (c) 2004 by Timothy A. Davis,
% Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.
% email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.
% web: http://www.cise.ufl.edu/research/sparse/amd
% --------------------------------------------------------------------------
%
%    Acknowledgements: This work was supported by the National Science
%       Foundation, under grants ASC-9111263, DMS-9223088, and CCR-0203270.
%
%    See also COLMMD, COLAMD, COLPERM, SYMAMD, SYMMMD, SYMRCM.

more on
help amd
more off
error ('amd mexFunction not found!  Type "amd_make" in MATLAB to compile amd');
