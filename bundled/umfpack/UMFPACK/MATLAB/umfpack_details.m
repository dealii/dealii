function umfpack_details
%UMFPACK_DETAILS details on all the options for using umfpack2 in MATLAB
%
% Factor or solve a sparse linear system, returning either the solution x to
% Ax=b or A'x'=b', the factorization LU=PAQ, or LU=P(R\A)Q.  A must be sparse.
% For the solve, A must be square and b must be a dense n-by-1 vector.  For LU
% factorization, A can be rectangular.  In both cases, A and/or b can be real
% or complex.
%
% UMFPACK analyzes the matrix and selects one of three strategies to factorize
% the matrix.  It first finds a set of k initial pivot entries of zero Markowitz
% cost.  This forms the first k rows and columns of L and U.  The remaining
% submatrix S is then analyzed, based on the symmetry of the nonzero pattern of
% the submatrix and the values on the diagaonal.  The strategies include:
%
%       (1) unsymmetric:  use a COLAMD pre-ordering, a column elimination tree
%           post-ordering, refine the column ordering during factorization,
%           and make no effort at selecting pivots on the diagonal.
%       (2) 2-by-2:  like the symmetric strategy (see below), except that local
%           row permutations are first made to attempt to place large entries
%           on the diagonal.
%       (3) symmetric:  use an AMD pre-ordering on the matrix S+S', an
%           elimination tree post-ordering, do not refine the column ordering
%           during factorization, and attempt to select pivots on the diagonal.
%
% Each of the following uses of umfpack2 (except for "Control = umfpack2") is
% stand-alone.  That is, no call to umfpack2 is required for any subsequent
% call.  In each usage, the Info output argument is optional.
%
% Example:
%
% [x, Info] = umfpack2 (A, '\', b) ;
% [x, Info] = umfpack2 (A, '\', b, Control) ;
% [x, Info] = umfpack2 (A, Qinit, '\', b, Control) ;
% [x, Info] = umfpack2 (A, Qinit, '\', b) ;
%
%       Solves Ax=b (similar to x = A\b in MATLAB).
%
% [x, Info] = umfpack2 (b, '/', A) ;
% [x, Info] = umfpack2 (b, '/', A, Control) ;
% [x, Info] = umfpack2 (b, '/', A, Qinit) ;
% [x, Info] = umfpack2 (b, '/', A, Qinit, Control) ;
%
%       Solves A'x'=b' (similar to x = b/A in MATLAB).
%
% [L, U, P, Q, R, Info] = umfpack2 (A) ;
% [L, U, P, Q, R, Info] = umfpack2 (A, Control) ;
% [L, U, P, Q, R, Info] = umfpack2 (A, Qinit) ;
% [L, U, P, Q, R, Info] = umfpack2 (A, Qinit, Control) ;
%
%       Returns the LU factorization of A.  P and Q are returned as permutation
%       matrices.  R is a diagonal sparse matrix of scale factors for the rows
%       of A, L is lower triangular, and U is upper triangular.  The
%       factorization is L*U = P*(R\A)*Q.  You can turn off scaling by setting
%       Control (17) to zero (in which case R = speye (m)), or by using the
%       following syntaxes (in which case Control (17) is ignored):
%
% [L, U, P, Q] = umfpack2 (A) ;
% [L, U, P, Q] = umfpack2 (A, Control) ;
% [L, U, P, Q] = umfpack2 (A, Qinit) ;
% [L, U, P, Q] = umfpack2 (A, Qinit, Control) ;
%
%       Same as above, except that no row scaling is performed.  The Info array
%       is not returned, either.
%
% [P1, Q1, Fr, Ch, Info] = umfpack2 (A, 'symbolic') ;
% [P1, Q1, Fr, Ch, Info] = umfpack2 (A, 'symbolic', Control) ;
% [P1, Q1, Fr, Ch, Info] = umfpack2 (A, Qinit, 'symbolic') ;
% [P1, Q1, Fr, Ch, Info] = umfpack2 (A, Qinit, 'symbolic', Control);
%
%       Performs only the fill-reducing column pre-ordering (including the
%       elimination tree post-ordering) and symbolic factorization.  Q1 is the
%       initial column permutation (either from colamd, amd, or the input
%       ordering Qinit), possibly followed by a column elimination tree post-
%       ordering or a symmetric elimination tree post-ordering, depending on
%       the strategy used.
%
%       For the unsymmetric strategy, P1 is the row ordering induced by Q1
%       (row-merge order). For the 2-by-2 strategy, P1 is the row ordering that
%       places large entries on the diagonal of P1*A*Q1.  For the symmetric
%       strategy, P1 = Q1.
%
%       Fr is a (nfr+1)-by-4 array containing information about each frontal
%       matrix, where nfr <= n is the number of frontal matrices.  Fr (:,1) is
%       the number of pivot columns in each front, and Fr (:,2) is the parent
%       of each front in the supercolumn elimination tree.  Fr (k,2) is zero if
%       k is a root.  The first Fr (1,1) columns of P1*A*Q1 are the pivot
%       columns for the first front, the next Fr (2,1) columns of P1*A*Q1
%       are the pivot columns for the second front, and so on.
%
%       For the unsymmetric strategy, Fr (:,3) is the row index of the first
%       row in P1*A*Q1 whose leftmost nonzero entry is in a pivot column for
%       the kth front.  Fr (:,4) is the leftmost descendent of the kth front.
%       Rows in the range Fr (Fr (k,4),3) to Fr (k+1,3)-1 form the entire set
%       of candidate pivot rows for the kth front (some of these will typically
%       have been selected as pivot rows of fronts Fr (k,3) to k-1, before the
%       factorization reaches the kth front.  If front k is a leaf node, then
%       Fr (k,4) is k.
%
%       Ch is a (nchains+1)-by-3 array containing information about each "chain"
%       (unifrontal sequence) of frontal matrices, and where nchains <= nfr
%       is the number of chains.  The ith chain consists of frontal matrices.
%       Chain (i,1) to Chain (i+1,1)-1, and the largest front in chain i is
%       Chain (i,2)-by-Chain (i,3).
%
%       This use of umfpack2 is not required to factor or solve a linear system
%       in MATLAB.  It analyzes the matrix A and provides information only.
%       The MATLAB statement "treeplot (Fr (:,2)')" plots the column elimination
%       tree.
%
% Control = umfpack2 ;
%
%       Returns a 20-by-1 vector of default parameter settings for umfpack2.
%
% umfpack_report (Control, Info) ;
%
%       Prints the current Control settings, and Info
%
% det = umfpack2 (A, 'det') ;
% [det dexp] = umfpack2 (A, 'det') ;
%
%       Computes the determinant of A.  The 2nd form returns the determinant
%       in the form det*10^dexp, where det is in the range +/- 1 to 10,
%       which helps to avoid overflow/underflow when dexp is out of range of
%       normal floating-point numbers.
%
% If present, Qinit is a user-supplied 1-by-n permutation vector.  It is an
% initial fill-reducing column pre-ordering for A; if not present, then colamd
% or amd are used instead.  If present, Control is a user-supplied 20-by-1
% array.   Control and Info are optional; if Control is not present, defaults
% are used.  If a Control entry is NaN, then the default is used for that entry.
%
%
% Copyright 1995-2007 by Timothy A. Davis, University of Florida.
% All Rights Reserved.
% UMFPACK is available under alternate licenses, contact T. Davis for details.
%
% UMFPACK License:
%
%   Your use or distribution of UMFPACK or any modified version of
%   UMFPACK implies that you agree to this License.
%
%   This library is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation; either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
%   USA
%
%   Permission is hereby granted to use or copy this program under the
%   terms of the GNU LGPL, provided that the Copyright, this License,
%   and the Availability of the original version is retained on all copies.
%   User documentation of any code that uses this code or any modified
%   version of this code must cite the Copyright, this License, the
%   Availability note, and "Used by permission." Permission to modify
%   the code and to distribute modified code is granted, provided the
%   Copyright, this License, and the Availability note are retained,
%   and a notice that the code was modified is included.
%
% Availability: http://www.cise.ufl.edu/research/sparse/umfpack
%
% See also umfpack, umfpack2, umfpack_make, umfpack_report,
%    umfpack_demo, and umfpack_simple.

more on
help umfpack_details
more off

