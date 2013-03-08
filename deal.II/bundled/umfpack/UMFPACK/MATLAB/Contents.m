%Contents of the UMFPACK sparse matrix toolbox:
%
% umfpack2        - computes x=A\b, x=A/b, or lu (A) for a sparse matrix A
% umfpack_make    - to compile umfpack2 for use in MATLAB
% umfpack_install - to compile and install umfpack2 and amd2 for use in MATLAB
% umfpack_details - details on all the options for using umfpack2 in MATLAB
% umfpack_report  - prints optional control settings and statistics
% umfpack_demo    - a lenghty demo
% umfpack_simple  - a simple demo
% umfpack_btf     - factorize A using a block triangular form
% umfpack_solve   - x = A\b or x = b/A
% lu_normest      - estimates norm (L*U-A, 1) without forming L*U-A
% luflop          - given L and U, computes # of flops required to compute them
% umfpack_test    - for testing umfpack2 (requires UFget)
%
% Example:
%   x = umfpack2 (A, '\', b) ;   % same as x = A\b, if A square and unsymmetric
%
% See also these built-in functions:
% umfpack           the built-in version of UMFPACK
% amd               symmetric minimum degree ordering
% colamd            unsymmetric column approx minimum degree ordering
% symamd            symmetric approx minimum degree ordering, based on colamd
%
% NOTE: UMFPACK is a built-in function in MATLAB, used in x=A\b.  This is the
% source code for the built-in umfpack, but the MATLAB function has been renamed
% here to umfpack2, to avoid a filename clash with itself.
%
% Copyright 1995-2007 by Timothy A. Davis.
% All Rights Reserved.  Type umfpack_details for License.

help Contents
