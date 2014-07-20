function x = umfpack_btf (A, b, Control)
%UMFPACK_BTF factorize A using a block triangular form
%
% Example:
%   x = umfpack_btf (A, b, Control)
%
% solve Ax=b by first permuting the matrix A to block triangular form via dmperm
% and then using UMFPACK to factorize each diagonal block.  Adjacent 1-by-1
% blocks are merged into a single upper triangular block, and solved via
% MATLAB's \ operator.  The Control parameter is optional (Type umfpack_details
% and umfpack_report for details on its use).  A must be square.
%
% See also umfpack, umfpack2, umfpack_details, dmperm

% Copyright 1995-2007 by Timothy A. Davis.

if (nargin < 2)
    help umfpack_btf
    error ('Usage: x = umfpack_btf (A, b, Control)') ;
end

[m n] = size (A) ;
if (m ~= n)
    help umfpack_btf
    error ('umfpack_btf:  A must be square') ;
end
m1 = size (b,1) ;
if (m1 ~= n)
    help umfpack_btf
    error ('umfpack_btf:  b has the wrong dimensions') ;
end

if (nargin < 3)
    Control = umfpack2 ;
end

%-------------------------------------------------------------------------------
% find the block triangular form
%-------------------------------------------------------------------------------

[p,q,r] = dmperm (A) ;
nblocks = length (r) - 1 ;

%-------------------------------------------------------------------------------
% solve the system
%-------------------------------------------------------------------------------

if (nblocks == 1 | sprank (A) < n)					    %#ok

    %---------------------------------------------------------------------------
    % matrix is irreducible or structurally singular
    %---------------------------------------------------------------------------

    x = umfpack_solve (A, '\', b, Control) ;

else

    %---------------------------------------------------------------------------
    % A (p,q) is in block triangular form
    %---------------------------------------------------------------------------

    b = b (p,:) ;
    A = A (p,q) ;
    x = zeros (size (b)) ;

    %---------------------------------------------------------------------------
    % merge adjacent singletons into a single upper triangular block
    %---------------------------------------------------------------------------

    [r, nblocks, is_triangular] = merge_singletons (r) ;

    %---------------------------------------------------------------------------
    % solve the system: x (q) = A\b
    %---------------------------------------------------------------------------

    for k = nblocks:-1:1

	% get the kth block
        k1 = r (k) ;
        k2 = r (k+1) - 1 ;

	% solve the system
	x (k1:k2,:) = solver (A (k1:k2, k1:k2), b (k1:k2,:), ...
	    is_triangular (k), Control) ;

        % off-diagonal block back substitution
        b (1:k1-1,:) = b (1:k1-1,:) - A (1:k1-1, k1:k2) * x (k1:k2,:) ;

    end

    x (q,:) = x ;

end

%-------------------------------------------------------------------------------
% merge_singletons
%-------------------------------------------------------------------------------

function [r, nblocks, is_triangular] = merge_singletons (r)
%
% Given r from [p,q,r] = dmperm (A), where A is square, return a modified r that
% reflects the merger of adjacent singletons into a single upper triangular
% block.  is_triangular (k) is 1 if the kth block is upper triangular.  nblocks
% is the number of new blocks.

nblocks = length (r) - 1 ;
bsize = r (2:nblocks+1) - r (1:nblocks) ;
t = [0 (bsize == 1)] ;
z = (t (1:nblocks) == 0 & t (2:nblocks+1) == 1) | t (2:nblocks+1) == 0 ;
y = [(find (z)) nblocks+1] ;
r = r (y) ;
nblocks = length (y) - 1 ;
is_triangular = y (2:nblocks+1) - y (1:nblocks) > 1 ;

%-------------------------------------------------------------------------------
% solve Ax=b, but check for small and/or triangular systems
%-------------------------------------------------------------------------------

function x = solver (A, b, is_triangular, Control)
if (is_triangular)
    % back substitution only
    x = A \ b ;
elseif (size (A,1) < 4)
    % a very small matrix, solve it as a dense linear system
    x = full (A) \ b ;
else
    % solve it as a sparse linear system
    x = umfpack_solve (A, '\', b, Control) ;
end
