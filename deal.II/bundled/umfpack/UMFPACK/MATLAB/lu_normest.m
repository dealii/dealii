function rho = lu_normest (A, L, U)
%LU_NORMEST estimates norm (L*U-A, 1) without forming L*U-A
%
% Example:
%
%       rho = lu_normest (A, L, U)
%
% which estimates the computation of the 1-norm:
%
%       rho = norm (A-L*U, 1)
%
% Authors:  William W. Hager, Math Dept., Univ. of Florida
%       Timothy A. Davis, CISE Dept., Univ. of Florida
%       Gainesville, FL, 32611, USA.
%       based on normest1, contributed on November, 1997
%
% This code can be quite easily adapted to estimate the 1-norm of any
% matrix E, where E itself is dense or not explicitly represented, but the
% computation of E (and E') times a vector is easy.  In this case, our matrix
% of interest is:
%
%       E = A-L*U
%
% That is, L*U is the LU factorization of A, where A, L and U
% are sparse.  This code works for dense matrices A and L too,
% but it would not be needed in that case, since E is easy to compute
% explicitly.  For sparse A, L, and U, computing E explicitly would be quite
% expensive, and thus normest (A-L*U) would be prohibitive.
%
% For a detailed description, see Davis, T. A. and Hager, W. W.,
% Modifying a sparse Cholesky factorization, SIAM J. Matrix Analysis and
% Applications, 1999, vol. 20, no. 3, 606-627.
%
% See also normest

% The three places that the matrix-vector multiply E*x is used are highlighted.
% Note that E is never formed explicity.

% Copyright 1995-2007 by William W. Hager and Timothy A. Davis

[m n] = size (A) ;

if (m ~= n)
    % pad A, L, and U with zeros so that they are all square
    if (m < n)
        U = [ U ; (sparse (n-m,n)) ] ;
        L = [ L , (sparse (m,n-m)) ; (sparse (n-m,n)) ] ;
        A = [ A ; (sparse (n-m,n)) ] ;
    else
        U = [ U , (sparse (n,m-n)) ; (sparse (m-n,m)) ] ;
        L = [ L , (sparse (m,m-n)) ] ;
        A = [ A , (sparse (m,m-n)) ] ;
    end
end

[m n] = size (A) ;	    %#ok

notvisited = ones (m, 1) ;  % nonvisited(j) is zero if j is visited, 1 otherwise
rho = 0 ;    % the global rho

for trial = 1:3 % {

   x = notvisited ./ sum (notvisited) ;
   rho1 = 0 ;    % the current rho for this trial

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%% COMPUTE Ex1 = E*x EFFICIENTLY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Ex1 = (A*x) - L*(U*x) ;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   rho2 = norm (Ex1, 1) ;

   while rho2 > rho1 % {

        rho1 = rho2 ;
        y = 2*(Ex1 >= 0) - 1 ;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPUTE z = E'*y EFFICIENTLY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        z = (A'*y) - U'*(L'*y) ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [zj, j] = max (abs (z .* notvisited)) ;
        j = j (1) ;
        if (abs (z (j)) > z'*x) % {
            x = zeros (m, 1) ;
            x (j) = 1 ;
            notvisited (j) = 0 ;
        else % } {
            break ;
        end % }

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPUTE Ex1 = E*x EFFICIENTLY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ex1 = (A*x) - L*(U*x) ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rho2 = norm (Ex1, 1) ;

    end % }

    rho = max (rho, rho1) ;

end % }
