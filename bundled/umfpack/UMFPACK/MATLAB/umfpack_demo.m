function umfpack_demo (c)
%UMFPACK_DEMO a lenghty demo
%
% A demo of UMFPACK for MATLAB.
%
% Example:
%   umfpack_demo
%
% See also umfpack, umfpack2, umfpack_make, umfpack_details, umfpack_report,
% and umfpack_simple.

% Copyright 1995-2007 by Timothy A. Davis.

%-------------------------------------------------------------------------------
% get default control parameters
%-------------------------------------------------------------------------------

control = umfpack2 ;
if (nargin < 1)
    fprintf ('\nEnter the printing level for UMFPACK''s output statistics:\n') ;
    fprintf ('0: none, 1: errors only, 2: statistics, 4: print some outputs\n');
    c = input ('5: print all output [default is 1]: ', 's') ;
    c = str2double (c) ;
end
if (isempty (c))
    c = 1 ;
end
control (1) = c ;

%-------------------------------------------------------------------------------
% solve a simple system
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('Factor and solve a small system, Ax=b, using default parameters\n') ;
if (control (1) > 1)
    fprintf ('(except for verbose printing enabled)\n') ;
end

load west0067
A = spconvert (west0067) ;

n = size (A, 1) ;
b = rand (n, 1) ;

fprintf ('Solving Ax=b via UMFPACK:\n') ;
xu = umfpack2 (A, '\', b, control) ;

fprintf ('Solving Ax=b via MATLAB:\n') ;
xm = A\b ;

fprintf ('Difference between UMFPACK and MATLAB solution: %g\n', ...
    norm (xu - xm, Inf)) ;

%-------------------------------------------------------------------------------
% spy the results
%-------------------------------------------------------------------------------

figure (1)
clf

subplot (2,3,1)
spy (A)
title ('The matrix A') ;

subplot (2,3,2)
[P1, Q1, Fr, Ch, Info] = umfpack2 (A, 'symbolic') ;			    %#ok
treeplot (Fr (1:end-1,2)') ;
title ('Supernodal column elimination tree') ;

subplot (2,3,3)
spy (P1 * A * Q1)
title ('A, with initial row and column order') ;

subplot (2,3,4)
fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('\nFactorizing [L, U, P, Q, R] = umfpack2 (A)\n') ;
[L, U, P, Q, R] = umfpack2 (A) ;
spy (P*A*Q)
title ('A, with final row/column order') ;

fprintf ('\nP * (R\\A) * Q - L*U should be zero:\n') ;
fprintf ('norm (P*(R\\A)*Q - L*U, 1) = %g (exact) %g (estimated)\n', ...
    norm (P * (R\A) * Q - L*U, 1), lu_normest (P * (R\A) * Q,  L, U)) ;

fprintf ('\nSolution to Ax=b via UMFPACK factorization:\n') ;
fprintf ('x = Q * (U \\ (L \\ (P * (R \\ b))))\n') ;
xu = Q * (U \ (L \ (P * (R \ b)))) ;

fprintf ('\nUMFPACK flop count: %d\n', luflop (L, U)) ;

subplot (2,3,5)
spy (spones (L) + spones (U))
title ('UMFPACK LU factors') ;

subplot (2,3,6)
fprintf ('\nFactorizing [L, U, P] = lu (A (:, q))\n') ;
fprintf ('If you are using a version of MATLAB prior to V6.0, then the\n') ;
fprintf ('following statement (q = colamd (A)) may fail.  Either download\n');
fprintf ('colamd from http://www.cise.ufl.edu/research/sparse, upgrade to\n') ;
fprintf ('MATLAB V6.0 or later, or replace the statement with\n') ;
fprintf ('q = colmmd (A) ;\n') ;
try
    q = colamd (A) ;
catch
    fprintf ('\n *** colamd not found, using colmmd instead *** \n') ;
    q = colmmd (A) ;
end
[L, U, P] = lu (A (:,q)) ;
spy (spones (L) + spones (U))
title ('MATLAB LU factors') ;

fprintf ('\nSolution to Ax=b via MATLAB factorization:\n') ;
fprintf ('x = U \\ (L \\ (P * b)) ; x (q) = x ;\n') ;
xm = U \ (L \ (P * b)) ;
xm (q) = xm ;

fprintf ('Difference between UMFPACK and MATLAB solution: %g\n', ...
    norm (xu - xm, Inf)) ;

fprintf ('\nMATLAB LU flop count: %d\n', luflop (L, U)) ;

%-------------------------------------------------------------------------------
% solve A'x=b
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('Solve A''x=b:\n') ;

fprintf ('Solving A''x=b via UMFPACK:\n') ;
xu = umfpack2 (b', '/', A, control) ;
xu = xu' ;

fprintf ('Solving A''x=b via MATLAB:\n') ;
xm = (b'/A)' ;

fprintf ('Difference between UMFPACK and MATLAB solution: %g\n', ...
    norm (xu - xm, Inf)) ;

%-------------------------------------------------------------------------------
% factor A' and then solve Ax=b using the factors of A'
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('Compute C = A'', and compute the LU factorization of C.\n') ;
fprintf ('Factorizing A'' can sometimes be better than factorizing A itself\n');
fprintf ('(less work and memory usage).  Solve C''x=b; the solution is the\n') ;
fprintf ('same as the solution to Ax=b for the original A.\n');

C = A' ;

% factorize C (P,Q) = L*U
[L, U, P, Q, R, info] = umfpack2 (C, control) ;				    %#ok

fprintf ('\nP * (R\\C) * Q - L*U should be zero:\n') ;
fprintf ('norm (P*(R\\C)*Q - L*U, 1) = %g (exact) %g (estimated)\n', ...
    norm (P * (R\C) * Q - L*U, 1), lu_normest (P * (R\C) * Q,  L, U)) ;

fprintf ('\nSolution to Ax=b via UMFPACK, using the factors of C:\n') ;
fprintf ('x = R \\ (P'' * (L'' \\ (U'' \\ (Q'' * b)))) ;\n') ;
xu = R \ (P' * (L' \ (U' \ (Q' * b)))) ;

fprintf ('Solution to Ax=b via MATLAB:\n') ;
xm = A\b ;

fprintf ('Difference between UMFPACK and MATLAB solution: %g\n', ...
    norm (xu - xm, Inf)) ;

%-------------------------------------------------------------------------------
% solve Ax=B
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('\nSolve AX=B, where B is n-by-10, and sparse\n') ;
B = sprandn (n, 10, 0.05) ;
XU = umfpack_solve (A, '\', B, control) ;
XM = A\B ;

fprintf ('Difference between UMFPACK and MATLAB solution: %g\n', ...
    norm (XU - XM, Inf)) ;

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('\nSolve AX=B, where B is n-by-10, and sparse, using umfpack_btf\n') ;
XU = umfpack_btf (A, B, control) ;

fprintf ('Difference between UMFPACK and MATLAB solution: %g\n', ...
    norm (XU - XM, Inf)) ;

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('\nSolve A''X=B, where B is n-by-10, and sparse\n') ;
XU = umfpack_solve (B', '/', A, control) ;
XM = B'/A ;

fprintf ('Difference between UMFPACK and MATLAB solution: %g\n', ...
    norm (XU - XM, Inf)) ;

%-------------------------------------------------------------------------------
% compute the determinant
%-------------------------------------------------------------------------------

fprintf ('\n--------------------------------------------------------------\n') ;
fprintf ('det(A): %g  UMFPACK determinant: %g\n', det (A), umfpack2 (A, 'det'));
