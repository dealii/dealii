function umfpack_report (Control, Info)
%UMFPACK_REPORT prints optional control settings and statistics
%
%   Example:
%       umfpack_report (Control, Info) ;
%
% Prints the current Control settings for umfpack2, and the statistical
% information returned by umfpack2 in the Info array.  If Control is
% an empty matrix, then the default control settings are printed.
%
% Control is 20-by-1, and Info is 90-by-1.  Not all entries are used.
%
% Alternative usages:
%
%       umfpack_report ([ ], Info) ;    print the default control parameters
%                                       and the Info array.
%       umfpack_report (Control) ;      print the control parameters only.
%       umfpack_report ;                print the default control parameters
%                                       and an empty Info array.
%
% See also umfpack, umfpack2, umfpack_make, umfpack_details,
% umfpack_demo, and umfpack_simple.

% Copyright 1995-2007 by Timothy A. Davis.

%-------------------------------------------------------------------------------
% get inputs, use defaults if input arguments not present
%-------------------------------------------------------------------------------

% The contents of Control and Info are defined in umfpack.h
if (nargin < 1)
    Control = [] ;
end
if (nargin < 2)
    Info = [] ;
end
if (isempty (Control))
    Control = umfpack2 ;
end
if (isempty (Info))
    Info = [ 0 (-ones (1, 89)) ] ;
end

%-------------------------------------------------------------------------------
% control settings
%-------------------------------------------------------------------------------

fprintf ('\nUMFPACK:  Control settings:\n\n') ;
fprintf ('    Control (1): print level: %d\n', Control (1)) ;
fprintf ('    Control (2): dense row parameter:    %g\n', Control (2)) ;
fprintf ('       "dense" rows have    > max (16, (%g)*16*sqrt(n_col)) entries\n', Control (2)) ;
fprintf ('    Control (3): dense column parameter: %g\n', Control (3)) ;
fprintf ('       "dense" columns have > max (16, (%g)*16*sqrt(n_row)) entries\n', Control (3)) ;
fprintf ('    Control (4): pivot tolerance: %g\n', Control (4)) ;
fprintf ('    Control (5): max block size for dense matrix kernels: %d\n', Control (5)) ;
prstrat ('    Control (6): strategy: %g ', Control (6)) ;
fprintf ('    Control (7): initial allocation ratio: %g\n', Control (7)) ;
fprintf ('    Control (8): max iterative refinement steps: %d\n', Control (8)) ;
fprintf ('    Control (13): 2-by-2 pivot tolerance: %g\n', Control (13)) ;
fprintf ('    Control (14): Q fixed during numeric factorization: %g ', Control (14)) ;
if (Control (14) > 0)
    fprintf ('(yes)\n') ;
elseif (Control (14) < 0)
    fprintf ('(no)\n') ;
else
    fprintf ('(auto)\n') ;
end
fprintf ('    Control (15): AMD dense row/column parameter: %g\n', Control (15)) ;
fprintf ('       "dense" rows/columns in A+A'' have > max (16, (%g)*sqrt(n)) entries.\n', Control (15)) ;
fprintf ('        Only used if the AMD ordering is used.\n') ;
fprintf ('    Control (16): diagonal pivot tolerance: %g\n', Control (16)) ;
fprintf ('        Only used if diagonal pivoting is attempted.\n') ;

fprintf ('    Control (17): scaling option: %g ', Control (17)) ;
if (Control (17) == 0)
    fprintf ('(none)\n') ;
elseif (Control (17) == 2)
    fprintf ('(scale the matrix by\n') ;
    fprintf ('        dividing each row by max. abs. value in each row)\n') ;
else
    fprintf ('(scale the matrix by\n') ;
    fprintf ('        dividing each row by sum of abs. values in each row)\n') ;
end

fprintf ('    Control (18): frontal matrix allocation ratio: %g\n', Control (18)) ;
fprintf ('    Control (19): drop tolerance: %g\n', Control (19)) ;
fprintf ('    Control (20): AMD and COLAMD aggressive absorption: %g ', Control (20)) ;
yes_no (Control (20)) ;

% compile-time options:

fprintf ('\n  The following options can only be changed at compile-time:\n') ;

if (Control (9) == 1)
    fprintf ('    Control (9): compiled to use the BLAS\n') ;
else
    fprintf ('    Control (9): compiled without the BLAS\n') ;
    fprintf ('        (you will not get the best possible performance)\n') ;
end

if (Control (10) == 1)
    fprintf ('    Control (10): compiled for MATLAB\n') ;
elseif (Control (10) == 2)
    fprintf ('    Control (10): compiled for MATLAB\n') ;
else
    fprintf ('    Control (10): not compiled for MATLAB\n') ;
    fprintf ('        Printing will be in terms of 0-based matrix indexing,\n') ;
    fprintf ('        not 1-based as is expected in MATLAB.  Diary output may\n') ;
    fprintf ('        not be properly recorded.\n') ;
end

if (Control (11) == 2)
    fprintf ('    Control (11): uses POSIX times ( ) to get CPU time and wallclock time.\n') ;
elseif (Control (11) == 1)
    fprintf ('    Control (11): uses getrusage to get CPU time.\n') ;
else
    fprintf ('    Control (11): uses ANSI C clock to get CPU time.\n') ;
    fprintf ('        The CPU time may wrap around, type "help cputime".\n') ;
end

if (Control (12) == 1)
    fprintf ('    Control (12): compiled with debugging enabled\n') ;
    fprintf ('        ###########################################\n') ;
    fprintf ('        ### This will be exceedingly slow! ########\n') ;
    fprintf ('        ###########################################\n') ;
else
    fprintf ('    Control (12): compiled for normal operation (no debugging)\n') ;
end

%-------------------------------------------------------------------------------
% Info:
%-------------------------------------------------------------------------------

if (nargin == 1)
    return
end

status = Info (1) ;
fprintf ('\nUMFPACK status:  Info (1): %d, ', status) ;

if (status == 0)
    fprintf ('OK\n') ;
elseif (status == 1)
    fprintf ('WARNING  matrix is singular\n') ;
elseif (status == -1)
    fprintf ('ERROR    out of memory\n') ;
elseif (status == -3)
    fprintf ('ERROR    numeric LU factorization is invalid\n') ;
elseif (status == -4)
    fprintf ('ERROR    symbolic LU factorization is invalid\n') ;
elseif (status == -5)
    fprintf ('ERROR    required argument is missing\n') ;
elseif (status == -6)
    fprintf ('ERROR    n <= 0\n') ;
elseif (status <= -7 & status >= -12 | status == -14)			    %#ok
    fprintf ('ERROR    matrix A is corrupted\n') ;
elseif (status == -13)
    fprintf ('ERROR    invalid system\n') ;
elseif (status == -15)
    fprintf ('ERROR    invalid permutation\n') ;
elseif (status == -911)
    fprintf ('ERROR    internal error!\n') ;
    fprintf ('Please report this error to Tim Davis (davis@cise.ufl.edu)\n') ;
else
    fprintf ('ERROR    unrecognized error.  Info array corrupted\n') ;
end

fprintf ('    (a -1 means the entry has not been computed):\n') ;

fprintf ('\n  Basic statistics:\n') ;
fprintf ('    Info (2):  %d, # of rows of A\n', Info (2)) ;
fprintf ('    Info (17): %d, # of columns of A\n', Info (17)) ;
fprintf ('    Info (3): %d, nnz (A)\n', Info (3)) ;
fprintf ('    Info (4): %d, Unit size, in bytes, for memory usage reported below\n', Info (4)) ;
fprintf ('    Info (5): %d, size of int (in bytes)\n', Info (5)) ;
fprintf ('    Info (6): %d, size of UF_long (in bytes)\n', Info (6)) ;
fprintf ('    Info (7): %d, size of pointer (in bytes)\n', Info (7)) ;
fprintf ('    Info (8): %d, size of numerical entry (in bytes)\n', Info (8)) ;

fprintf ('\n  Pivots with zero Markowitz cost removed to obtain submatrix S:\n') ;
fprintf ('    Info (57): %d, # of pivots with one entry in pivot column\n', Info (57)) ;
fprintf ('    Info (58): %d, # of pivots with one entry in pivot row\n', Info (58)) ;
fprintf ('    Info (59): %d, # of rows/columns in submatrix S (if square)\n', Info (59)) ;
fprintf ('    Info (60): ') ;
if (Info (60) > 0)
    fprintf ('submatrix S square and diagonal preserved\n') ;
elseif (Info  (60) == 0)
    fprintf ('submatrix S not square or diagonal not preserved\n') ;
else
    fprintf ('\n') ;
end
fprintf ('    Info (9):  %d, # of "dense" rows in S\n', Info (9)) ;
fprintf ('    Info (10): %d, # of empty rows in S\n', Info (10)) ;
fprintf ('    Info (11): %d, # of "dense" columns in S\n', Info (11)) ;
fprintf ('    Info (12): %d, # of empty columns in S\n', Info (12)) ;
fprintf ('    Info (34): %g, symmetry of pattern of S\n', Info (34)) ;
fprintf ('    Info (35): %d, # of off-diagonal nonzeros in S+S''\n', Info (35)) ;
fprintf ('    Info (36): %d, nnz (diag (S))\n', Info (36)) ;

fprintf ('\n  2-by-2 pivoting to place large entries on diagonal:\n') ;
fprintf ('    Info (52): %d, # of small diagonal entries of S\n', Info (52)) ;
fprintf ('    Info (53): %d, # of unmatched small diagonal entries\n', Info (53)) ;
fprintf ('    Info (54): %g, symmetry of P2*S\n', Info (54)) ;
fprintf ('    Info (55): %d, # of off-diagonal entries in (P2*S)+(P2*S)''\n', Info (55)) ;
fprintf ('    Info (56): %d, nnz (diag (P2*S))\n', Info (56)) ;

fprintf ('\n  AMD results, for strict diagonal pivoting:\n') ;
fprintf ('    Info (37): %d, est. nz in L and U\n', Info (37)) ;
fprintf ('    Info (38): %g, est. flop count\n', Info (38)) ;
fprintf ('    Info (39): %g, # of "dense" rows in S+S''\n', Info (39)) ;
fprintf ('    Info (40): %g, est. max. nz in any column of L\n', Info (40)) ;

fprintf ('\n  Final strategy selection, based on the analysis above:\n') ;
prstrat ('    Info (19): %d, strategy used ', Info (19)) ;
fprintf ('    Info (20): %d, ordering used ', Info (20)) ;
if (Info (20) == 0)
    fprintf ('(COLAMD on A)\n') ;
elseif (Info (20) == 1)
    fprintf ('(AMD on A+A'')\n') ;
elseif (Info (20) == 2)
    fprintf ('(provided by user)\n') ;
else
    fprintf ('(undefined ordering option)\n') ;
end
fprintf ('    Info (32): %d, Q fixed during numeric factorization: ', Info (32)) ;
yes_no (Info (32)) ;
fprintf ('    Info (33): %d, prefer diagonal pivoting: ', Info (33)) ;
yes_no (Info (33)) ;

fprintf ('\n  symbolic analysis time and memory usage:\n') ;
fprintf ('    Info (13): %d, defragmentations during symbolic analysis\n', Info (13)) ;
fprintf ('    Info (14): %d, memory used during symbolic analysis (Units)\n', Info (14)) ;
fprintf ('    Info (15): %d, final size of symbolic factors (Units)\n', Info (15)) ;
fprintf ('    Info (16): %.2f, symbolic analysis CPU time (seconds)\n', Info (16)) ;
fprintf ('    Info (18): %.2f, symbolic analysis wall clock time (seconds)\n', Info (18)) ;

fprintf ('\n  Estimates computed in the symbolic analysis:\n') ;
fprintf ('    Info (21): %d, est. size of LU factors (Units)\n', Info (21)) ;
fprintf ('    Info (22): %d, est. total peak memory usage (Units)\n', Info (22)) ;
fprintf ('    Info (23): %d, est. factorization flop count\n', Info (23)) ;
fprintf ('    Info (24): %d, est. nnz (L)\n', Info (24)) ;
fprintf ('    Info (25): %d, est. nnz (U)\n', Info (25)) ;
fprintf ('    Info (26): %d, est. initial size, variable-part of LU (Units)\n', Info (26)) ;
fprintf ('    Info (27): %d, est. peak size, of variable-part of LU (Units)\n', Info (27)) ;
fprintf ('    Info (28): %d, est. final size, of variable-part of LU (Units)\n', Info (28)) ;
fprintf ('    Info (29): %d, est. max frontal matrix size (# of entries)\n', Info (29)) ;
fprintf ('    Info (30): %d, est. max # of rows in frontal matrix\n', Info (30)) ;
fprintf ('    Info (31): %d, est. max # of columns in frontal matrix\n', Info (31)) ;

fprintf ('\n  Computed in the numeric factorization (estimates shown above):\n') ;
fprintf ('    Info (41): %d, size of LU factors (Units)\n', Info (41)) ;
fprintf ('    Info (42): %d, total peak memory usage (Units)\n', Info (42)) ;
fprintf ('    Info (43): %d, factorization flop count\n', Info (43)) ;
fprintf ('    Info (44): %d, nnz (L)\n', Info (44)) ;
fprintf ('    Info (45): %d, nnz (U)\n', Info (45)) ;
fprintf ('    Info (46): %d, initial size of variable-part of LU (Units)\n', Info (46)) ;
fprintf ('    Info (47): %d, peak size of variable-part of LU (Units)\n', Info (47)) ;
fprintf ('    Info (48): %d, final size of variable-part of LU (Units)\n', Info (48)) ;
fprintf ('    Info (49): %d, max frontal matrix size (# of numerical entries)\n', Info (49)) ;
fprintf ('    Info (50): %d, max # of rows in frontal matrix\n', Info (50)) ;
fprintf ('    Info (51): %d, max # of columns in frontal matrix\n', Info (51)) ;

fprintf ('\n  Computed in the numeric factorization (no estimates computed a priori):\n') ;
fprintf ('    Info (61): %d, defragmentations during numeric factorization\n', Info (61)) ;
fprintf ('    Info (62): %d, reallocations during numeric factorization\n', Info (62)) ;
fprintf ('    Info (63): %d, costly reallocations during numeric factorization\n', Info (63)) ;
fprintf ('    Info (64): %d, integer indices in compressed pattern of L and U\n', Info (64)) ;
fprintf ('    Info (65): %d, numerical values stored in L and U\n', Info (65)) ;
fprintf ('    Info (66): %.2f, numeric factorization CPU time (seconds)\n', Info (66)) ;
fprintf ('    Info (76): %.2f, numeric factorization wall clock time (seconds)\n', Info (76)) ;
if (Info (66) > 0.05 & Info (43) > 0)					    %#ok
fprintf ('    mflops in numeric factorization phase: %.2f\n', 1e-6 * Info (43) / Info (66)) ;
end
fprintf ('    Info (67): %d, nnz (diag (U))\n', Info (67)) ;
fprintf ('    Info (68): %g, reciprocal condition number estimate\n', Info (68)) ;
fprintf ('    Info (69): %g, matrix was ', Info (69)) ;
if (Info (69) == 0)
    fprintf ('not scaled\n') ;
elseif (Info (69) == 2)
    fprintf ('scaled (row max)\n') ;
else
    fprintf ('scaled (row sum)\n') ;
end
fprintf ('    Info (70): %g, min. scale factor of rows of A\n', Info (70)) ;
fprintf ('    Info (71): %g, max. scale factor of rows of A\n', Info (71)) ;
fprintf ('    Info (72): %g, min. abs. on diagonal of U\n', Info (72)) ;
fprintf ('    Info (73): %g, max. abs. on diagonal of U\n', Info (73)) ;
fprintf ('    Info (74): %g, initial allocation parameter used\n', Info (74)) ;
fprintf ('    Info (75): %g, # of forced updates due to frontal growth\n', Info (75)) ;
fprintf ('    Info (77): %d, # of off-diaogonal pivots\n', Info (77)) ;
fprintf ('    Info (78): %d, nnz (L), if no small entries dropped\n', Info (78)) ;
fprintf ('    Info (79): %d, nnz (U), if no small entries dropped\n', Info (79)) ;
fprintf ('    Info (80): %d, # of small entries dropped\n', Info (80)) ;

fprintf ('\n  Computed in the solve step:\n') ;
fprintf ('    Info (81): %d, iterative refinement steps taken\n', Info (81)) ;
fprintf ('    Info (82): %d, iterative refinement steps attempted\n', Info (82)) ;
fprintf ('    Info (83): %g, omega(1), sparse-backward error estimate\n', Info (83)) ;
fprintf ('    Info (84): %g, omega(2), sparse-backward error estimate\n', Info (84)) ;
fprintf ('    Info (85): %d, solve flop count\n', Info (85)) ;
fprintf ('    Info (86): %.2f, solve CPU time (seconds)\n', Info (86)) ;
fprintf ('    Info (87): %.2f, solve wall clock time (seconds)\n', Info (87)) ;

fprintf ('\n    Info (88:90): unused\n\n') ;

%-------------------------------------------------------------------------------

function prstrat (fmt, strategy)
% prstrat print the ordering strategy
fprintf (fmt, strategy) ;
if (strategy == 1)
    fprintf ('(unsymmetric)\n') ;
    fprintf ('        Q = COLAMD (A), Q refined during numerical\n') ;
    fprintf ('        factorization, and no attempt at diagonal pivoting.\n') ;
elseif (strategy == 2)
    fprintf ('(symmetric, with 2-by-2 pivoting)\n') ;
    fprintf ('        P2 = row permutation to place large values on the diagonal\n') ;
    fprintf ('        Q = AMD (P2*A+(P2*A)''), Q not refined during numeric factorization,\n') ;
    fprintf ('        and diagonal pivoting attempted.\n') ;
elseif (strategy == 3)
    fprintf ('(symmetric)\n') ;
    fprintf ('        Q = AMD (A+A''), Q not refined during numeric factorization,\n') ;
    fprintf ('        and diagonal pivoting (P=Q'') attempted.\n') ;
else
    % strategy = 0 ;
    fprintf ('(auto)\n') ;
end

%-------------------------------------------------------------------------------

function yes_no (s)
% yes_no print yes or no
if (s == 0)
    fprintf ('(no)\n') ;
else
    fprintf ('(yes)\n') ;
end
