function f = luflop (L, U)						    %#ok
%LUFLOP given L and U, computes # of flops required to compute them
%
% Example:
% f = luflop (L, U)
%
% Given an LU factorization, compute how many flops took to compute it.  This
% is the same as (assuming U has a zero-free diagonal):
%
%   Lnz = full (sum (spones (L))) - 1 ;
%   Unz = full (sum (spones (U')))' - 1 ;
%   f = 2*Lnz*Unz + sum (Lnz) ;
%
% except that no extra workspace is allocated for spones (L) and spones (U).
% L and U must be sparse.
%
% Note: the above expression has a subtle undercount when exact numerical
% cancelation occurs.  Try [L,U,P] = lu (sparse (ones (10))) and then
% luflop (L,U).
%
% See also LU

% Copyright 1995-2007 by Timothy A. Davis.

help luflop
error ('luflop mexFunction not found!  Use umfpack_make to compile luflop.') ;
