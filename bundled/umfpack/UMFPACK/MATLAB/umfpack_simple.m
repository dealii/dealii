%UMFPACK_SIMPLE a simple demo
%
% Example:
%   umfpack_simple
%
% Copyright 1995-2007 by Timothy A. Davis.
%
% UMFPACK License:
%
%     Your use or distribution of UMFPACK or any modified version of
%     UMFPACK implies that you agree to this License.  UMFPACK is
%     is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 2.1 of the License, or (at your option) any later version.

% Availability: http://www.cise.ufl.edu/research/sparse/umfpack
%
% See also: umfpack, umfpack2, umfpack_details

help umfpack_simple

format short

A = [
 2  3  0  0  0
 3  0  4  0  6
 0 -1 -3  2  0
 0  0  1  0  0
 0  4  2  0  1
] ;
fprintf ('A = \n') ; disp (A) ;

A = sparse (A) ;

b = [8 45 -3 3 19]' ;
fprintf ('b = \n') ; disp (b) ;

fprintf ('Solution to Ax=b via UMFPACK:\n') ;
fprintf ('x1 = umfpack2 (A, ''\\'', b)\n') ;

x1 = umfpack2 (A, '\', b) ;
fprintf ('x1 = \n') ; disp (x1) ;

fprintf ('Solution to Ax=b via MATLAB:\n') ;
fprintf ('x2 = A\\b\n') ;

x2 = A\b ;
fprintf ('x2 = \n') ; disp (x2) ;

fprintf ('norm (x1-x2) should be small: %g\n', norm (x1-x2)) ;

fprintf ('Type ''umfpack_demo'' for a full demo of UMFPACK\n') ;
