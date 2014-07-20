function amd_make
% AMD_MAKE:  compiles the AMD mexFunction for MATLAB
%
% --------------------------------------------------------------------------
% AMD Version 1.1 (Jan. 21, 2004), Copyright (c) 2004 by Timothy A. Davis,
% Patrick R. Amestoy, and Iain S. Duff.  See ../README for License.
% email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.
% web: http://www.cise.ufl.edu/research/sparse/amd
% --------------------------------------------------------------------------
%
% See also: amd, amd_demo

help amd_make
fprintf ('Compiling the AMD mexFunction:\n') ;
cmd = sprintf ('mex -inline -O -output amd -I..%sInclude amd_mex.c', filesep) ;
files = {'amd_order', 'amd_dump', 'amd_postorder', 'amd_post_tree', ...
    'amd_aat', 'amd_2', 'amd_1', 'amd_defaults', 'amd_control', 'amd_info', ...
    'amd_valid' } ;
for i = 1 : length (files)
    cmd = sprintf ('%s ..%sSource%s%s.c', cmd, filesep, filesep, files {i}) ;
end
fprintf ('%s\n', cmd) ;
try
    eval (cmd) ;
catch
    fprintf ('Compilation not successful.\n') ;
end

input ('\nHit enter to run the AMD demo\n') ;
more on
amd_demo
more off
