function umfpack_make (lapack)
%UMFPACK_MAKE to compile umfpack2 for use in MATLAB
%
% Compiles the umfpack2 mexFunction and then runs a simple demo.
%
% Example:
%   umfpack_make				% use default LAPACK and BLAS
%   umfpack_make ('lcc_lib/libmwlapack.lib')	% try this if umfpack_make fails
%
% See also: umfpack, umfpack2, umfpack_details, umfpack_report, umfpack_demo,
% and umfpack_simple.

% Copyright 1995-2007 by Timothy A. Davis.

details = 0 ;

d = '' ;
if (~isempty (strfind (computer, '64')))
    d = ' -largeArrayDims' ;
end

[v,pc] = getversion ;
fprintf ('Compiling UMFPACK for MATLAB Version %g\n', v) ;

if (pc)
    obj = 'obj' ;
else
    obj = 'o' ;
end

kk = 0 ;

%-------------------------------------------------------------------------------
% BLAS option
%-------------------------------------------------------------------------------

if (nargin < 1)
    if (pc)
	if (v < 6.5)
	    % MATLAB 6.1 and earlier: use the version supplied here
	    lapack = 'lcc_lib/libmwlapack.lib' ;
	    fprintf ('Using %s.  If this fails with dgemm and others\n',lapack);
	    fprintf ('undefined, then edit umfpack_make.m and modify the') ;
	    fprintf (' statement:\nlapack = ''%s'' ;\n', lapack) ;
	else
	    lapack = 'libmwlapack.lib' ;
	end
    else
	% For other systems, mex should find lapack on its own, but this has
	% been broken in MATLAB R2007a; the following is now required.
	lapack = '-lmwlapack' ;
    end
end

%-------------------------------------------------------------------------------
% -DNPOSIX option (for sysconf and times timer routines)
%-------------------------------------------------------------------------------

posix = '' ;

% if (~pc)
%     msg = [ ...
%    '--------------------------------------------------------------\n', ...
%    '\nUMFPACK can use the POSIX routines sysconf () and times ()\n', ...
%    'to provide CPU time and wallclock time statistics.  If you do not\n', ...
%    'have a POSIX-compliant operating system, then UMFPACK won''t\n', ...
%    'compile.  If you don''t know which option to pick, try the\n', ...
%    'default.  If you get an error saying that sysconf and/or times\n', ...
%    'are not defined, then recompile with the non-POSIX option.\n', ...
%    '\nPlease select one of the following options:\n', ...
%    '    1:  use POSIX sysconf and times routines (default)\n', ...
%    '    2:  do not use POSIX routines\n'] ;
%    fprintf (msg) ;
%    posix = str2num (input (': ', 's')) ;
%    if (isempty (posix))
%	posix = 1 ;
%    end
%    if (posix == 2)
%        fprintf ('\nNot using POSIX sysconf and times routines.\n') ;
%        posix = ' -DNPOSIX' ;
%    else
%        fprintf ('\nUsing POSIX sysconf and times routines.\n') ;
%        posix = '' ;
%    end
% end

%-------------------------------------------------------------------------------
% mex command
%-------------------------------------------------------------------------------

umfdir = '../Source/' ;
amddir = '../../AMD/Source/' ;
incdir = ' -I../Include -I../Source -I../../AMD/Include -I../../UFconfig' ;
mx = sprintf ('mex -O%s%s%s ', posix, incdir, d) ;
% fprintf ('compile options:\n%s\n', mx) ;

%-------------------------------------------------------------------------------
% source files
%-------------------------------------------------------------------------------

% non-user-callable umf_*.[ch] files:
umfch = { 'assemble', 'blas3_update', ...
        'build_tuples', 'create_element', ...
        'dump', 'extend_front', 'garbage_collection', ...
        'get_memory', 'init_front', 'kernel', ...
        'kernel_init', 'kernel_wrapup', ...
        'local_search', 'lsolve', 'ltsolve', ...
        'mem_alloc_element', 'mem_alloc_head_block', ...
        'mem_alloc_tail_block', 'mem_free_tail_block', ...
        'mem_init_memoryspace', ...
        'report_vector', 'row_search', 'scale_column', ...
        'set_stats', 'solve', 'symbolic_usage', 'transpose', ...
        'tuple_lengths', 'usolve', 'utsolve', 'valid_numeric', ...
        'valid_symbolic', 'grow_front', 'start_front', '2by2', ...
	'store_lu', 'scale' } ;

% non-user-callable umf_*.[ch] files, int versions only (no real/complex):
umfint = { 'analyze', 'apply_order', 'colamd', 'free', 'fsize', ...
        'is_permutation', 'malloc', 'realloc', 'report_perm', ...
	'singletons' } ;

% non-user-callable and user-callable amd_*.[ch] files (int versions only):
amdsrc = { 'aat', '1', '2', 'dump', 'postorder', 'post_tree', 'defaults', ...
        'order', 'control', 'info', 'valid', 'preprocess', 'global' } ;

% user-callable umfpack_*.[ch] files (real/complex):
user = { 'col_to_triplet', 'defaults', 'free_numeric', ...
        'free_symbolic', 'get_numeric', 'get_lunz', ...
        'get_symbolic', 'get_determinant', 'numeric', 'qsymbolic', ...
        'report_control', 'report_info', 'report_matrix', ...
        'report_numeric', 'report_perm', 'report_status', ...
        'report_symbolic', 'report_triplet', ...
        'report_vector', 'solve', 'symbolic', ...
        'transpose', 'triplet_to_col', 'scale' ...
	'load_numeric', 'save_numeric', 'load_symbolic', 'save_symbolic' } ;

% user-callable umfpack_*.[ch], only one version
generic = { 'timer', 'tictoc', 'global' } ;

M = cell (0) ;

%-------------------------------------------------------------------------------
% Create the umfpack2 and amd2 mexFunctions for MATLAB (int versions only)
%-------------------------------------------------------------------------------

for k = 1:length(umfint)
    [M, kk] = make (M, '%s -DDLONG -c %sumf_%s.c', 'umf_%s.%s', ...
	'umf_%s_%s.%s', mx, umfint {k}, umfint {k}, 'm', obj, umfdir, ...
	kk, details) ;
end

rules = { [mx ' -DDLONG'] , [mx ' -DZLONG'] } ;
kinds = { 'md', 'mz' } ;

for what = 1:2

    rule = rules {what} ;
    kind = kinds {what} ;

    [M, kk] = make (M, '%s -DCONJUGATE_SOLVE -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s.%s', rule, 'ltsolve', 'lhsolve', kind, obj, umfdir, ...
	kk, details) ;

    [M, kk] = make (M, '%s -DCONJUGATE_SOLVE -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s.%s', rule, 'utsolve', 'uhsolve', kind, obj, umfdir, ...
	kk, details) ;

    [M, kk] = make (M, '%s -DDO_MAP -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s_map_nox.%s', rule, 'triplet', 'triplet', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -DDO_VALUES -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s_nomap_x.%s', rule, 'triplet', 'triplet', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -c %sumf_%s.c', 'umf_%s.%s',  ...
        'umf_%s_%s_nomap_nox.%s', rule, 'triplet', 'triplet', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -DDO_MAP -DDO_VALUES -c %sumf_%s.c', 'umf_%s.%s', ...
        'umf_%s_%s_map_x.%s', rule, 'triplet', 'triplet', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -DFIXQ -c %sumf_%s.c', 'umf_%s.%s', ...
	'umf_%s_%s_fixq.%s', rule, 'assemble', 'assemble', kind, obj, ...
	umfdir, kk, details) ;

    [M, kk] = make (M, '%s -DDROP -c %sumf_%s.c', 'umf_%s.%s', ...
	'umf_%s_%s_drop.%s', rule, 'store_lu', 'store_lu', kind, obj, ...
	umfdir, kk, details) ;

    for k = 1:length(umfch)
        [M, kk] = make (M, '%s -c %sumf_%s.c', 'umf_%s.%s', 'umf_%s_%s.%s', ...
            rule, umfch {k}, umfch {k}, kind, obj, umfdir, kk, details) ;
    end

    [M, kk] = make (M, '%s -DWSOLVE -c %sumfpack_%s.c', 'umfpack_%s.%s', ...
        'umfpack_%s_w%s.%s', rule, 'solve', 'solve', kind, obj, umfdir, ...
	kk, details) ;

    for k = 1:length(user)
        [M, kk] = make (M, '%s -c %sumfpack_%s.c', 'umfpack_%s.%s', ...
            'umfpack_%s_%s.%s', rule, user {k}, user {k}, kind, obj, ...
	    umfdir, kk, details) ;
    end
end

for k = 1:length(generic)
    [M, kk] = make (M, '%s -c %sumfpack_%s.c', 'umfpack_%s.%s', ...
	'umfpack_%s_%s.%s', mx, generic {k}, generic {k}, 'm', obj, ...
	umfdir, kk, details) ;
end

%----------------------------------------
% AMD routines (int only)
%----------------------------------------

for k = 1:length(amdsrc)
    [M, kk] = make (M, '%s -DDLONG -c %samd_%s.c', 'amd_%s.%s', ...
	'amd_%s_%s.%s', mx, amdsrc {k}, amdsrc {k}, 'm', obj, amddir, ...
	kk, details) ;
end

%----------------------------------------
% compile the umfpack2 mexFunction
%----------------------------------------

C = sprintf ('%s -output umfpack2 umfpackmex.c', mx) ;
for i = 1:length (M)
    C = [C ' ' (M {i})] ;   %#ok
end
C = [C ' ' lapack] ;
kk = cmd (C, kk, details) ;

%----------------------------------------
% delete the object files
%----------------------------------------

for i = 1:length (M)
    rmfile (M {i}) ;
end

%----------------------------------------
% compile the luflop mexFunction
%----------------------------------------

cmd (sprintf ('%s -output luflop luflopmex.c', mx), kk, details) ;

fprintf ('\nUMFPACK successfully compiled\n') ;

%===============================================================================
% end of umfpack_make
%===============================================================================


%-------------------------------------------------------------------------------

function rmfile (file)
% rmfile:  delete a file, but only if it exists
if (length (dir (file)) > 0)						    %#ok
    delete (file) ;
end

%-------------------------------------------------------------------------------

function cpfile (src, dst)
% cpfile:  copy the src file to the filename dst, overwriting dst if it exists
rmfile (dst)
if (length (dir (src)) == 0)	%#ok
    fprintf ('File does not exist: %s\n', src) ;
    error ('File does not exist') ;
end
copyfile (src, dst) ;

%-------------------------------------------------------------------------------

function mvfile (src, dst)
% mvfile:  move the src file to the filename dst, overwriting dst if it exists
cpfile (src, dst) ;
rmfile (src) ;

%-------------------------------------------------------------------------------

function kk = cmd (s, kk, details)
%CMD: evaluate a command, and either print it or print a "."
if (details)
    fprintf ('%s\n', s) ;
else
    if (mod (kk, 60) == 0)
	fprintf ('\n') ;
    end
    kk = kk + 1 ;
    fprintf ('.') ;
end
eval (s) ;

%-------------------------------------------------------------------------------

function [M, kk] = make (M, s, src, dst, rule, file1, file2, kind, obj, ...
    srcdir, kk, details)
% make:  execute a "make" command for a source file
kk = cmd (sprintf (s, rule, srcdir, file1), kk, details) ;
src = sprintf (src, file1, obj) ;
dst = sprintf (dst, kind, file2, obj) ;
mvfile (src, dst) ;
M {end + 1} = dst ;


%-------------------------------------------------------------------------------
function [v,pc] = getversion
% determine the MATLAB version, and return it as a double.
% only the primary and secondary version numbers are kept.
% MATLAB 7.0.4 becomes 7.0, version 6.5.2 becomes 6.5, etc.
v = version ;
t = find (v == '.') ;
if (length (t) > 1)
    v = v (1:(t(2)-1)) ;
end
v = str2double (v) ;
try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
catch
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
end
