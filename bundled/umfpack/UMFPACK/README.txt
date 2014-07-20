UMFPACK Version 5.0.2:  a set of routines solving sparse linear systems via LU
    factorization.  Requires three other packages:  the BLAS (dense matrix
    operations), AMD (sparse matrix minimum degree ordering), and UFconfig.
    Includes a C-callable and MATLAB interface, and a basic FORTRAN 77
    interface to a subset of the C-callable routines.  Requires AMD Version
    2.0 or later.

The AMD, UFconfig, and UMFPACK directories must all reside in the same parent
directory.

Quick start (Unix, or Windows with Cygwin):

    To compile, test, and install both UMFPACK and AMD, the UMFPACK and AMD
    directories must be in the same parent directory.  To configure, edit
    the UFconfig/UFconfig.mk file (otherwise, you may get warnings that the
    BLAS (dgemm, etc) are not found).  You may use UMFPACK_CONFIG = -DNBLAS in
    the UFconfig/UFconfig.mk file, to avoid using the BLAS, but UMFPACK will be
    slow.  Next, cd to this directory (UMFPACK) and type "make".  To compile
    and run a FORTRAN demo program for Harwell/Boeing matrices, type "make hb".
    To compile a FORTRAN main program that calls the 32-bit C-callable UMFPACK
    library, type "make fortran".  When done, type "make clean" to remove
    unused *.o files (keeps the compiled libraries and demo programs).  See
    the User Guide (Doc/UserGuide.pdf), or ../UFconfig/UFconfig.mk for more
    details (including options for compiling in 64-bit mode).

Quick start (for MATLAB users):

    To compile, test, and install the UMFPACK mexFunction, cd to the
    UMFPACK/MATLAB directory and type umfpack_make at the MATLAB prompt.

    NOTE: DO NOT ATTEMPT TO USE THIS CODE IN 64-BIT MATLAB (v7.3).
    It is not yet ported to that version of MATLAB.

--------------------------------------------------------------------------------

UMFPACK, Copyright (c) 1995-2006 by Timothy A.  Davis.  All Rights Reserved.
UMFPACK is available under alternate licences; contact T. Davis for details.

UMFPACK License:

    Your use or distribution of UMFPACK or any modified version of
    UMFPACK implies that you agree to this License.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
    USA

    Permission is hereby granted to use or copy this program under the
    terms of the GNU LGPL, provided that the Copyright, this License,
    and the Availability of the original version is retained on all copies.
    User documentation of any code that uses this code or any modified
    version of this code must cite the Copyright, this License, the
    Availability note, and "Used by permission." Permission to modify
    the code and to distribute modified code is granted, provided the
    Copyright, this License, and the Availability note are retained,
    and a notice that the code was modified is included.

Availability:

    http://www.cise.ufl.edu/research/sparse/umfpack

    UMFPACK (including versions 2.2.1 and earlier, in FORTRAN) is available at
    http://www.cise.ufl.edu/research/sparse.  MA38 is available in the Harwell
    Subroutine Library.  This version of UMFPACK includes a modified form of
    COLAMD Version 2.0, originally released on Jan. 31, 2000, also available at
    http://www.cise.ufl.edu/research/sparse.  COLAMD V2.0 is also incorporated
    as a built-in function in MATLAB version 6.1, by The MathWorks, Inc.
    (http://www.mathworks.com).  COLAMD V1.0 appears as a column-preordering
    in SuperLU (SuperLU is available at http://www.netlib.org).
    UMFPACK v4.0 is a built-in routine in MATLAB 6.5.
    UMFPACK v4.3 is a built-in routine in MATLAB 7.1.

--------------------------------------------------------------------------------

Refer to ../AMD/README for the License for AMD, which is a separate
package for ordering sparse matrices that is required by UMFPACK.
UMFPACK v4.5 cannot use AMD v1.1 or earlier.  UMFPACK 5.0
requires AMD v2.0 or later.

--------------------------------------------------------------------------------

This is the UMFPACK README.txt file.  It is a terse overview of UMFPACK.
Refer to the User Guide (Doc/UserGuide.pdf) for how to install and use UMFPACK,
or to the Quick Start Guide, QuickStart.pdf.

Description:

    UMFPACK is a set of routines for solving unsymmetric sparse linear systems,
    Ax=b, using the Unsymmetric MultiFrontal method.  Written in ANSI/ISO C,
    with a MATLAB (Version 6.0 or later) interface.

    For best performance, UMFPACK requires an optimized BLAS library.  It can
    also be compiled without any BLAS at all.  UMFPACK requires AMD Version 2.0.

Authors:

    Timothy A. Davis (davis@cise.ufl.edu), University of Florida.

    Includes a modified version of COLAMD V2.0, by Stefan I. Larimore and
    Timothy A. Davis, University of Florida.  The COLAMD algorithm was developed
    in collaboration with John Gilbert, Xerox Palo Alto Research Center, and
    Esmond Ng, Lawrence Berkeley National Laboratory.

    Includes AMD, by Timothy A. Davis, Patrick R. Amestoy, and Iain S. Duff.

    UMFPACK Version 2.2.1 (MA38 in the Harwell Subroutine Library) is
    co-authored with Iain S. Duff, Rutherford Appleton Laboratory.

Acknowledgements:

    This work was supported by the National Science Foundation, under
    grants DMS-9504974, DMS-9803599, and CCR-0203270.

    Portions of this work were done while on sabbatical at Stanford University
    and Lawrence Berkeley National Laboratory (with funding from the SciDAC
    program).  I would like to thank Gene Golub, Esmond Ng, and Horst Simon
    for making this sabbatical possible.

    I would also like to thank the many researchers who provided sparse
    matrices from a wide range of domains and used earlier versions of UMFPACK/
    MA38 in their applications, and thus assisted in the practical development
    of the algorithm (see http://www.cise.ufl.edu/research/sparse, future
    contributions of matrices are always welcome).

    The MathWorks, Inc., provided a pre-release of MATLAB V6 which allowed me
    to release the first umfpack mexFunction (v3.0) about 6 months earlier than
    I had originally planned.  They also supported the extension of UMFPACK to
    complex, singular, and rectangular matrices (UMFPACK v4.0).

    Penny Anderson (The MathWorks, Inc.), Anshul Gupta (IBM), and Friedrich
    Grund (WAIS) assisted in porting UMFPACK to different platforms.  Penny
    Anderson also incorporated UMFPACK v4.0 into MATLAB, for lu, backslash (\),
    and forward slash (/).

    David Bateman (Motorola) wrote the initial version of the packed complex
    input option, and umfpack_get_determinant.

--------------------------------------------------------------------------------
Files and directories in the UMFPACK distribution:
--------------------------------------------------------------------------------

    ----------------------------------------------------------------------------
    Subdirectories of the UMFPACK directory:
    ----------------------------------------------------------------------------

    Doc		documentation
    Source	primary source code
    Include	include files for use in your code that calls UMFPACK
    Demo	demo programs.  also serves as test of the UMFPACK installation.
    MATLAB	UMFPACK mexFunction for MATLAB, and supporting m-files
    Lib		where the compiled C-callable UMFPACK library is placed.

    ----------------------------------------------------------------------------
    Files in the UMFPACK directory:
    ----------------------------------------------------------------------------

    Makefile	top-level Makefile for GNU make or original make.
		Windows users would require Cygwin to use "make"

    README.txt	this file

    ----------------------------------------------------------------------------
    Doc directory: documentation
    ----------------------------------------------------------------------------

    ChangeLog			change log
    License			the UMFPACK License
    Makefile			for creating the documentation
    QuickStart.tex		Quick Start guide (source)
    QuickStart.pdf		Quick Start guide (PDF)
    UserGuide.bib		User Guide (references)
    UserGuide.sed1		sed script for processing UserGuide.stex
    UserGuide.sed2		sed script for processing UserGuide.stex
    UserGuide.stex		User Guide (LaTeX)
    UserGuide.pdf		User Guide (PDF)

    ----------------------------------------------------------------------------
    Source directory:
    ----------------------------------------------------------------------------

    cholmod_blas.h		an exact copy of CHOLMOD/Include/cholmod_blas.h

    umfpack_col_to_triplet.c	convert col form to triplet
    umfpack_defaults.c		set Control defaults
    umfpack_free_numeric.c	free Numeric object
    umfpack_free_symbolic.c	free Symbolic object
    umfpack_get_determinant.c	compute determinant from Numeric object
    umfpack_get_lunz.c		get nz's in L and U
    umfpack_get_numeric.c	get Numeric object
    umfpack_get_symbolic.c	get Symbolic object
    umfpack_load_numeric.c	load Numeric object from file
    umfpack_load_symbolic.c	load Symbolic object from file
    umfpack_numeric.c		numeric factorization
    umfpack_qsymbolic.c		symbolic factorization, user Q
    umfpack_report_control.c	print Control settings
    umfpack_report_info.c	print Info statistics
    umfpack_report_matrix.c	print col or row-form sparse matrix
    umfpack_report_numeric.c	print Numeric object
    umfpack_report_perm.c	print permutation
    umfpack_report_status.c	print return status
    umfpack_report_symbolic.c	print Symbolic object
    umfpack_report_triplet.c	print triplet matrix
    umfpack_report_vector.c	print dense vector
    umfpack_save_numeric.c	save Numeric object to file
    umfpack_save_symbolic.c	save Symbolic object to file
    umfpack_scale.c		scale a vector
    umfpack_solve.c		solve a linear system
    umfpack_symbolic.c		symbolic factorization
    umfpack_tictoc.c		timer
    umfpack_timer.c		timer
    umfpack_transpose.c		transpose a matrix
    umfpack_triplet_to_col.c	convert triplet to col form

    umf_config.h		configuration file (BLAS, memory, timer)
    umf_internal.h		definitions internal to UMFPACK
    umf_version.h		version definitions (int/UF_long, real/complex)

    umf_2by2.[ch]
    umf_analyze.[ch]		symbolic factorization of A'*A
    umf_apply_order.[ch]	apply column etree postorder
    umf_assemble.[ch]		assemble elements into current front
    umf_blas3_update.[ch]	rank-k update.  Uses level-3 BLAS
    umf_build_tuples.[ch]	construct tuples for elements
    umf_colamd.[ch]		COLAMD pre-ordering, modified for UMFPACK
    umf_create_element.[ch]	create a new element
    umf_dump.[ch]		debugging routines, not normally active
    umf_extend_front.[ch]	extend the current frontal matrix
    umf_free.[ch]		free memory
    umf_fsize.[ch]		determine largest front in each subtree
    umf_garbage_collection.[ch]	compact Numeric->Memory
    umf_get_memory.[ch]		make Numeric->Memory bigger
    umf_grow_front.[ch]		make current frontal matrix bigger
    umf_init_front.[ch]		initialize a new frontal matrix
    umf_is_permutation.[ch]	checks the validity of a permutation vector
    umf_kernel.[ch]		the main numeric factorization kernel
    umf_kernel_init.[ch]	initializations for umf_kernel
    umf_kernel_wrapup.[ch]	wrapup for umf_kernel
    umf_local_search.[ch]	local row and column pivot search
    umf_lsolve.[ch]		solve Lx=b
    umf_ltsolve.[ch]		solve L'x=b and L.'x=b
    umf_malloc.[ch]		malloc some memory
    umf_mem_alloc_element.[ch]		allocate element in Numeric->Memory
    umf_mem_alloc_head_block.[ch]	alloc. block at head of Numeric->Memory
    umf_mem_alloc_tail_block.[ch]	alloc. block at tail of Numeric->Memory
    umf_mem_free_tail_block.[ch]	free block at tail of Numeric->Memory
    umf_mem_init_memoryspace.[ch]	initialize Numeric->Memory
    umf_realloc.[ch]		realloc memory
    umf_report_perm.[ch]	print a permutation vector
    umf_report_vector.[ch]	print a double vector
    umf_row_search.[ch]		look for a pivot row
    umf_scale.[ch]		scale the pivot column
    umf_scale_column.[ch]	move pivot row & column into place, log P and Q
    umf_set_stats.[ch]		set statistics (final or estimates)
    umf_singletons.[ch]		find all zero-cost pivots
    umf_solve.[ch]		solve a linear system
    umf_start_front.[ch]	start a new frontal matrix for one frontal chain
    umf_store_lu.[ch]		store LU factors of current front
    umf_symbolic_usage.[ch]	determine memory usage for Symbolic object
    umf_transpose.[ch]		transpose a matrix in row or col form
    umf_triplet.[ch]		convert triplet to column form
    umf_tuple_lengths.[ch]	determine the tuple list lengths
    umf_usolve.[ch]		solve Ux=b
    umf_utsolve.[ch]		solve U'x=b and U.'x=b
    umf_valid_numeric.[ch]	checks the validity of a Numeric object
    umf_valid_symbolic.[ch]	check the validity of a Symbolic object

    ----------------------------------------------------------------------------
    Include directory:
    ----------------------------------------------------------------------------

    umfpack.h			include file for user programs.  Includes all of
				the following files.  This serves are source-
				code level documenation.  These files are also
				used to construct the User Guide.

    umfpack_col_to_triplet.h
    umfpack_defaults.h
    umfpack_free_numeric.h
    umfpack_free_symbolic.h
    umfpack_get_determinant.h
    umfpack_get_lunz.h
    umfpack_get_numeric.h
    umfpack_get_symbolic.h
    umfpack_load_numeric.h
    umfpack_load_symbolic.h
    umfpack_numeric.h
    umfpack_qsymbolic.h
    umfpack_report_control.h
    umfpack_report_info.h
    umfpack_report_matrix.h
    umfpack_report_numeric.h
    umfpack_report_perm.h
    umfpack_report_status.h
    umfpack_report_symbolic.h
    umfpack_report_triplet.h
    umfpack_report_vector.h
    umfpack_save_numeric.h
    umfpack_save_symbolic.h
    umfpack_scale.h
    umfpack_solve.h
    umfpack_symbolic.h
    umfpack_tictoc.h
    umfpack_timer.h
    umfpack_transpose.h
    umfpack_triplet_to_col.h

    umfpack_wsolve.h		note that there is no umfpack_wsolve.c.  The
				umfpack_*_wsolve routines are created from the
				umfpack_solve.c file.

    ----------------------------------------------------------------------------
    Demo directory:
    ----------------------------------------------------------------------------

    Makefile			for GNU make or original make

    umfpack_simple.c		a simple demo
    umpack_xx_demo.c		template to create the demo codes below

    umfpack_di_demo.sed		for creating umfpack_di_demo.c
    umfpack_dl_demo.sed		for creating umfpack_dl_demo.c
    umfpack_zi_demo.sed		for creating umfpack_zi_demo.c
    umfpack_zl_demo.sed		for creating umfpack_zl_demo.c

    umfpack_di_demo.c		a full demo (real/int version)
    umfpack_dl_demo.c		a full demo (real/UF_long version)
    umfpack_zi_demo.c		a full demo (complex/int version)
    umfpack_zl_demo.c		a full demo (complex/UF_long version)

    umfpack_di_demo.out		umfpack_di_demo output
    umfpack_dl_demo.out		umfpack_dl_demo output
    umfpack_zi_demo.out		umfpack_zi_demo output
    umfpack_zl_demo.out		umfpack_zl_demo output

    umf4.c			a demo (real/int) for Harwell/Boeing matrices
    umf4.out			output of "make hb"
    HB				directory of sample Harwell/Boeing matrices
    readhb.f			reads HB matrices, keeps zero entries
    readhb_nozeros.f		reads HB matrices, removes zero entries
    readhb_size.f		reads HB matrix dimension, nnz
    tmp				empty directory for umf4.c demo

    umf4_f77wrapper.c		a simple FORTRAN interface for UMFPACK.
				compile with "make fortran"
    umf4hb.f			a demo of the FORTRAN interface
    umf4hb.out			output of "make fortran"

    umf4_f77zwrapper.c		a simple FORTRAN interface for the complex
				UMFPACK routines.  compile with "make fortran"
    umf4zhb.f			a demo of the FORTRAN interface (complex)
    umf4zhb.out			output of umf4zhb with HB/qc324.cua

    umf4hb64.f			64-bit version of umf4hb.f

    simple_compile		a single command that compiles the double/int
				version of UMFPACK (useful prototype for
				Microsoft Visual Studio project)

    ----------------------------------------------------------------------------
    MATLAB directory:
    ----------------------------------------------------------------------------

    Contents.m			for "help umfpack" listing of toolbox contents
    GNUmakefile			a nice Makefile, for GNU make
    Makefile			an ugly Unix Makefile (for older make's)

    lu_normest.m		1-norm estimate of A-L*U (by Hager & Davis).
    luflop.m			for "help luflop"
    luflopmex.c			luflop mexFunction, for computing LU flop count
    umfpack.m			for "help umfpack"
    umfpack_btf.m		solve Ax=b using umfpack and dmperm
    umfpack_demo.m		a full umfpack demo
    umfpack_details.m		the details of how to use umfpack
    umfpack_make.m		compile the umfpack mexFunction within MATLAB
    umfpack_report.m		report statistics
    umfpack_simple.m		a simple umfpack demo
    umfpack_solve.m		x=A\b or b/A for arbitrary b
    umfpack_test.m		extensive test, requires UF sparse matrices
    umfpackmex.c		the umfpack mexFunction
    west0067.mat		sparse matrix for umfpack_demo.m

    umfpack_demo.m.out		output of umfpack_demo.m
    umfpack_simple.m.out	output of umfpack_simple

    lcc_lib/lapacksyms.def	LAPACK definitions for lcc compiler (Windows)
    lcc_lib/libmwlapack.lib	LAPACK definitions for lcc compiler (Windows)

    ----------------------------------------------------------------------------
    Lib directory:  libumfpack.a library placed here
    ----------------------------------------------------------------------------

    GNUmakefile			a nice Makefile, for GNU make
    Makefile			an ugly Unix Makefile (for older make's)
    libumfpack.def		UMPFACK definitions for Windows
