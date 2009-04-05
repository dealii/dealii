############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries = $(lib-lac.g)  \
            $(lib-base.g)
libraries_o = $(lib-lac.o)  \
            $(lib-base.o)

default: run-tests 

############################################################

tests_x = vector-vector vector_memory \
	full_matrix \
	householder \
	tridiagonal_matrix \
	tridiagonal_matrix_sym \
	sparsity_pattern \
	sparsity_pattern_* \
	chunk_sparsity_pattern_* \
	compressed_sparsity_pattern_* \
	compressed_simple_sparsity_pattern_* \
	compressed_set_sparsity_pattern_* \
	constraints_zero* \
	constraint_graph_zero \
	sparse_matrices matrices \
	block_sparsity_pattern_01 \
	block_indices \
	block_vector \
	block_vector_* \
	block_matrices \
	is_block_matrix \
	matrix_lib \
	matrix_out \
	pointer_matrix \
	pointer_matrix_vector \
	solver \
	eigen \
	filtered_matrix \
	sparse_ilu* \
	lapack_fill \
	block_minres \
	identity_matrix* \
	trace \
	print_formatted_* \
	bicgstab_early \
	complex* \
	vector_view*


# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))

############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests

.PHONY: default
