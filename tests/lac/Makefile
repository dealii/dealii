############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
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

tests_x = vector-vector \
	full_matrix tridiagonal_matrix \
	sparsity_pattern sparse_matrices \
	block_vector block_vector_iterator block_matrices \
	matrix_lib matrix_out \
	solver eigen \
	sparse_ilu lapack


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