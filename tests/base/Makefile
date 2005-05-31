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

libraries += $(lib-lac.g)  \
             $(lib-base.g)

default: run-tests 

############################################################

threads.cc : threads.pl
	@echo =====regenerating== $@
	@$(PERL) $^ > $@

tests_x =  logtest \
	   reference \
	   quadrature_* qprojector \
	   path_search \
	   slice_vector \
	   table \
	   tensor \
	   symmetric_tensor_* \
	   full_tensor_* \
           timer \
	   threads \
	   polynomial* \
	   auto_derivative_function \
	   anisotropic_* \
           hierarchical \
	   data_out_base \
	   function_parser

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))


############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
