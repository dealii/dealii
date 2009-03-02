############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009 by the deal.II authors
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

tests_x =  geometry_info_* \
	   logtest \
	   log_* \
	   reference \
	   quadrature_* qprojector \
	   path_search \
	   slice_vector \
	   table \
	   tensor \
	   tensor_* \
	   symmetric_tensor_* \
	   full_tensor_* \
           timer \
	   polynomial* \
	   auto_derivative_function \
	   anisotropic_* \
           hierarchical \
	   data_out_base \
	   data_out_base_dx \
	   data_out_base_eps \
	   data_out_base_gmv \
	   data_out_base_gnuplot \
	   data_out_base_povray \
	   data_out_base_tecplot \
	   data_out_base_vtk \
	   function_* \
           functions_* \
	   utilities_* \
	   is_finite* \
	   number_traits_*

ifeq ($(enable-multithreading),no)
  tests_x += threads
endif

ifeq (tecio, $(findstring tecio, $LIBS))
  tests_x += data_out_base_tecplot_bin
endif

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))


############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
