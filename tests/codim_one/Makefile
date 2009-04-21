############################################################
# $Id: Makefile 16047 2008-05-07 16:51:50Z heltai $
# Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries = $(lib-deal2-1d.g) \
            $(lib-deal2-2d.g) \
            $(lib-deal2-3d.g) \
            $(lib-lac.g)  \
            $(lib-base.g)

default: run-tests 

############################################################

tests_x = grid* \
	mapping* \
	interpolation* \
	surface* \
	data_out* \
	projection*\
	integrate*\
	gradients*\
	tensor_matrix_conversion*\
	bem*



# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))


############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
