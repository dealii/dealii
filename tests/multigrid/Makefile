############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006, 2007 by the deal.II authors
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

tests_x = cycles dof_* count_* boundary_* renumbering_* \
	transfer* \
	smoother_block \
	mg_output \
	mg_renumbered*

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))

############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
