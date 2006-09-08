############################################################
# $Id: Makefile 13843 2006-09-06 21:29:34Z kanschat $
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
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
libraries_o = $(lib-deal2-1d.o) \
            $(lib-deal2-2d.o) \
            $(lib-deal2-3d.o) \
            $(lib-lac.o)  \
            $(lib-base.o)

default: run-tests 

############################################################

tests_x = umfpack*


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