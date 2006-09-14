############################################################
# $Id: Makefile 13841 2006-09-06 21:04:19Z kanschat $
# Copyright (C) 2006 by the deal.II authors
############################################################

############################################################
# Test if test suite reports the expected results
#
# This directory contains tests failing at different stages.
# When running all tests, this should be reflected in the results.
#
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

default: run-tests 

############################################################

tests_x = no_compile_* no_link_* no_run_* no_diff_* ok_*

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))

############################################################
# phony header file to force compilation whenever
# make is invoked
############################################################

phony.h:
	touch phony.h


include ../Makefile.rules
include Makefile.depend
include Makefile.tests

.PHONY: default phony.h
