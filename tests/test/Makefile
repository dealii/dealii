############################################################
# $Id: Makefile 13841 2006-09-06 21:04:19Z kanschat $
# Copyright (C) 2006, 2010 by the deal.II authors
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

# all .cc-files are tests:
tests = $(basename $(wildcard *.cc))

############################################################
# phony header file to force compilation whenever
# make is invoked
############################################################

phony.h:
	touch phony.h

############################################################
# Make sure even numbers are are always recompiled
############################################################

no_*_?[02468]/obj.g.o:     phony.h

############################################################
# Another odd target that sets al OK files to the OK status
# Try 'make report' after doing this
############################################################

all-ok:
	@for test in $(tests); do echo 'diff ok' > $$test/OK ; done


include ../Makefile.rules
include Makefile.depend
include Makefile.tests

.PHONY: default phony.h
