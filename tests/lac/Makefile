############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
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

# all .cc-files are tests:
tests = $(basename $(wildcard *.cc))

############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests

.PHONY: default
