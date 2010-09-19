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

# all .cc-files are tests:
tests = $(basename $(wildcard *.cc))

############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
