# Generated automatically from Makefile.in by configure.
############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries = $(lib-lac.g)  \
            $(lib-base.g)
lib-1d    = $(lib-deal2-1d.g)
lib-2d    = $(lib-deal2-2d.g)
lib-3d    = $(lib-deal2-3d.g)

default: run-tests 

############################################################

transfer.exe      : transfer.g.$(OBJEXT) $(lib-2d)           $(libraries)
mg1.exe           : mg1.g.$(OBJEXT)      $(lib-2d)           $(libraries)

tests = transfer

############################################################


include ../Makefile.rules
include Makefile.depend
