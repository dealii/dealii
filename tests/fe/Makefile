# Generated automatically from Makefile.in by configure.
############################################################
# $Id$
# Copyright (C) 2000, 2001 by the deal.II authors
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

derivatives.exe         : derivatives.go         $(libraries)
fe_data_test.exe        : fe_data_test.go        $(libraries)
mapping.exe             : mapping.go             $(libraries)
mapping_c1.exe          : mapping_c1.go          $(libraries)
shapes.exe              : shapes.go              $(libraries)
numbering.exe           : numbering.go           $(libraries)
mapping_q1_eulerian.exe : mapping_q1_eulerian.go $(libraries)

tests = derivatives fe_data_test mapping mapping_c1 shapes \
        numbering mapping_q1_eulerian

############################################################


include ../Makefile.rules
include Makefile.depend
