# Generated automatically from Makefile.in by configure.
############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002 by the deal.II authors
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
fe_traits.exe           : fe_traits.go           $(libraries)
internals.exe           : internals.go           $(libraries)
mapping.exe             : mapping.go             $(libraries)
mapping_c1.exe          : mapping_c1.go          $(libraries)
mapping_q1_eulerian.exe : mapping_q1_eulerian.go $(libraries)
non_primitive_1.exe     : non_primitive_1.go     $(libraries)
numbering.exe           : numbering.go           $(libraries)
shapes.exe              : shapes.go              $(libraries)

tests = derivatives fe_data_test fe_traits internals \
        mapping mapping_c1 mapping_q1_eulerian  \
        non_primitive_1 numbering shapes 

############################################################


include ../Makefile.rules
include Makefile.depend
