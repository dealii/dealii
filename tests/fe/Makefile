# Generated automatically from Makefile.in by configure.
############################################################
# Makefile,v 1.15 2002/06/13 12:51:13 hartmann Exp
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
fe_tools_test.exe       : fe_tools_test.go       $(libraries)
fe_traits.exe           : fe_traits.go           $(libraries)
internals.exe           : internals.go           $(libraries)
mapping.exe             : mapping.go             $(libraries)
mapping_c1.exe          : mapping_c1.go          $(libraries)
mapping_q1_eulerian.exe : mapping_q1_eulerian.go $(libraries)
nedelec.exe             : nedelec.go             $(libraries)
nedelec_2.exe           : nedelec_2.go           $(libraries)
nedelec_3.exe           : nedelec_3.go           $(libraries)
non_primitive_1.exe     : non_primitive_1.go     $(libraries)
non_primitive_2.exe     : non_primitive_2.go     $(libraries)
numbering.exe           : numbering.go           $(libraries)
shapes.exe              : shapes.go              $(libraries)
transfer.exe		: transfer.go		 $(libraries)
up_and_down.exe         : up_and_down.go	 $(libraries)

tests = fe_data_test fe_traits fe_tools_test mapping \
        mapping_c1 shapes derivatives numbering mapping_q1_eulerian \
        transfer internals \
        non_primitive_1 non_primitive_2 \
        nedelec nedelec_2 nedelec_3 up_and_down

############################################################


include ../Makefile.rules
include Makefile.depend
