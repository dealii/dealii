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

derivatives.exe         : derivatives.g.$(OBJEXT)         $(libraries)
fe_data_test.exe        : fe_data_test.g.$(OBJEXT)        $(libraries)
fe_tools_test.exe       : fe_tools_test.g.$(OBJEXT)       $(libraries)
fe_traits.exe           : fe_traits.g.$(OBJEXT)           $(libraries)
internals.exe           : internals.g.$(OBJEXT)           $(libraries)
mapping.exe             : mapping.g.$(OBJEXT)             $(libraries)
mapping_c1.exe          : mapping_c1.g.$(OBJEXT)          $(libraries)
mapping_q1_eulerian.exe : mapping_q1_eulerian.g.$(OBJEXT) $(libraries)
nedelec.exe             : nedelec.g.$(OBJEXT)             $(libraries)
nedelec_2.exe           : nedelec_2.g.$(OBJEXT)           $(libraries)
nedelec_3.exe           : nedelec_3.g.$(OBJEXT)           $(libraries)
non_primitive_1.exe     : non_primitive_1.g.$(OBJEXT)     $(libraries)
non_primitive_2.exe     : non_primitive_2.g.$(OBJEXT)     $(libraries)
numbering.exe           : numbering.g.$(OBJEXT)           $(libraries)
shapes.exe              : shapes.g.$(OBJEXT)              $(libraries)
transfer.exe		: transfer.g.$(OBJEXT)	          $(libraries)
up_and_down.exe         : up_and_down.g.$(OBJEXT)         $(libraries)

tests = fe_data_test fe_traits fe_tools_test mapping \
        mapping_c1 shapes derivatives numbering mapping_q1_eulerian \
        transfer internals \
        non_primitive_1 non_primitive_2 \
        nedelec nedelec_2 nedelec_3 up_and_down

############################################################


include ../Makefile.rules
include Makefile.depend
