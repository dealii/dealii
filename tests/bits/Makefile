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

anna_1.exe              : anna_1.g.$(OBJEXT)              $(libraries)
anna_2.exe              : anna_2.g.$(OBJEXT)              $(libraries)
anna_3.exe              : anna_3.g.$(OBJEXT)              $(libraries)
anna_5.exe              : anna_5.g.$(OBJEXT)              $(libraries)
gerold_1.exe            : gerold_1.g.$(OBJEXT)            $(libraries)
roy_1.exe               : roy_1.g.$(OBJEXT)               $(libraries)

tests = anna_1 anna_2 anna_3 anna_5 \
	gerold_1 \
	roy_1

############################################################


include ../Makefile.rules
include Makefile.depend
