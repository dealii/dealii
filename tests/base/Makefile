# Generated automatically from Makefile.in by configure.
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

libraries = $(lib-lac.g)  \
            $(lib-base.g)

default: run-tests 

############################################################

logtest.exe        : logtest.go            $(lib-base.g)
polynomial_test.exe: polynomial_test.go    $(lib-base.g)
polynomial1d.exe   : polynomial1d.go       $(lib-base.g)
quadrature_test.exe: quadrature_test.go    $(libraries)
reference.exe      : reference.go abort.go $(libraries)
tensor.exe         : tensor.go             $(libraries)
timer.exe          : timer.go              $(libraries)
auto_derivative_function.exe        : auto_derivative_function.go            $(libraries)


tests = logtest reference quadrature_test tensor timer polynomial1d polynomial_test auto_derivative_function

############################################################


include ../Makefile.rules
include Makefile.depend
