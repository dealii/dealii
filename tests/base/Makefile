# Generated automatically from Makefile.in by configure.
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

libraries += $(lib-lac.g)  \
             $(lib-base.g)

default: run-tests 

############################################################

threads.cc : threads.pl
	@echo =====regenerating== $@
	@$(PERL) $^ > $@

data_out_base.exe  : data_out_base.g.$(OBJEXT)                      $(lib-base.g)
logtest.exe        : logtest.g.$(OBJEXT)                            $(lib-base.g)
polynomial_test.exe: polynomial_test.g.$(OBJEXT)                    $(lib-base.g)
polynomial1d.exe   : polynomial1d.g.$(OBJEXT)                       $(lib-base.g)
quadrature_test.exe: quadrature_test.g.$(OBJEXT)                    $(lib-base.g)
quadrature_selector.exe: quadrature_selector.g.$(OBJEXT)            $(lib-base.g)
reference.exe      : reference.g.$(OBJEXT)                          $(lib-base.g)
table.exe          : table.g.$(OBJEXT)                              $(lib-base.g)
tensor.exe         : tensor.g.$(OBJEXT)                             $(lib-base.g)
timer.exe          : timer.g.$(OBJEXT)                              $(lib-base.g)
threads.exe        : threads.g.$(OBJEXT)                            $(lib-base.g)
anisotropic_1.exe  : anisotropic_1.g.$(OBJEXT)                      $(lib-base.g)
anisotropic_2.exe  : anisotropic_2.g.$(OBJEXT)                      $(lib-base.g)
hierarchical.exe   : hierarchical.g.$(OBJEXT)                       $(lib-base.g)
bdm.exe            : bdm.g.$(OBJEXT)                                $(lib-base.g)
auto_derivative_function.exe : auto_derivative_function.g.$(OBJEXT) $(libraries) $(lib-base.g)

tests =  logtest reference quadrature_test quadrature_selector table tensor \
          timer threads polynomial1d polynomial_test \
	  auto_derivative_function anisotropic_1 anisotropic_2 \
          hierarchical data_out_base

############################################################


include ../Makefile.rules
include Makefile.depend
