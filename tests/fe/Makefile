############################################################
# Makefile,v 1.15 2002/06/13 12:51:13 hartmann Exp
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
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
fe_tools.exe            : fe_tools.g.$(OBJEXT)            $(libraries)
fe_tools_test.exe       : fe_tools_test.g.$(OBJEXT)       $(libraries)
function.exe            : function.g.$(OBJEXT)            $(libraries)
traits.exe              : traits.g.$(OBJEXT)              $(libraries)
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
dgq_1.exe               : dgq_1.g.$(OBJEXT)               $(libraries)
q_1.exe                 : q_1.g.$(OBJEXT)                 $(libraries)
q_2.exe                 : q_2.g.$(OBJEXT)                 $(libraries)
q_3.exe                 : q_3.g.$(OBJEXT)                 $(libraries)
q_4.exe                 : q_4.g.$(OBJEXT)                 $(libraries)
rt_1.exe                : rt_1.g.$(OBJEXT)                $(libraries)
rt_2.exe                : rt_2.g.$(OBJEXT)                $(libraries)
rt_3.exe                : rt_3.g.$(OBJEXT)                $(libraries)
rt_4.exe                : rt_4.g.$(OBJEXT)                $(libraries)
rt_5.exe                : rt_5.g.$(OBJEXT)                $(libraries)
rt_6.exe                : rt_6.g.$(OBJEXT)                $(libraries)
rtn_1.exe               : rtn_1.g.$(OBJEXT)               $(libraries)
rtdiff.exe              : rtdiff.g.$(OBJEXT)              $(libraries)
dgp_monomial_1.exe	: dgp_monomial_1.g.$(OBJEXT)      $(libraries)
dgp_monomial_2.exe	: dgp_monomial_2.g.$(OBJEXT)      $(libraries)
system_1.exe            : system_1.g.$(OBJEXT)            $(libraries)
shapes.exe              : shapes.g.$(OBJEXT)              $(libraries)
transfer.exe		: transfer.g.$(OBJEXT)	          $(libraries)
up_and_down.exe         : up_and_down.g.$(OBJEXT)         $(libraries)

tests = fe_data_test traits fe_tools fe_tools_test mapping \
        mapping_c1 shapes derivatives function numbering mapping_q1_eulerian \
        transfer internals \
	dgq_1 \
	q_1 q_2 q_3 q_4 \
	rt_1 rt_2 rt_3 rt_4 rt_5 rt_6 rtdiff rtn_1 \
	dgp_monomial_1 dgp_monomial_2 \
	system_1 \
        non_primitive_1 non_primitive_2 \
        nedelec nedelec_2 nedelec_3 up_and_down

############################################################


include ../Makefile.rules
include Makefile.depend
