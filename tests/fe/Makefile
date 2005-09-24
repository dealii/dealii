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

%.exe: %.g.$(OBJEXT) $(libraries)
	@echo =====linking======= $@
	@$(CXX) $(LDFLAGS) -o $@ ../abort.o $^ $(LIBS)

############################################################

tests = fe_data_test traits fe_tools fe_tools_test mapping \
        mapping_c1 shapes derivatives function numbering mapping_q1_eulerian \
        transfer internals \
	dgq_1 \
	q_1 q_2 q_3 q_4 interpolate_q1 \
	rt_1 rt_2 rt_3 rt_4 rt_5 rt_6 rtdiff rtn_1 interpolate_rtn \
	dgp_monomial_1 dgp_monomial_2 \
	system_1 interpolate_system \
        non_primitive_1 non_primitive_2 \
        nedelec nedelec_2 nedelec_3 up_and_down

############################################################


include ../Makefile.rules
include Makefile.depend
