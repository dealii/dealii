############################################################
# Makefile,v 1.15 2002/06/13 12:51:13 hartmann Exp
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
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

tests_x=fe_data_test traits fe_tools fe_tools_* mapping \
        mapping_c1 shapes derivatives function numbering mapping_*_eulerian \
        transfer internals derivatives_face \
	dgq_1 \
	q_* \
	interpolate_q1 \
	rt_* \
	abf_* \
	fe_sys_* \
	rtdiff \
	rtn_* \
	interpolate_rt interpolate_rtn \
	dgp_monomial_1 dgp_monomial_2 \
	system_index system_1 interpolate_system \
        non_primitive_* \
        nedelec* \
	up_and_down \
	check_derivatives \
	injection_* \
	copy_* \
	cell_similarity_*

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))

############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
