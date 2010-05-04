############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
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

tests_x = hp_dof_handler \
	  crash_* \
	  q_collection_* \
	  mapping_collection_* \
	  fe_collection_* \
	  fe_nothing_* \
	  continuous_* \
	  n_dofs \
	  random \
	  step-* \
	  n_active_fe_indices \
	  integrate_difference \
	  hp_*_dof_identities_* \
	  hp_hanging_nodes_* \
	  interpolate_* \
	  hp_constraints_* \
	  matrices* \
	  vectors* \
	  get_active_fe_indices \
	  create_* \
	  renumber_* \
	  dof_renumbering_* \
	  compare*

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))

############################################################


include ../Makefile.rules
include Makefile.depend
include Makefile.tests
