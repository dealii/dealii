# $Id$
# Copyright W. Bangerth, 2003


# deal.II top level directory
D=../..

include $D/common/Make.global_options



make_dependencies : make_dependencies.cc
	@echo ============================ Compiling $@
	@$(CXX) $(CXXFLAGS.o) $< -o $@
