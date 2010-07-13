# $Id$
# Copyright W. Bangerth, 2003


# deal.II top level directory
D=../..

include $D/common/Make.global_options

default: make_dependencies$(EXEEXT) expand_instantiations$(EXEEXT) \
         report_features$(EXEEXT)


# rules to generate the program. let is be dependent on the global flags for
# the following reason: if you go to a different system, but keep with the
# same directory, you get into trouble because make_dependencies is not
# deleted by "make clean"; following this, we generate empty Makefile.dep's
# because we can't run the executable. So better be careful and regenerate
# make_dependencies whenever we think it may be beneficial
make_dependencies$(EXEEXT) : make_dependencies.cc $D/common/Make.global_options
	@echo ============================ Compiling $@
	@$(CXX) $(CXXFLAGS.o) $(LDFLAGS) $< -o $@

expand_instantiations$(EXEEXT) : expand_instantiations.cc $D/common/Make.global_options
	@echo ============================ Compiling $@
	@$(CXX) $(CXXFLAGS.o) $(LDFLAGS) $< -o $@

report_features$(EXEEXT) : report_features.cc $D/common/Make.global_options
	@echo ============================ Compiling $@
	@$(CXX) $(CXXFLAGS.o) $(LDFLAGS) $< -o $@

clean:
	-rm -f make_dependencies
	-rm -f make_dependencies.o
	-rm -f expand_instantiations
	-rm -f expand_instantiations.o

.PHONY: clean
