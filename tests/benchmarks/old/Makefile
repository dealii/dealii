############################################################
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

include ../Makefile.paths
include $D/common/Make.global_options
debug-mode = on

libraries = $(lib-deal2-2d.o) $(lib-deal2-3d.o) $(lib-lac.o) $(lib-base.o)
libraries.g = $(lib-deal2-2d.g) $(lib-deal2-3d.g) $(lib-lac.g) $(lib-base.g)

default: run-tests 

############################################################

tests_x =  dof_handler_*

# from above list of regular expressions, generate the real set of
# tests
expand = $(shell echo $(addsuffix .cc,$(1)) \
           | $(PERL) -pi -e 's/\.cc//g;')
tests = $(call expand,$(tests_x))


############################################################
# First how to create executables, including all necessary
# flags:
############################################################

flags     = $(CXXFLAGS.g)

ifeq ($(findstring gcc,$(GXX_VERSION)),gcc)
flags += -Wno-missing-noreturn
endif

%/obj.g.$(OBJEXT) : %.cc
	@echo =====debug========= $<
	@echo Compiling > $*/status
	@$(CXX) $(flags) -c $< -o $@

%/obj.$(OBJEXT) : %.cc
	@echo =====optimized===== $<
	@echo Compiling > $*/status
	@$(CXX) $(CXXFLAGS.o) -c $< -o $@

%/exe : %/obj.$(OBJEXT) $(libraries)
	@echo =====linking======= $@
	@echo Linking > $*/status
	@$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

%/g.exe : %/obj.g.$(OBJEXT) $(libraries.g)
	@echo =====linking======= $@
	@echo Linking > $*/status
	@$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

%/output : %/exe
	@echo =====running======= $<
	@$<

############################################################
# After all these general rules, here is the target to be
# executed by make: for each entry in the list $(tests)
# perform a check.
############################################################
run-tests: $(tests:%=%/output)

build: $(tests:%=%/exe)

refs: $(tests:%=%/ref)


Makefile.depend: $(shell echo *.cc)
	@echo =====Dependencies== $@
	@$D/common/scripts/make_dependencies $(INCLUDE) $^ > $@ \
		|| (rm -f $@ ; false)
	@$(PERL) -pi -e 's#\.((g\.)?$(OBJEXT):)#/obj.\1#g;' $@


-include Makefile.depend
