############################################################
# $Id$
# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007 by the deal.II authors
############################################################

############################################################
# Include general settings for including DEAL libraries
############################################################

D=../../deal.II
include $D/common/Make.global_options
debug-mode = on

libraries += $(lib-deal2-2d.g) $(lib-deal2-3d.g) $(lib-lac) $(lib-base)

default: run-tests 

############################################################

tests_x =  dof_handler_size

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

######################################################################
# Don't put $(libraries) into this line since it is already in $^
######################################################################
%/exe : %/obj.$(OBJEXT)
	@echo =====linking======= $@
	@echo Linking > $*/status
	@$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)


############################################################
# After all these general rules, here is the target to be
# executed by make: for each entry in the list $(tests)
# perform a check.
############################################################
run-tests: $(tests:%=%/output)

build: $(tests:%=%/exe)

refs: $(tests:%=%/ref)


Makefile.tests: Makefile $(shell echo *.cc)
	@echo =====Targets======= $@
	@for i in $(tests) ; do \
		echo "$$i/exe : $$i/obj.\$$(OBJEXT) \$$(libraries)"; \
	 done \
		> $@




include Makefile.depend
include Makefile.tests
