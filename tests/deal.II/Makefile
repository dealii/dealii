# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998

# use debug mode by default for all testcases
debug-mode = on

###############################################################################
# Internals

D = ../..

include ../../common/Make.global_options




# get lists of files we need
cc-files    = $(shell echo *.cc)
lib-h-files = $(shell echo ../../deal_II/include/*/*.h) \
              $(shell echo ../../lac/include/*/*.h) \
              $(shell echo ../../base/include/*/*.h)
executables = $(cc-files:.cc=.testcase)
output-files= $(cc-files:.cc=.output)

# list of libraries needed to link with
libs.g   = $D/lac/lib/liblac.g.a $D/base/lib/libbase.g.a
libs     = $D/lac/lib/liblac.a   $D/base/lib/libbase.a


# check whether we use debug mode or not
ifeq ($(debug-mode),on)
libraries = $(libs.g)
flags     = $(CXXFLAGS.g) -L$D/deal.II/lib
lib-tag   = "deal_II_libraries.g="
endif

ifeq ($(debug-mode),off)
libraries = $(libs)
flags     = $(CXXFLAGS.o) -L$D/deal.II/lib
lib-tag   = "deal_II_libraries="
endif

%.go : %.cc Makefile
	@echo =====debug========= $<
	@$(CXX) $(flags) -c $< -o $@
%.o : %.cc Makefile
	@echo =====optimized===== $<
	@$(CXX) $(flags) -c $< -o $@



# rule to make executables
%.testcase: %.cc $(libraries)
	@echo ============================ Compiling Testcase:   $<
	@$(CXX) $(flags) $< -o $@ $(shell grep $(lib-tag) $< | perl -pe 's/.*=//;') $(libraries)




#################################### Here comes the interesting part

# make a rule for each executable
targets = $(addsuffix .target, $(executables))

# define the action of the targets: execute it. do so by re-getting
# the right filename by dropping the suffix added for the make-rule
%.target : %
	@echo ============================ Executing Testcase:   $<
	@./$< > $(subst testcase,output,$<)
	@diff $(subst testcase,output,$<) $(subst testcase,expect,$<)

# this is the main target
# make dependence so that it first wants to compile all files then
# execute them
check: $(executables) $(targets)



clean:
	-rm -f *.o *.go *~ Make.depend Makefile.dep $(executables) $(output-files) tmp/*



.PHONY: check clean


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
#
#use perl to generate rules for the .go files as well
#as to make rules not for tria.o and the like, but
#rather for libnumerics.a(tria.o)
Makefile.dep: $(cc-files) $(h-files) $(lib-h-files) Makefile
	@echo ============================ Remaking Makefile
	@perl ../../common/scripts/Make_dep.pl ./Obj $(INCLUDE) $(cc-files) | \
		perl -pi -e 's?^./Obj.a\((.*)\.o\):?$$1.testcase:?g;'   \
		> Makefile.dep


include Makefile.dep

