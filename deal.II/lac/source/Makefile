# $Id$

D=../..

cc-files = $(shell echo *.cc)
o-files  = $(addprefix ../lib/o/, $(cc-files:.cc=.o))
go-files = $(addprefix ../lib/go/, $(cc-files:.cc=.go))
h-files  = $(shell echo ../include/*/*.h)

forward-declarations = $D/lac/include/lac/forward-declarations.h

include $D/common/Make.global_options


../lib/go/%.go :
	@echo ============================ Compiling with debugging information:   $<
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS.g) -c $< -o $@
../lib/o/%.o :
	@echo ============================ Compiling with optimization:   $<
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS.o) -c $< -o $@


lib: ../lib/liblac.g.a ../lib/liblac.a
libg: ../lib/liblac.g.a

../lib/liblac.a: $(forward-declarations) $(o-files)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(o-files)

../lib/liblac.g.a: $(forward-declarations) $(go-files)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(go-files)

clean:
	rm -f ../lib/o/*.o ../lib/o/*.go *~ Makefile.dep ../lib/lib*





.PHONY: lib lib.a lib.g.a clean


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
#
#see the Makefile of the deal.II directory for a thorough
#description of what happens here
../lib/Makefile.dep: $(cc-files) $(forward-declarations) $(h-files) Makefile
ifneq (1,${recursive-make-fwd-decl})
	@echo ============================ Remaking Makefile
	@perl $D/common/scripts/Make_dep.pl ../lib/liblac $(INCLUDE) $(cc-files)               \
		| perl -p -e 's!^.*\((.*\.(g?o))\):!../lib/$$2/$$1:!g;' \
		> ../lib/Makefile.dep
else
	@echo ============================ Not remaking Makefile since in nested mode
endif


include ../lib/Makefile.dep

include $D/common/Makescripts/forward_declarations

