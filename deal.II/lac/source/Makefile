# $Id$

D=../..

cc-files = $(wildcard *.cc)
o-files  = $(cc-files:.cc=.o)
go-files = $(cc-files:.cc=.go)
h-files  = $(wildcard ../include/*/*.h)

forward-declarations = $D/lac/include/lac/forward-declarations.h

include $D/deal.II/Make.global_options


%.go : %.cc Makefile
	@echo ============================ Compiling with debugging information:   $<
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o : %.cc Makefile
	@echo ============================ Compiling with optimization:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@


lib: lib.g.a lib.a


lib.a: $(forward-declarations) ../lib/liblac.a($(o-files))

lib.g.a: $(forward-declarations) ../lib/liblac.g.a($(go-files))

clean:
	rm -f *.o *.go *~ Makefile.dep ../lib/lib*





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
	@perl ../Make_dep.pl ../lib/liblac $(INCLUDE) $(cc-files) \
		> ../lib/Makefile.dep
else
	@echo ============================ Not remaking Makefile since in nested mode
endif


include ../lib/Makefile.dep

include $D/deal.II/Make.forward-declarations
