# $Id$

D=../..

cc-files = $(wildcard *.cc)
o-files  = $(cc-files:.cc=.o)
go-files = $(cc-files:.cc=.go)
h-files  = $(wildcard ../include/*/*.h)

forward-declarations = ../include/base/forward-declarations.h

include $D/deal.II/Make.global_options
include $D/deal.II/Make.forward-declarations



%.go : %.cc Makefile
	@echo ======== DEBUG =================== $<
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o : %.cc Makefile
	@echo ======== OPTIMIZE ================ $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@


lib: lib.g.a lib.a


lib.a: ../lib/libbase.a($(o-files))

lib.g.a: ../lib/libbase.g.a($(go-files))

clean:
	rm -f *.o *.go *~ Makefile.dep ../lib/lib*



.PHONY: lib lib.a lib.g.a clean


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
#
#use perl to generate rules for the .go files as well
#as to make rules not for tria.o and the like, but
#rather for libbasic.a(tria.o)
../lib/Makefile.dep: $(cc-files) $(filter-out $(forward-declarations),$(h-files)) Makefile
	@echo ============================ Remaking Makefile
	@perl ../Make_dep.pl ../lib/libbase $(INCLUDE) $(cc-files) \
		> ../lib/Makefile.dep


include ../lib/Makefile.dep
