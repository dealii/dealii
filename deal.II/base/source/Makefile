# $Id$

D=../..

include $D/deal.II/Make.global_options

cc-files = $(wildcard *.cc)
o-files  = $(cc-files:.cc=.o)
go-files = $(cc-files:.cc=.go)
h-files  = $(wildcard ../include/*/*.h)


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

# rule how to make the file containing all the forward declarations
../include/base/forward-declarations.h: $(filter-out %forward-declarations.h,$(h-files))
	@echo ============================ Generating $@
	@perl ../../deal.II/Make_forward_declarations $(filter-out %forward-declarations.h,$(h-files)) > $@


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
../lib/Makefile.dep: $(cc-files) $(h-files)
	@echo ============================ Remaking Makefile
	@perl ../Make_dep.pl ../lib/libbase $(INCLUDE) $(cc-files) \
		> ../lib/Makefile.dep


include ../lib/Makefile.dep
