# $Id$

CXX       =c++
CXXFLAGS.g=-DDEBUG -g -Wall -pedantic -Wconversion -Winline -Woverloaded-virtual -I../../include
CXXFLAGS  =-O3 -Wuninitialized -finline-functions -ffast-math -funroll-loops -felide-constructors -fnonnull-objects -I../../include

ifeq ($(shell uname),SunOS)
CXXFLAGS  := -V2.7.2.3 $(CXXFLAGS)
CXXFLAGS.g:= -V2.7.2.3 $(CXXFLAGS.g)
endif


cc-files = $(wildcard *.cc)
o-files  = $(cc-files:.cc=.o)
go-files = $(cc-files:.cc=.go)
h-files  = $(wildcard ../include/*/*.h)


%.go : %.cc
	@echo ============================ Compiling with debugging information:   $<
	@echo $(CXX) ... -c $< -o $@
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o : %.cc
	@echo ============================ Compiling with optimization:   $<
	@echo $(CXX) ... -c $< -o $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@


lib: lib.g.a lib.a

lib.a: ../../lib/libbasic.a($(o-files))

lib.g.a: ../../lib/libbasic.g.a($(go-files))

clean:
	rm -f *.o *.go *~ Makefile.dep



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
Makefile.dep: $(cc-files) $(h-files)
	@echo ============================ Remaking Makefile
	@$(CXX) -I../../include -DDEBUG -M *.cc > Makefile.dep.o
	@cp Makefile.dep.o Makefile.dep.go
	@perl -pi -e 's!^(.*\.o):!../../lib/libbasic.a(\1):!g' Makefile.dep.o
	@perl -pi -e 's!\.o:!.go:!g' Makefile.dep.go
	@perl -pi -e 's!^(.*\.go):!../../lib/libbasic.g.a(\1):!g' Makefile.dep.go
	@cat Makefile.dep.o Makefile.dep.go > Makefile.dep
	@rm Makefile.dep.o Makefile.dep.go

include Makefile.dep
