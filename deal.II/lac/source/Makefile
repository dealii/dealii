# $Id$

CXX       = c++
INCLUDE   = -I../include -I../../base/include
CXXFLAGS.g= -DDEBUG -g -Wall -W -pedantic -Wconversion \
            -Winline -Woverloaded-virtual -fno-rtti \
            $(INCLUDE)
CXXFLAGS  =-O3 -Wuninitialized -finline-functions -ffast-math \
           -funroll-loops -felide-constructors -fnonnull-objects \
           -fno-rtti $(INCLUDE)

ifeq ($(shell uname),Linux)
CXX       = /home/wolf/bin/gcc/bin/c++
endif

ifeq ($(shell uname),SunOS)
CXX       = g++
endif


cc-files = $(wildcard *.cc)
o-files  = $(cc-files:.cc=.o)
go-files = $(cc-files:.cc=.go)
h-files  = $(wildcard ../include/*/*.h)


%.go : %.cc Makefile
	@echo ============================ Compiling with debugging information:   $<
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o : %.cc Makefile
	@echo ============================ Compiling with optimization:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@


lib: lib.g.a lib.a


lib.a: ../lib/liblac.a($(o-files))

lib.g.a: ../lib/liblac.g.a($(go-files))

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
../lib/Makefile.dep: $(cc-files) $(h-files)
	@echo ============================ Remaking Makefile
	@perl ../Make_dep.pl ../lib/liblac $(INCLUDE) $(cc-files) \
		> ../lib/Makefile.dep


include ../lib/Makefile.dep
