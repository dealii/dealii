# $Id$

CXX       = c++
INCLUDE   = -I../../include -I../../../lac/include \
            -I../../../base/include -I../../../mia/include
CXXFLAGS.g= -DDEBUG -g -Wall -W -pedantic -Wconversion \
            -Winline -Woverloaded-virtual -fno-rtti -fno-exceptions \
            $(INCLUDE)
CXXFLAGS  =-O3 -Wuninitialized -finline-functions -ffast-math \
           -funroll-loops -felide-constructors -fnonnull-objects \
           -fno-rtti -fno-exceptions $(INCLUDE)

ifeq ($(shell uname),Linux)
CXX       = /home/wolf/bin/gcc-2.8.1/bin/c++
endif

ifeq ($(shell uname),SunOS)
CXX       = /usr/local/gcc/gcc-2.8.1/bin/c++
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

lib.a: ../../lib/libnumerics.a($(o-files))

lib.g.a: ../../lib/libnumerics.g.a($(go-files))

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
#rather for libnumerics.a(tria.o)
Makefile.dep: $(cc-files) $(h-files)
	@echo ============================ Remaking Makefile
	@perl ../../Make_dep.pl ../../lib/libnumerics $(INCLUDE) $(cc-files) \
		> Makefile.dep


include Makefile.dep
