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
INCLUDE2.8  =  -I/usr/local/source/libstdc++-2.8.0/libstdc++ \
               -I/usr/local/source/libstdc++-2.8.0/libstdc++/stl \
               -I/usr/local/source/libstdc++-2.8.0/libio
CXXFLAGS  := $(CXXFLAGS) $(INCLUDE2.8)
CXXFLAGS.g:= $(CXXFLAGS.g) $(INCLUDE2.8)
endif

LIBPATH   = -L../../lib -L../../../lac/lib
LIBFILES.g= ../../lib/libbasic.g.a ../../lib/libgrid.g.a \
            ../../lib/libfe.g.a ../../lib/libnumerics.g.a \
            ../../../lac/lib/liblac.g.a
LIBFILES  = $(LIBFILES.g:.g.a=.a)
LIBS.g    = $(LIBPATH) -lbasic.g -lgrid.g -lfe.g -lnumerics.g -llac.g
LIBS      = $(LIBS.g:.g=)


cc-files = $(wildcard *.cc)
o-files  = $(cc-files:.cc=.o)
go-files = $(cc-files:.cc=.go)
h-files  = $(wildcard ../include/*/*.h)




%.go : %.cc poisson.h
	@echo ============================ Compiling with debugging information:   $<
	@echo $(CXX) ... -c $< -o $@
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o : %.cc poisson.h
	@echo ============================ Compiling with optimization:   $<
	@echo $(CXX) ... -c $< -o $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@


poisson: $(go-files) $(LIBFILES.g)
	@echo ============================ Linking poisson
	@$(CXX) -g $(CXXFLAGS.g) $(go-files) -o $@ ../../../mia/control.o $(LIBS.g) -lg++


run: poisson
	poisson poisson.prm
	gnuplot make_ps


clean:
	rm -f *.o *.go *~ poisson gnuplot.* *.eps



.PHONY: run clean
