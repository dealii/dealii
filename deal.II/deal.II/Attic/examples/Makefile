# $Id$

CXX=c++
INCLUDE   = -I../include -I../../lac/include -I../../base/include
CXXFLAGS.g= -DDEBUG -g -Wall -W -pedantic -Wconversion \
            -Winline -Woverloaded-virtual -fno-rtti -fno-exceptions
CXXFLAGS  =-O3 -Wuninitialized -finline-functions -ffast-math \
           -funroll-loops -felide-constructors -fnonnull-objects \
           -fno-rtti -fno-exceptions

LIBPATH   = -L../lib -L../../lac/lib
LIBFILES.g= ../lib/libbasic.g.a ../lib/libgrid.g.a \
            ../lib/libfe.g.a ../lib/libnumerics.g.a \
            ../../lac/lib/liblac.g.a
LIBFILES  = $(LIBFILES.g:.g.a=.a)
LIBS.g    = $(LIBPATH) -lgrid.g -lfe.g -lnumerics.g -llac.g -lbasic.g
LIBS      = $(LIBS.g:.g=)


ifeq ($(shell uname),Linux)
CXX       = /home/wolf/bin/gcc-2.8.1/bin/c++
endif

ifeq ($(shell uname),SunOS)
CXX       = /usr/local/gcc/gcc-2.8.1/bin/c++
endif

cc-files = grid/grid_test.cc dof/dof_test.cc poisson/poisson.cc convergence/convergence.cc
o-files  = $(cc-files:.cc=.o)
h-files  = $(wildcard ../include/*.h)

examples: $(cc-files:.cc=)

%.o : %.cc
	@echo ================= Compiling $<
	@$(CXX) $(CXXFLAGS.g) $(INCLUDE) -o $@ -g -c $<

grid/grid_test: grid/grid_test.o $(LIBFILES)
	@echo ================= Linking $@
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS.g)

dof/dof_test: dof/dof_test.o $(LIBFILES.g)
	@echo ================= Linking $@
	@$(CXX) $(CXXFLAGS.g) -o $@ $< $(LIBS.g) -lg++

poisson/poisson:
	@cd poisson ; make

convergence/convergence: convergence/convergence.o $(LIBFILES.g)
	@echo ================= Linking $@
	@$(CXX) $(CXXFLAGS.g) -o $@ $< $(LIBS) ../../mia/control.o


run: run_grid_test run_dof_test run_poisson_test run_convergence_test

run_grid_test:
	cd grid ; grid_test 4 ; mv grid.1 grid.4
	cd grid ; grid_test 3 ; mv grid.1 grid.3
	cd grid ; grid_test 2 ; mv grid.1 grid.2
	cd grid ; grid_test 1
	cd grid ; gnuplot make_ps

run_dof_test:
	cd dof ; dof_test dof_test.prm
	cd dof ; gnuplot make_ps

run_poisson_test:
	cd poisson ; make run

run_convergence_test:
	cd convergence ; convergence
	cd convergence ; gnuplot make_ps

clean:
	cd grid ; rm -f grid.[1234] *.eps *.o *~ grid_test
	cd dof  ; rm -f grid.* sparsity.* *.o *~ dof_test
	cd poisson ; make clean
	cd convergence ; rm -f *.o convergence *~ *.eps gnuplot.*


.PHONY: run run_grid_test run_dof_test run_poisson_test run_convergence_test clean


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
Makefile.dep: $(cc-files) $(h-files)
	@echo ================= Making $@
	@$(CXX) $(INCLUDE) -M $(cc-files) > Makefile.dep

include Makefile.dep
