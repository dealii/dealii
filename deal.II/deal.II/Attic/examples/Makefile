# $Id$

CXX=c++
CXXFLAGS.g= -DDEBUG -g -Wall -I../include
CXXFLAGS  = -O3 -I../include

LIBS.g    = -L../lib -lbasic.g -lgrid.g -lfe.g -lnumerics.g -ldeal.g
LIBS      = -L../lib -lbasic -lgrid -lfe -lnumerics  -ldeal

ifeq ($(shell uname),SunOS)
CXXFLAGS  := -V2.7.2.3 $(CXXFLAGS)
CXXFLAGS.g:= -V2.7.2.3 $(CXXFLAGS.g)
endif


cc-files = grid/grid_test.cc
o-files  = $(cc-files:.cc=.o)
h-files  = $(wildcard ../include/*.h)

examples: grid/grid_test dof/dof_test


grid/grid_test: grid/grid_test.o ../lib/libgrid.a
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

dof/dof_test: dof/dof_test.o ../lib/libgrid.g.a ../lib/libdeal.g.a
	$(CXX) $(CXXFLAGS.g) -o $@ $< $(LIBS.g) -lg++


run: run_grid_test run_dof_test

run_grid_test:
	cd grid ; grid_test 4 ; mv grid.1 grid.4
	cd grid ; grid_test 3 ; mv grid.1 grid.3
	cd grid ; grid_test 2 ; mv grid.1 grid.2
	cd grid ; grid_test 1
	cd grid ; gnuplot make_ps

run_dof_test:
	cd dof ; dof_test dof_test.prm
	cd dof ; gnuplot make_ps

clean:
	cd grid ; rm -f grid.[1234] *.eps *.o *~ grid_test
	cd dof  ; rm -f grid.* sparsity.* *.o *~ dof_test


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
Makefile.dep: $(cc-files) $(h-files)
	$(CXX) -I../include -M $(cc-files) > Makefile.dep

include Makefile.dep
