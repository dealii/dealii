# $Id$

CXX=c++
INCPATH   = -I../include -I../../lac/include -I../../base/include
CXXFLAGS.g= -DDEBUG -g -Wall $(INCPATH)
CXXFLAGS  = -O3 $(INCPATH)

LIBPATH   = -L../lib -L../../lac/lib
LIBS.g    = $(LIBPATH) -lbasic.g -lgrid.g -lfe.g -lnumerics.g -llac.g
LIBS      = $(LIBPATH) -lbasic -lgrid -lfe -lnumerics 

ifeq ($(shell uname),SunOS)
CXXFLAGS  := -V2.7.2.3 $(CXXFLAGS)
CXXFLAGS.g:= -V2.7.2.3 $(CXXFLAGS.g)
endif


cc-files = grid/grid_test.cc dof/dof_test.cc poisson/poisson.cc
o-files  = $(cc-files:.cc=.o)
h-files  = $(wildcard ../include/*.h)

examples: $(cc-files:.cc=)

%.o : %.cc
	@echo ================= Compiling $<
	@$(CXX) $(CXXFLAGS.g) -o $@ -g -c $<

grid/grid_test: grid/grid_test.o
	@echo ================= Linking $@
	@$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS)

dof/dof_test: dof/dof_test.o
	@echo ================= Linking $@
	@$(CXX) $(CXXFLAGS.g) -o $@ $< $(LIBS.g) -lg++

poisson/poisson: poisson/poisson.o
	@echo ================= Linking $@
	@$(CXX) $(CXXFLAGS.g) -o $@ $< ../../mia/control.o $(LIBS.g) -lg++


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
	@echo ================= Making $@
	@$(CXX) $(INCPATH) -M $(cc-files) > Makefile.dep

include Makefile.dep
