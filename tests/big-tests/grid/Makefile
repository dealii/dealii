# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998

# Template for makefiles for the examples subdirectory. In principle,
# everything should be done automatically if you set the target file
# here correctly:
target   = grid_test

# All dependencies between files should be updated by the included
# file Makefile.dep if necessary. Object files are compiled into
# the archives ./Obj.a and ./Obj.g.a. By default, the debug version
# is used to link. It you don't like that, change the following
# variable to "off"
debug-mode = on

# If you want your program to be linked with extra object or library
# files, specify them here:
user-libs =

# To run the program, use "make run"; to give parameters to the program,
# give the parameters to the following variable:
run-parameters  = 2 1

# To execute additional action apart from running the program, fill
# in this list:
additional-run-action = grid_test 2 2 ; grid_test 2 3 ; grid_test 2 4 ; grid_test 3 1 ; grid_test 3 2 ; grid_test 3 3 ; grid_test 3 4

# To specify which files are to be deleted by "make clean" (apart from
# the usual ones: object files, executables, backups, etc), fill in the
# following list
delete-files = results/[123]d*




###############################################################################
# Internals

#deal include base path
D = ../../..

include $D/common/Make.global_options


# get lists of files we need
cc-files    = $(wildcard *.cc)
o-files     = $(cc-files:.cc=.o)
go-files    = $(cc-files:.cc=.go)
h-files     = $(wildcard *.h)
lib-h-files = $(wildcard ../../include/*/*.h)

# list of libraries needed to link with
libs     = ./Obj.a   ../../lib/libdeal_II_2d.g.a ../../lib/libdeal_II_3d.a   ../../../lac/lib/liblac.a   ../../../base/lib/libbase.a
libs.g   = ./Obj.g.a ../../lib/libdeal_II_2d.g.a ../../lib/libdeal_II_3d.g.a ../../../lac/lib/liblac.g.a ../../../base/lib/libbase.g.a


# check whether we use debug mode or not
ifeq ($(debug-mode),on)
libraries = $(libs.g)
flags     = $(CXXFLAGS.g)
endif

ifeq ($(debug-mode),off)
libraries = $(libs)
flags     = $(CXXFLAGS)
endif



# make rule for the target
$(target) : $(libraries)
	@echo ============================ Linking $@
	@$(CXX) $(flags) -o $@ $(libraries) $(user-libs)

# rule how to run the program
run: $(target)
	$(target) $(run-parameters)
	$(additional-run-action)


# rule to make object files
%.go : %.cc
	@echo ============================ Compiling with debugging information:   $<
	@echo $(CXX) ... -c $< -o $@
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o : %.cc
	@echo ============================ Compiling with optimization:   $<
	@echo $(CXX) ... -c $< -o $@
	@$(CXX) $(CXXFLAGS) -c $< -o $@


# rules which files the libraries depend upon
Obj.a: ./Obj.a($(o-files))
Obj.g.a: ./Obj.g.a($(go-files))


clean:
	rm -f *.o *.go *~ Makefile.dep Obj.a Obj.g.a $(target) $(delete-files)



.PHONY: clean


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
#
#use perl to generate rules for the .go files as well
#as to make rules not for tria.o and the like, but
#rather for libnumerics.a(tria.o)
Makefile.dep: $(cc-files) $(h-files) $(lib-h-files)
	@echo ============================ Remaking Makefile
	@perl $D/common/scripts/Make_dep.pl ./Obj $(INCLUDE) $(cc-files) \
		> Makefile.dep


include Makefile.dep

