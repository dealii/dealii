# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998

#this file uses many small rules to allow
#   make -j
#on multiprocessor machines. Feel free to use this feature,
#it gives you quicker compile time and annoys those FEATFLOW
#people ;-)
#
#I think I have removed any non-supported cross-references, so you
#should be able to compile with an arbitrary number of processes.

deal_II_dimension=2


cc-files = $(wildcard */*.cc)
h-files  = $(wildcard ../include/*/*.h)

# strip subdirectories from cc file names, preprend a "../lib/[g]o"
# and change the suffix. We place the object files with or without
# debug info in different directories, since we want to use the
# -frepo switch to gcc; this, however, creates a .rpo file which
# would be the same for both debug and non-debug versions, if we
# don't care about different directories
tmp1     = $(notdir $(cc-files))
o-files  = $(addprefix ../lib/$(deal_II_dimension)d/o/,  $(tmp1:.cc=_$(deal_II_dimension)d.o) )
go-files = $(addprefix ../lib/$(deal_II_dimension)d/go/, $(tmp1:.cc=_$(deal_II_dimension)d.go))


CXX       = c++
INCLUDE   = -I../include -I../../lac/include \
            -I../../base/include -I../../mia/include
CXXFLAGS.g= -DDEBUG -g -Wall -W -pedantic -Wconversion \
            -Winline -Woverloaded-virtual \
            $(INCLUDE) -Ddeal_II_dimension=$(deal_II_dimension)
CXXFLAGS  =-O3 -Wuninitialized -finline-functions -ffast-math \
           -felide-constructors -fnonnull-objects \
           $(INCLUDE) \
           -Ddeal_II_dimension=$(deal_II_dimension)

%.go :
	@echo ============================ Compiling with debugging information:   $<
	@$(CXX) $(CXXFLAGS.g) -c $< -o $@
%.o :
	@echo ============================ Compiling with optimization:   $<
	@$(CXX) $(CXXFLAGS) -c $< -o $@






all: ../lib/libdeal_II.g.a ../lib/libdeal_II.a

../lib/libdeal_II.g.a: $(go-files)
	@echo ============================ Updating library:   $@
	@ar ruv $@ $(go-files)

../lib/libdeal_II.a: $(o-files)
	@echo ============================ Updating library:   $@
	@ar ruv $@ $(o-files)


clean:
	@cd basic   ; rm -f *.o *.go *~
	@cd grid    ; rm -f *.o *.go *~
	@cd fe      ; rm -f *.o *.go *~
	@cd numerics; rm -f *.o *.go *~



.PHONY: all clean


#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
#
#use perl to generate rules for the .go files as well
#as to make rules not for tria.o and the like, but
#rather for libgrid.a(tria.o)
#
#Later addition: the perl script was made to compile
#right into a library. Since this is not what we want
#anymore, but since 1) I am too lazy to change the
#script and 2) others may use the script as it was,
#we don't change the script but rather pipe its
#result through another perl script to get what we
#want.
#
#Next addition: replace .g?o by _2d.g?o or whatever
#the dimension is, because we want to have unique
#object file names in the libraries in order to write
#programs which use the 2d *and* the 3d library;
#if the files within these libraries had the same names,
#we may get into trouble.
../lib/Makefile.dep: $(cc-files) $(h-files) Makefile ../Make.global_options
	@echo ============================ Remaking Makefile
	@perl ../Make_dep.pl ../lib/libgrid $(INCLUDE) $(cc-files) \
		| perl -p -e 's!^.*\(.*/(.*)\):!../lib/$(deal_II_dimension)d/$$1:!g;' \
		| perl -pe 's!(/[^/]+\.(o|go)):!/$$2$$1:!g;' \
		| perl -pe 's!(\.g?o)!_$(deal_II_dimension)d$$1!g;' \
		> ../lib/Makefile.dep


include ../lib/Makefile.dep
