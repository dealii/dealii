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

D=../..

cc-files = $(wildcard */*.cc)
h-files  = $(wildcard ../include/*/*.h)

forward-declarations = $D/deal.II/include/basic/forward-declarations.h

include $D/deal.II/Make.global_options




# strip subdirectories from cc file names, preprend a "../lib/[g]o"
# and change the suffix. We place the object files with or without
# debug info in different directories, since we want to use the
# -frepo switch to gcc; this, however, creates a .rpo file which
# would be the same for both debug and non-debug versions, if we
# don't care about different directories
tmp1        = $(notdir $(cc-files))
o-files-1d  = $(addprefix ../lib/1d/o/,  $(tmp1:.cc=_1d.o) )
go-files-1d = $(addprefix ../lib/1d/go/, $(tmp1:.cc=_1d.go))
o-files-2d  = $(addprefix ../lib/2d/o/,  $(tmp1:.cc=_2d.o) )
go-files-2d = $(addprefix ../lib/2d/go/, $(tmp1:.cc=_2d.go))
o-files-3d  = $(addprefix ../lib/3d/o/,  $(tmp1:.cc=_3d.o) )
go-files-3d = $(addprefix ../lib/3d/go/, $(tmp1:.cc=_3d.go))





default: 2d
all: 1d 2d 3d

1d: ../lib/libdeal_II_1d.g.a ../lib/libdeal_II_1d.a

2d: ../lib/libdeal_II_2d.g.a ../lib/libdeal_II_2d.a

3d: ../lib/libdeal_II_3d.g.a ../lib/libdeal_II_3d.a

# Debug version only

debug: 2dg 1dg 3dg

1dg: ../lib/libdeal_II_1d.g.a
2dg: ../lib/libdeal_II_2d.g.a
3dg: ../lib/libdeal_II_3d.g.a


# rules how to generate object files from source files. note that
# there are some files which needs exceptional rules; these are
# listed immediately below
../lib/1d/go/%.go :
	@echo ==============1d======debug============= $<
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=1 -c $< -o $@
../lib/1d/o/%.o :
	@echo ==============1d======optimized========= $<
	@$(CXX) $(CXXFLAGS) -Ddeal_II_dimension=1 -c $< -o $@

../lib/2d/go/%.go :
	@echo ==============2d======debug============= $<
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=2 -c $< -o $@
../lib/2d/o/%.o :
	@echo ==============2d======optimized========= $<
	@$(CXX) $(CXXFLAGS) -Ddeal_II_dimension=2 -c $< -o $@

../lib/3d/go/%.go :
	@echo ==============3d======debug============= $<
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=3 -c $< -o $@
../lib/3d/o/%.o :
	@echo ==============3d======optimized========= $<
	@$(CXX) $(CXXFLAGS) -Ddeal_II_dimension=3 -c $< -o $@

# special rules for exceptional files
../lib/3d/o/q1_mapping.jacobians_3d.o:
	@echo "==============3d======opt==(special)====" $<
	@$(CXX) $(filter-out -O2 -Wuninitialized,$(CXXFLAGS)) -Ddeal_II_dimension=3 -c $< -o $@



# rules how to make the libraries themselves
../lib/libdeal_II_1d.g.a: $(go-files-1d)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(go-files-1d)
../lib/libdeal_II_1d.a: $(o-files-1d)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(o-files-1d)

../lib/libdeal_II_2d.g.a: $(go-files-2d)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(go-files-2d)
../lib/libdeal_II_2d.a: $(o-files-2d)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(o-files-2d)

../lib/libdeal_II_3d.g.a: $(go-files-3d)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(go-files-3d)
../lib/libdeal_II_3d.a: $(o-files-3d)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(o-files-3d)


clean:
	@cd basic   ; rm -f *.o *.go *~
	@cd grid    ; rm -f *.o *.go *~
	@cd fe      ; rm -f *.o *.go *~
	@cd numerics; rm -f *.o *.go *~



.PHONY: default 1d 2d 3d all clean


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
#
#
#additional clue: if the forward declarations file does
#not exist, then the script generating the dependencies
#will not find it and will not print it into the dependencies.
#but then also 'make' will not know that it needs to be
#generated before compilation of files start, which leads
#to a cpp error. we therefore have to make sure the forward
#declaration file exists before we even start to build the
#dependencies; since it will not be in the $(h-files) in the
#critical case, we need to add it by hand. next problem:
#the generation of the forward declarations requires a
#nested call to 'make', which would then rebuild the
#Makefile.dep, which this rule won't see, so we end up
#generating it twice. therefore: do nothing if in a nested
#make
../lib/Makefile.dep: $(cc-files) $(h-files) $(forward-declarations) Makefile ../Make.global_options
ifneq (1,${recursive-make-fwd-decl})
	@echo ============================ Remaking Makefile
	@echo "# This Makefile was automatically generated by ../source/Makefile" \
		> ../lib/Makefile.dep
	@echo "# See there for more information." >> ../lib/Makefile.dep
	@perl ../Make_dep.pl ../lib/libgrid $(INCLUDE) $(cc-files) \
		| perl -p -e 's!^.*\(.*/(.*)\):!../lib/DIM_PLACEHOLDER/$$1:!g;' \
		| perl -pe 's!(/[^/]+\.(o|go)):!/$$2$$1:!g;' \
		| perl -pe 's!(\.g?o)!_DIM_PLACEHOLDER$$1!g;' \
	        | perl -pe 's!^\.\./lib/DIM_PLACEHOLDER/(.*):!../lib/1d/$$1 ../lib/2d/$$1 ../lib/3d/$$1:!g;' \
	        | perl -pe 's!(\.\./lib/1d/[^ ]+)_DIM_PLACEHOLDER!$$1_1d!g;' \
	        | perl -pe 's!(\.\./lib/2d/[^ ]+)_DIM_PLACEHOLDER!$$1_2d!g;' \
	        | perl -pe 's!(\.\./lib/3d/[^ ]+)_DIM_PLACEHOLDER!$$1_3d!g;' \
		>> ../lib/Makefile.dep
	@echo ============================ Done
else
	@echo ============================ Not remaking Makefile since in nested mode
endif


include ../lib/Makefile.dep

include $D/deal.II/Make.forward-declarations

