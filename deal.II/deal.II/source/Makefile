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

forward-declarations = ../include/basic/forward-declarations.h

include $D/deal.II/Make.global_options
include $D/deal.II/Make.forward-declarations




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


# rule how to make the file containing all the forward declarations
# this rule is rather complicated, it works as follows:
#
# - if this is the outermost 'make', then create a file $@.old.
#   if the forward declaration file $@ already exists, then copy
#   it over, otherwise leave it empty. this way, we always have
#   a file to run diff on
# - write the forward declarations for the presently available
#   header files to $@.new
# - we would like to compare the output of the script generating
#   the forward declarations with the old file, or the empty
#   string if there is no such file; this is why we created $@.old
# - if the files differ, we would like to copy over the old file.
# - this can either be done using shell programming using 'if [...]'
#   or using conditionals inside 'make'. the first way is not 
#   portable over different shells, however, and I also do not know
#   how to write it. the second way has a problem also: conditionals
#   are evaluated at the time the Makefile is read first, so we
#   can't compare against $@.old, since this file does not exist at
#   that time.
# - we therefore invoke a second 'make', which can compare the two
#   files against each other in a preconditional clause and
#   copy over the file if necessary.
# - delete the two temporary files
$(forward-declarations): $(filter-out $(forward-declarations),$(h-files))
ifeq (0,${MAKELEVEL})
	@echo ============================ Checking $@
	@touch $@.old
	-@cp $@ $@.old
	@perl ../Make_forward_declarations $(filter-out $(forward-declarations),$(h-files)) > $@.new
	@make --silent $@
	@rm $@.old $@.new
else
  ifneq (,$(shell diff $(forward-declarations).old $(forward-declarations).new))
	@echo ============================ Generating $@
	@cp $@.new $@
  else
	@echo "============================          $@ needs no change"
  endif
endif




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
../lib/3d/o/fe_linear_mapping.jacobians_3d.o:
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
../lib/Makefile.dep: $(cc-files) $(filter-out $(forward-declarations),$(h-files)) Makefile ../Make.global_options
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


include ../lib/Makefile.dep
