# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999



# deal.II top level directory
D=..


# create lists of file names. note that we need several output
# files for data_out_base.cc, which is why we use the several
# intermediate o-files-base* rules
cc-files     = $(shell echo source/*.cc)
o-files-base1= $(addprefix lib/o/, $(patsubst source/%,%,$(cc-files:.cc=.o)))
o-files-base2= $(filter-out %data_out_base.o, $(o-files-base1)) \
               $(patsubst %.o,%_1d.o,$(filter %data_out_base.o, $(o-files-base1))) \
               $(patsubst %.o,%_2d.o,$(filter %data_out_base.o, $(o-files-base1))) \
               $(patsubst %.o,%_3d.o,$(filter %data_out_base.o, $(o-files-base1))) \
               $(patsubst %.o,%_4d.o,$(filter %data_out_base.o, $(o-files-base1)))
o-files-base3= $(sort $(o-files-base2))
o-files      = $(o-files-base3)
go-files     = $(patsubst lib/o/%,lib/go/%,$(o-files-base3:.o=.go))
h-files      = $(shell echo include/base/*.h)

forward-declarations = $D/base/include/base/forward-declarations.h

# get options which hold for all files of the project
include $D/common/Make.global_options


# urgh, this is a weird hack, since this function uses the $@ variable
# which is only set in rules generating output files.
#
# the two variables here do the following: maybe-define-dimension1
# filters out the 1, 2, or 3 if it is contained within some file
# name sufficed by _1d, _2d, or _3d. It contains three parts,
# but it is obvious that at most one branch will succeed.
# This variable will then contain several spaces and maybe one
# digit, which are removed by the 'strip' command
#
# maybe-define-dimension will substitute anything by
# -Ddata_out_dimension=anything, if there is an 'anything'. If not,
# i.e. if the filename contained no _*d, then nothing is
# substituted and the result of the variable is an empty string
maybe-define-dimension1 = $(strip $(findstring 1,$(findstring _1d.,$@)) \
                                  $(findstring 2,$(findstring _2d.,$@)) \
                                  $(findstring 3,$(findstring _3d.,$@)) \
                                  $(findstring 4,$(findstring _4d.,$@)))
maybe-define-dimension  = $(patsubst %,-Ddata_out_dimension=%,$(maybe-define-dimension1))

# production rules
lib/go/%.go :
	@echo =================== Compiling with debugging information:   $(<F) \
	      "   " $(patsubst -Ddata_out_dimension=%,dimension=%,$(maybe-define-dimension))
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS.g) $(maybe-define-dimension) -c $< -o $@
lib/o/%.o :
	@echo =================== Compiling with optimization:   $(<F) \
	      "   " $(patsubst -Ddata_out_dimension=%,dimension=%,$(maybe-define-dimension))
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS.o) $(maybe-define-dimension) -c $< -o $@


# rules
lib: lib/libbase.g.a lib/libbase.a
libg: lib/libbase.g.a

lib/libbase.a: $(forward-declarations) $(o-files)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(o-files)

lib/libbase.g.a: $(forward-declarations) $(go-files)
	@echo ======================================== Updating library:   $@
	@ar ruv $@ $(go-files)


clean:
	-rm -f *~ include/*~ include/*/*~ lib/lib* lib/o/*.o lib/go/*.go include/base/forward-declarations.h*



#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
#
#note that here we translate file names with data_out_base
#to three file names with the dimension for which we compile
#appended
lib/Makefile.dep: $(cc-files) $(forward-declarations) $(h-files) Makefile
	@echo ============================ Remaking Makefile
	@perl $D/common/scripts/Make_dep_2.pl $(INCLUDE) $(cc-files) \
		| perl -p -e 's!(lib/g?o/data_out_base)(.g?o)!$$1_1d$$2 $$1_2d$$2 $$1_3d$$2 $$1_4d$$2!g;' \
		> lib/Makefile.dep


# include all the dependencies
include lib/Makefile.dep

# include a rule which tells how to remake the forward declarations file
include $D/common/Makescripts/forward_declarations


.PHONY: clean
.PHONY: lib lib.a lib.g.a
