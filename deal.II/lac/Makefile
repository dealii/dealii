# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999


# deal.II top level directory
D=..


# create lists of file names. note that we need several output
# files for data_out_base.cc, which is why we use the several
# intermediate o-files-base* rules
cc-files     = $(shell echo source/*.cc)
o-files-base1= $(addprefix lib/o/, $(patsubst source/%,%,$(cc-files:.cc=.o)))
o-files-base2= $(filter-out %data_out_base.o, $(o-files-base1))
o-files-base3= $(sort $(o-files-base2))
o-files      = $(o-files-base3)
go-files     = $(patsubst lib/o/%,lib/go/%,$(o-files-base3:.o=.go))
h-files      = $(shell echo include/lac/*.h)

forward-declarations = $D/lac/include/lac/forward_declarations.h

# get options which hold for all files of the project
include $D/common/Make.global_options


# production rules
lib/go/%.go :
	@echo =====lac==============debug========= $(<F)
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS.g) -c $< -o $@
lib/o/%.o :
	@echo =====lac==============optimized===== $(<F)
	@$(CXX) $(CXXFLAGS) $(CXXFLAGS.o) -c $< -o $@


# rules
lib: lib/liblac.g.a lib/liblac.a
libg: lib/liblac.g.a

lib/liblac.a: $(forward-declarations) $(o-files)
	@echo =====lac==============optimized===== Updating library:   $@
	@ar ruv $@ $(o-files)

lib/liblac.g.a: $(forward-declarations) $(go-files)
	@echo =====lac==============debug========= Updating library:   $@
	@ar ruv $@ $(go-files)


clean:
	-rm -f *~ include/*~ include/*/*~ lib/lib* lib/o/*.o lib/go/*.go \
               include/base/forward_declarations.h*



#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
lib/Makefile.dep: $(cc-files) $(forward-declarations) $(h-files) Makefile
	@echo ============================ Remaking Makefile
	@$(PERL) $D/common/scripts/make_dependencies.pl $(INCLUDE) $(cc-files) \
		> lib/Makefile.dep

# include all the dependencies
include lib/Makefile.dep

# include a rule which tells how to remake the forward declarations file
include $D/common/scripts/forward_declarations


.PHONY: clean
.PHONY: lib lib.a lib.g.a
