# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998, 1999


# deal.II top level directory
D=..


# top level targets
default: 2dg 2d

1dg: lib/libdeal_II_1d.g.a 
1d : lib/libdeal_II_1d.a
2dg: lib/libdeal_II_2d.g.a 
2d : lib/libdeal_II_2d.a
3dg: lib/libdeal_II_3d.g.a 
3d : lib/libdeal_II_3d.a


# create lists of file names. note that we need several output
# files for data_out_base.cc, which is why we use the several
# intermediate o-files-base* rules
cc-files     = $(shell echo source/*/*.cc)
h-files      = $(shell echo include/*/*.h)

forward-declarations = $D/deal.II/include/grid/forward_declarations.h


# strip subdirectories from cc file names, preprend a "lib/[g]o"
# and change the suffix. We place the object files with or without
# debug info in different directories, since we want to use the
# -frepo switch to gcc; this, however, creates a .rpo file which
# would be the same for both debug and non-debug versions, if we
# don't care about different directories
tmp1        = $(notdir $(cc-files))
o-files-1d  = $(addprefix lib/1d/o/,  $(tmp1:.cc=_1d.o) )
go-files-1d = $(addprefix lib/1d/go/, $(tmp1:.cc=_1d.go))
o-files-2d  = $(addprefix lib/2d/o/,  $(tmp1:.cc=_2d.o) )
go-files-2d = $(addprefix lib/2d/go/, $(tmp1:.cc=_2d.go))
o-files-3d  = $(addprefix lib/3d/o/,  $(tmp1:.cc=_3d.o) )
go-files-3d = $(addprefix lib/3d/go/, $(tmp1:.cc=_3d.go))

o-files  = $(o-files-1d)  $(o-files-2d)  $(o-files-3d)
go-files = $(go-files-1d) $(go-files-2d) $(go-files-3d)

# get options which hold for all files of the project
include $D/common/Make.global_options


# rules how to generate object files from source files. note that
# there are some files which needs exceptional rules; these are
# listed immediately below
lib/1d/go/%.go :
	@echo =====deal.II====1d====debug========= $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=1 -c $< -o $@
lib/1d/o/%.o :
	@echo =====deal.II====1d====optimized===== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=1 -c $< -o $@

lib/2d/go/%.go :
	@echo =====deal.II====2d====debug========= $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=2 -c $< -o $@
lib/2d/o/%.o :
	@echo =====deal.II====2d====optimized===== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=2 -c $< -o $@

lib/3d/go/%.go :
	@echo =====deal.II====3d====debug========= $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.g) -Ddeal_II_dimension=3 -c $< -o $@
lib/3d/o/%.o :
	@echo =====deal.II====3d====optimized===== $(<:source/%=%)
	@$(CXX) $(CXXFLAGS.o) -Ddeal_II_dimension=3 -c $< -o $@

# special rules for exceptional files
lib/3d/o/q1_mapping.jacobians_3d.o:
	@echo =====deal.II====3d====opt==special== $(<:source/%=%)
	@$(CXX) $(filter-out -O2 -Wuninitialized,$(CXXFLAGS.o)) -Ddeal_II_dimension=3 -c $< -o $@


# rules how to make the libraries themselves
lib/libdeal_II_1d.g.a: $(go-files-1d)
	@echo =====deal.II====1d====debug========= Updating library:   $(@F)
	@ar ruv $@ $(go-files-1d)
lib/libdeal_II_1d.a: $(o-files-1d)
	@echo =====deal.II====1d====optimized===== Updating library:   $(@F)
	@ar ruv $@ $(o-files-1d)

lib/libdeal_II_2d.g.a: $(go-files-2d)
	@echo =====deal.II====2d====debug========= Updating library:   $(@F)
	@ar ruv $@ $(go-files-2d)
lib/libdeal_II_2d.a: $(o-files-2d)
	@echo =====deal.II====2d====optimized===== Updating library:   $(@F)
	@ar ruv $@ $(o-files-2d)

lib/libdeal_II_3d.g.a: $(go-files-3d)
	@echo =====deal.II====3d====debug========= Updating library:   $(@F)
	@ar ruv $@ $(go-files-3d)
lib/libdeal_II_3d.a: $(o-files-3d)
	@echo =====deal.II====3d====optimized===== Updating library:   $(@F)
	@ar ruv $@ $(o-files-3d)



#Rule to generate the dependency file. This file is
#automagically remade whenever needed, i.e. whenever
#one of the cc-/h-files changed. Make detects whether
#to remake this file upon inclusion at the bottom
#of this file.
lib/Makefile.dep: $(cc-files) $(forward-declarations) $(h-files) Makefile
	@echo ============================ Remaking Makefile
	@$(PERL) $D/common/scripts/make_dependencies.pl $(INCLUDE) $(cc-files) \
		| $(PERL) -p -e 's!^lib/(g?o/.*):!lib/DIM_PLACEHOLDER/$$1:!g;' \
		| $(PERL) -pe 's!(\.g?o):!_DIM_PLACEHOLDER$$1:!g;' \
	        | $(PERL) -pe 's!^lib/DIM_PLACEHOLDER/(.*):!lib/1d/$$1 lib/2d/$$1 lib/3d/$$1:!g;'\
	        | $(PERL) -pe 's!(lib/1d/[^ ]+)_DIM_PLACEHOLDER!$$1_1d!g;' \
	        | $(PERL) -pe 's!(lib/2d/[^ ]+)_DIM_PLACEHOLDER!$$1_2d!g;' \
	        | $(PERL) -pe 's!(lib/3d/[^ ]+)_DIM_PLACEHOLDER!$$1_3d!g;' \
		> lib/Makefile.dep

# include all the dependencies
include lib/Makefile.dep

# include a rule which tells how to remake the forward declarations file
include $D/common/scripts/forward_declarations




################### an additional rule for the examples

examples:
	cd examples ; $(MAKE)



################### some clean-up rules

clean: examples-clean lib-clean 

lib-clean:
	-rm -f *~ */*~ */*/*~ include/grid/forward_declarations.h*
	cd lib ; $(MAKE) clean

examples-clean:
	cd examples ; $(MAKE) clean

###################

.PHONY: all deal.II examples TAGS
.PHONY: clean lib-clean examples-clean source-clean
.PHONY: default 1d 2d 3d all clean
